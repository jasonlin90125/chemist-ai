import { useState } from 'react';
import { VisualMolecule } from '../types/molecule';
import { moleculeApi } from '../api/client';

export const useMolecule = () => {
    const [molecule, setMolecule] = useState<VisualMolecule | null>(null);
    const [originalMolfile, setOriginalMolfile] = useState<string | null>(null);

    const [status, setStatus] = useState<"IDLE" | "LOADING" | "DIFFING">("IDLE");
    const [selectedIndices, setSelectedIndices] = useState<number[]>([]);
    const [error, setError] = useState<string | null>(null);
    const [variantIdx, setVariantIdx] = useState(0);
    const [lastRequest, setLastRequest] = useState<{ action: string, parameters: any, indices: number[] } | null>(null);

    // Load initial molecule
    const loadInitial = async () => {
        try {
            const data = await moleculeApi.getInitial();
            setMolecule(data);
        } catch (e) {
            console.error("Failed to load molecule", e);
            setError("Failed to load initial molecule.");
        }
    };

    const updateSelectionState = (indices: number[]) => {
        setSelectedIndices(indices);
    };

    const requestEdit = async (prompt: string, manualIndices?: number[], liveMolfile?: string | null) => {
        setStatus("LOADING");
        setError(null);
        setVariantIdx(0); // Reset variant on new prompt

        const finalIndices = manualIndices || selectedIndices;
        const baseMolfile = liveMolfile || (molecule?.mol_block) || null;
        setOriginalMolfile(baseMolfile);

        try {
            const currentMolecule: VisualMolecule = {
                molecule_id: molecule?.molecule_id || "temp",
                atoms: [],
                bonds: [],
                mol_block: baseMolfile || ""
            };

            // SIMPLE ROUTING: Detect common commands
            const p = prompt.toLowerCase().trim();
            let proposal: VisualMolecule | null = null;

            const addMatch = p.match(/^add\s+([\w#*()=@-]+)(\s+.*)?$/);

            if (addMatch && finalIndices.length === 1) {
                const fragment = addMatch[1];
                const action = "add_substructure";
                const params = { fragment, variant_idx: 0 };
                setLastRequest({ action, parameters: params, indices: finalIndices });

                proposal = await moleculeApi.simpleEdit({
                    action,
                    current_molecule: currentMolecule,
                    selected_indices: finalIndices,
                    parameters: params
                });
            } else if (p === "remove" || p === "delete") {
                const action = "remove_atoms";
                const params = {};
                setLastRequest({ action, parameters: params, indices: finalIndices });

                proposal = await moleculeApi.simpleEdit({
                    action,
                    current_molecule: currentMolecule,
                    selected_indices: finalIndices,
                    parameters: params
                });
            } else if (p === "aromatize" || p === "dearomatize") {
                const action = p;
                const params = {};
                setLastRequest({ action, parameters: params, indices: finalIndices });
                proposal = await moleculeApi.simpleEdit({
                    action: p,
                    current_molecule: currentMolecule,
                    selected_indices: finalIndices,
                    parameters: params
                });
            }

            if (!proposal) {
                // If not simple, we could still track the LLM variant cycling in theory,
                // but for now let's focus on tools.
                setLastRequest(null);
                proposal = await moleculeApi.edit({
                    current_molecule: currentMolecule,
                    user_prompt: prompt,
                    selected_indices: finalIndices
                });
            }

            setMolecule(proposal);
            setStatus("DIFFING");
        } catch (e: any) {
            console.error(e);
            setError(e.response?.data?.detail || e.message || "Error");
            setStatus("IDLE");
        }
    };

    const cycleVariant = async (direction: number) => {
        if (!lastRequest || !originalMolfile) return;

        const newIdx = Math.max(0, variantIdx + direction);
        setVariantIdx(newIdx);
        setStatus("LOADING");

        try {
            const originalMolecule: VisualMolecule = {
                molecule_id: molecule?.molecule_id || "temp",
                atoms: [],
                bonds: [],
                mol_block: originalMolfile
            };

            const proposal = await moleculeApi.simpleEdit({
                action: lastRequest.action,
                current_molecule: originalMolecule,
                selected_indices: lastRequest.indices,
                parameters: { ...lastRequest.parameters, variant_idx: newIdx }
            });

            setMolecule(proposal);
            setStatus("DIFFING");
        } catch (e: any) {
            console.error(e);
            const msg = e.response?.data?.detail || e.message || "Unique variant cycling failed.";
            setError(msg);
            setStatus("DIFFING");
        }
    };

    const acceptChange = (): string | null => {
        if (!molecule) return null;
        const molBlock = molecule.mol_block || null;

        const cleanMolecule: VisualMolecule = {
            ...molecule,
            atoms: molecule.atoms.map(a => ({ ...a, diff_state: "EXISTING" })),
            bonds: molecule.bonds.map(b => ({ ...b, diff_state: "EXISTING" }))
        };

        setMolecule(cleanMolecule);
        setOriginalMolfile(null);
        setSelectedIndices([]);
        setLastRequest(null);
        setStatus("IDLE");
        return molBlock;
    };

    const rejectChange = (): string | null => {
        const revertTo = originalMolfile || null;
        if (revertTo && molecule) {
            setMolecule({
                ...molecule,
                mol_block: revertTo,
                atoms: [],
                bonds: []
            });
        }
        setOriginalMolfile(null);
        setLastRequest(null);
        setStatus("IDLE");
        return revertTo;
    };

    return {
        molecule,
        status,
        error,
        loadInitial,
        requestEdit,
        cycleVariant,
        acceptChange,
        rejectChange,
        handleSelection: updateSelectionState
    };
};
