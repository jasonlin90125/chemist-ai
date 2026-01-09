import { useState } from 'react';
import { VisualMolecule } from '../types/molecule';
import { moleculeApi } from '../api/client';
import { ClipboardEntry } from '../types/clipboard';
import { v4 as uuidv4 } from 'uuid';

export const useMolecule = () => {
    const [molecule, setMolecule] = useState<VisualMolecule | null>(null);
    const [originalMolfile, setOriginalMolfile] = useState<string | null>(null);

    const [status, setStatus] = useState<"IDLE" | "LOADING" | "DIFFING">("IDLE");
    const [selectedIndices, setSelectedIndices] = useState<number[]>([]);
    const [error, setError] = useState<string | null>(null);

    // Multi-proposal state
    const [proposals, setProposals] = useState<VisualMolecule[]>([]);
    const [currentProposalIdx, setCurrentProposalIdx] = useState(0);

    const [lastRequest, setLastRequest] = useState<{ action: string, parameters: any, indices: number[] } | null>(null);

    // Load initial molecule
    const loadInitial = async () => {
        try {
            const data = await moleculeApi.getInitial();
            setMolecule(data);
            return data;
        } catch (e) {
            console.error("Failed to load molecule", e);
            setError("Failed to load initial molecule.");
            return null;
        }
    };

    const updateSelectionState = (indices: number[]) => {
        setSelectedIndices(indices);
    };

    const requestEdit = async (prompt: string, manualIndices?: number[], liveMolfile?: string | null, selectedMaps: number[] = []) => {
        setStatus("LOADING");
        setError(null);
        setProposals([]);
        setCurrentProposalIdx(0);

        const finalIndices = manualIndices || selectedIndices;
        const baseMolfile = liveMolfile || (molecule?.mol_block) || null;
        setOriginalMolfile(baseMolfile);

        try {
            const molBlockString = typeof baseMolfile === 'string' ? baseMolfile : (baseMolfile ? String(baseMolfile) : "");

            const currentMolecule: VisualMolecule = {
                molecule_id: String(molecule?.molecule_id || "temp"),
                atoms: [],
                bonds: [],
                mol_block: molBlockString
            };

            const p = prompt.toLowerCase().trim();
            let results: VisualMolecule[] = [];
            let actionName = "AI Synthesis";

            const addMatch = p.match(/^add\s+([\w#*()=@-]+)(\s+.*)?$/);

            if (addMatch && finalIndices.length === 1) {
                const fragment = addMatch[1];
                actionName = `Add ${fragment}`;
                results = await moleculeApi.multiEdit({
                    action: "add_substructure",
                    current_molecule: currentMolecule,
                    selected_indices: [...finalIndices],
                    selected_maps: [...selectedMaps],
                    parameters: { fragment }
                });
            } else if (p === "remove" || p === "delete") {
                actionName = "Remove Atoms";
                const single = await moleculeApi.simpleEdit({
                    action: "remove_atoms",
                    current_molecule: currentMolecule,
                    selected_indices: [...finalIndices],
                    selected_maps: [...selectedMaps],
                    parameters: {}
                });
                results = [single];
            } else if (p === "aromatize" || p === "dearomatize") {
                actionName = p.charAt(0).toUpperCase() + p.slice(1);
                const single = await moleculeApi.simpleEdit({
                    action: p,
                    current_molecule: currentMolecule,
                    selected_indices: [...finalIndices],
                    selected_maps: [...selectedMaps],
                    parameters: {}
                });
                results = [single];
            }

            if (!results || results.length === 0) {
                const single = await moleculeApi.edit({
                    current_molecule: currentMolecule,
                    user_prompt: prompt,
                    selected_indices: [...finalIndices],
                    selected_maps: [...selectedMaps]
                });
                results = [single];
            }

            setProposals(results);
            setMolecule(results[0]);
            setStatus("DIFFING");
            return { action: actionName, proposalCount: results.length };
        } catch (e: any) {
            console.error("Edit Request Error:", e);

            let errorMessage = "An error occurred";

            if (e.response?.data) {
                const data = e.response.data;
                if (data.detail) {
                    if (typeof data.detail === 'string') {
                        errorMessage = data.detail;
                    } else if (Array.isArray(data.detail)) {
                        // Handle FastAPI validation errors
                        errorMessage = data.detail.map((err: any) => `${err.loc.join('.')}: ${err.msg}`).join('; ');
                    } else {
                        errorMessage = JSON.stringify(data.detail);
                    }
                } else if (data.message) {
                    errorMessage = String(data.message);
                } else {
                    errorMessage = JSON.stringify(data);
                }
            } else if (e.message) {
                errorMessage = e.message === "Network Error"
                    ? "Unable to connect to server. Please check if the backend is running."
                    : String(e.message);
            }

            setError(errorMessage);
            setStatus("IDLE");
            return null;
        }
    };

    const cycleVariant = (direction: number) => {
        if (proposals.length <= 1) return;

        let nextIdx = currentProposalIdx + direction;
        if (nextIdx < 0) nextIdx = proposals.length - 1;
        if (nextIdx >= proposals.length) nextIdx = 0;

        setCurrentProposalIdx(nextIdx);
        setMolecule(proposals[nextIdx]);
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
        setProposals([]);
        setCurrentProposalIdx(0);
        setStatus("IDLE");
        return molBlock;
    };

    const rejectChange = (): string | null => {
        const revertToBlock = originalMolfile || null;
        if (revertToBlock && molecule) {
            setMolecule({
                ...molecule,
                mol_block: revertToBlock,
                atoms: [],
                bonds: []
            });
        }
        setOriginalMolfile(null);
        setLastRequest(null);
        setProposals([]);
        setCurrentProposalIdx(0);
        setStatus("IDLE");
        return revertToBlock;
    };

    const revertTo = (entry: ClipboardEntry) => {
        setMolecule({
            molecule_id: `rev_${uuidv4().split('-')[0]}`,
            atoms: [],
            bonds: [],
            mol_block: entry.mol_block,
            smiles: entry.smiles,
            svg: entry.svg
        } as VisualMolecule);
        setOriginalMolfile(null);
        setLastRequest(null);
        setProposals([]);
        setCurrentProposalIdx(0);
        setStatus("IDLE");
    };

    return {
        molecule,
        status,
        error,
        proposalsCount: proposals.length,
        currentProposalIdx,
        loadInitial,
        requestEdit,
        cycleVariant,
        acceptChange,
        rejectChange,
        revertTo,
        handleSelection: updateSelectionState
    };
};
