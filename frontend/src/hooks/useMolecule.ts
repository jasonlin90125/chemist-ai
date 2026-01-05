import { useState } from 'react';
import { VisualMolecule } from '../types/molecule';
import { moleculeApi } from '../api/client';

export const useMolecule = () => {
    const [molecule, setMolecule] = useState<VisualMolecule | null>(null);
    const [originalMolecule, setOriginalMolecule] = useState<VisualMolecule | null>(null);

    const [status, setStatus] = useState<"IDLE" | "LOADING" | "DIFFING">("IDLE");
    const [selectedIndices, setSelectedIndices] = useState<number[]>([]);

    // Load initial molecule
    const loadInitial = async () => {
        try {
            const data = await moleculeApi.getInitial();
            setMolecule(data);
        } catch (e) {
            console.error("Failed to load molecule", e);
        }
    };

    const updateSelectionState = (indices: number[]) => {
        setSelectedIndices(indices);
        if (molecule) {
            setMolecule({
                ...molecule,
                atoms: molecule.atoms.map(a => ({
                    ...a,
                    ui_state: indices.includes(a.id) ? "SELECTED" : "DEFAULT"
                }))
            });
        }
    };


    const requestEdit = async (prompt: string, manualIndices?: number[]) => {
        if (!molecule) return;

        setStatus("LOADING");
        // Save current as original for "Reject" capability
        setOriginalMolecule(molecule);

        try {
            const proposal = await moleculeApi.edit({
                current_molecule: molecule,
                user_prompt: prompt,
                selected_indices: manualIndices || selectedIndices
            });

            setMolecule(proposal);
            setStatus("DIFFING");
        } catch (e) {
            console.error(e);
            setStatus("IDLE");
        }
    };

    const acceptChange = () => {
        if (!molecule) return;
        // Commit: Clear diff states
        const cleanMolecule: VisualMolecule = {
            ...molecule,
            atoms: molecule.atoms.map(a => ({ ...a, diff_state: "EXISTING" })),
            bonds: molecule.bonds.map(b => ({ ...b, diff_state: "EXISTING" }))
        };

        setMolecule(cleanMolecule);
        setOriginalMolecule(null);
        setSelectedIndices([]);
        setStatus("IDLE");
    };

    const rejectChange = () => {
        if (originalMolecule) {
            setMolecule(originalMolecule);
        }
        setOriginalMolecule(null);
        setStatus("IDLE");
    };

    return {
        molecule,
        status,
        loadInitial,
        requestEdit,
        acceptChange,
        rejectChange,
        handleSelection: updateSelectionState
    };
};
