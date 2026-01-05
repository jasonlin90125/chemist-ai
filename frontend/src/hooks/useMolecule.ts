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

    const handleSelection = (indices: number[]) => {
        // If we are in diffing mode, block selection
        if (status === 'DIFFING') return;

        // Toggle logic
        setSelectedIndices(prev => {
            // For simple single select or toggle:
            // If clicked is in list, remove it. Else add it.
            const clicked = indices[0]; // Assuming single click for now from canvas
            if (prev.includes(clicked)) {
                return prev.filter(i => i !== clicked);
            }
            return [...prev, clicked];
        });

        // Update UI state in molecule object for visual feedback
        if (molecule) {
            const newAtoms = molecule.atoms.map(a => ({
                ...a,
                ui_state: (selectedIndices.includes(a.id) || indices.includes(a.id)) && !(selectedIndices.includes(a.id) && indices.includes(a.id))
                    ? "SELECTED"
                    : "DEFAULT"
                // Logic above is flawed for toggle, let's simplify:
                // We'll re-render based on the new 'selectedIndices' inside the component or here.
                // Better: keep molecule pure data, let canvas derive selection visuals? 
                // Current implementation expects 'ui_state' in data. So we must map it.
            } as const));

            // We need to wait for the state update to settle to map correctly, 
            // but for now let's just cheat and assume we are adding.
            // A better way is to derived the visual molecule from base + selection state.
            // For MVP, let's just trust the backend resets states or we handle it on prompt send.

            setMolecule({ ...molecule, atoms: newAtoms });
        }
    };

    // Correction: The handleSelection logic above is messy. 
    // Let's rely on re-mapping the molecule whenever selectedIndices changes.
    // We'll do this in a useEffect inside the hook or just return a derived molecule.
    // For simplicity, let's just keep track of indices and let the UI render them as selected.
    // But our AtomNode expects `data.ui_state`. 
    // Let's update `molecule` state when `selectedIndices` changes.

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


    const requestEdit = async (prompt: string) => {
        if (!molecule) return;

        setStatus("LOADING");
        // Save current as original for "Reject" capability
        setOriginalMolecule(molecule);

        try {
            const proposal = await moleculeApi.edit({
                current_molecule: molecule,
                user_prompt: prompt,
                selected_indices: selectedIndices
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
