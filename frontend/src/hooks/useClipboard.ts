import { useState, useEffect } from 'react';
import { ClipboardEntry } from '../types/clipboard';
import { VisualMolecule } from '../types/molecule';
import { v4 as uuidv4 } from 'uuid';

export const useClipboard = () => {
    const [entries, setEntries] = useState<ClipboardEntry[]>([]);

    // Initialize from localStorage
    useEffect(() => {
        const saved = localStorage.getItem('chemist_clipboard');
        if (saved) {
            try {
                setEntries(JSON.parse(saved));
            } catch (e) {
                console.error("Failed to parse clipboard from localStorage", e);
            }
        }
    }, []);

    // Save to localStorage
    useEffect(() => {
        localStorage.setItem('chemist_clipboard', JSON.stringify(entries));
    }, [entries]);

    const addToClipboard = (mol: VisualMolecule, status: 'accepted' | 'rejected' | 'pending' = 'accepted') => {
        // Extract SMILES from mol_block if possible, or just use placeholder
        // In a real app we might want the backend to return SMILES for every proposal
        const newEntry: ClipboardEntry = {
            id: uuidv4(),
            smiles: mol.smiles || "Unknown SMILES",
            mol_block: mol.mol_block || "",
            svg: mol.svg,
            status,
            timestamp: Date.now()
        };
        setEntries(prev => [newEntry, ...prev]);
    };

    const clearClipboard = () => {
        setEntries([]);
        localStorage.removeItem('chemist_clipboard');
    };

    return {
        entries,
        addToClipboard,
        clearClipboard
    };
};
