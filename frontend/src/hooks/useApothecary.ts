import { useState, useEffect } from 'react';
import { ApothecaryEntry } from '../types/apothecary';
import { VisualMolecule } from '../types/molecule';
import { v4 as uuidv4 } from 'uuid';
import { moleculeApi } from '../api/client';

export const useApothecary = (onShowToast?: (msg: string, type?: 'success' | 'info' | 'warning' | 'error') => void) => {
    const [entries, setEntries] = useState<ApothecaryEntry[]>([]);

    // Initialize from localStorage
    useEffect(() => {
        const saved = localStorage.getItem('chemist_apothecary');
        if (saved) {
            try {
                setEntries(JSON.parse(saved));
            } catch (e) {
                console.error("Failed to parse apothecary from localStorage", e);
            }
        }
    }, []);

    // Save to localStorage
    useEffect(() => {
        localStorage.setItem('chemist_apothecary', JSON.stringify(entries));
    }, [entries]);

    const addToApothecary = (mol: VisualMolecule) => {
        // Check for duplicates by SMILES
        const isDuplicate = entries.some(e => e.smiles === mol.smiles);
        if (isDuplicate) {
            onShowToast?.("Already added to Apothecary.", "warning");
            return false;
        }

        const newEntry: ApothecaryEntry = {
            id: `MOL-${uuidv4().substring(0, 8).toUpperCase()}`,
            smiles: mol.smiles || "Unknown SMILES",
            mol_block: mol.mol_block || "",
            svg: mol.svg,
            timestamp: Date.now()
        };
        setEntries(prev => [newEntry, ...prev]);
        onShowToast?.("Molecule added to Apothecary", "success");
        return true;
    };

    const removeFromApothecary = (id: string) => {
        setEntries(prev => prev.filter(e => e.id !== id));
    };

    const clearApothecary = () => {
        setEntries([]);
        localStorage.removeItem('chemist_apothecary');
    };

    const exportToCSV = () => {
        const headers = ["ID", "SMILES", "Timestamp"];
        const rows = entries.map(e => [
            e.id,
            e.smiles,
            new Date(e.timestamp).toISOString()
        ]);

        const csvContent = [headers, ...rows].map(r => r.join(",")).join("\n");
        const blob = new Blob([csvContent], { type: 'text/csv' });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `apothecary_${new Date().toISOString().split('T')[0]}.csv`;
        a.click();
    };

    const exportToSMI = () => {
        const smiContent = entries.map(e => `${e.smiles} ${e.id}`).join("\n");
        const blob = new Blob([smiContent], { type: 'text/plain' });
        const url = window.URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `apothecary_${new Date().toISOString().split('T')[0]}.smi`;
        a.click();
    };

    const exportToSDF = async () => {
        try {
            const blob = await moleculeApi.exportSDF(
                entries.map(e => ({
                    mol_block: e.mol_block,
                    id: e.id
                }))
            );

            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = `apothecary_${new Date().toISOString().split('T')[0]}.sdf`;
            a.click();
            window.URL.revokeObjectURL(url);
        } catch (e) {
            console.error("Failed to export to SDF", e);
        }
    };

    return {
        entries,
        addToApothecary,
        removeFromApothecary,
        clearApothecary,
        exportToCSV,
        exportToSMI,
        exportToSDF
    };
};
