import { useEffect, useRef, forwardRef, useImperativeHandle } from 'react';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import 'ketcher-react/dist/index.css';
import { VisualMolecule } from '../../types/molecule';

// Define the Ketcher interface
interface Ketcher {
    setMolecule: (mol: string) => Promise<void>;
    getMolfile: (format?: string) => Promise<string>;
    getSmiles: () => Promise<string>;
    layout: () => void;
    editor: {
        selection: (selection?: {
            atoms?: number[];
            bonds?: number[];
        } | null) => {
            atoms?: number[];
            bonds?: number[];
        } | null;
    };
}

export interface KetcherEditorRef {
    getSelectedAtoms: () => number[];
    getMolfile: () => Promise<string | null>;
    setMolecule: (molBlock: string) => void;
    layout: () => void;
}

interface KetcherEditorProps {
    molecule: VisualMolecule | null;
    onInit?: (ketcher: Ketcher) => void;
}

const structServiceProvider = new StandaloneStructServiceProvider();

export const KetcherEditor = forwardRef<KetcherEditorRef, KetcherEditorProps>(
    ({ molecule, onInit }, ref) => {
        const ketcherRef = useRef<Ketcher | null>(null);

        useImperativeHandle(ref, () => ({
            getSelectedAtoms: () => {
                if (!ketcherRef.current) return [];
                const selection = ketcherRef.current.editor.selection();
                return selection?.atoms || [];
            },
            getMolfile: async () => {
                if (!ketcherRef.current) return null;
                try {
                    // Force v3000 for better compatibility with RDKit
                    return await ketcherRef.current.getMolfile('v3000');
                } catch (e) {
                    console.error("Failed to get molfile:", e);
                    return null;
                }
            },
            setMolecule: (molBlock: string) => {
                if (ketcherRef.current) {
                    ketcherRef.current.setMolecule(molBlock);
                }
            },
            layout: () => {
                if (ketcherRef.current) {
                    ketcherRef.current.layout();
                }
            }
        }));

        const applySelections = (ketcher: Ketcher, mol: VisualMolecule) => {
            // Give Ketcher a moment to finish rendering the molecule
            setTimeout(() => {
                try {
                    const selectedAtoms = mol.atoms
                        .filter(a => a.ui_state === "SELECTED")
                        .map(a => a.id);

                    const addedAtoms = mol.atoms
                        .filter(a => a.diff_state === "ADDED")
                        .map(a => a.id);

                    const finalAtomSelection = [...new Set([...selectedAtoms, ...addedAtoms])];

                    // Also select bonds where both source and target atoms are selected
                    const atomSet = new Set(finalAtomSelection);
                    const finalBondSelection = mol.bonds
                        .map((b, idx) => ({ bond: b, idx }))
                        .filter(({ bond }) => atomSet.has(bond.source) && atomSet.has(bond.target))
                        .map(({ idx }) => idx);

                    if (finalAtomSelection.length > 0) {
                        ketcher.editor.selection({
                            atoms: finalAtomSelection,
                            bonds: finalBondSelection
                        });
                    } else {
                        ketcher.editor.selection(null);
                    }
                } catch (e) {
                    console.warn("Failed to apply selection to Ketcher:", e);
                }
            }, 150);
        };

        const handleInit = (ketcher: any) => {
            ketcherRef.current = ketcher;
            window.ketcher = ketcher; // For debugging

            if (onInit) {
                onInit(ketcher);
            }

            // Initial load
            if (molecule?.mol_block) {
                ketcher.setMolecule(molecule.mol_block).then(() => {
                    applySelections(ketcher, molecule);
                });
            }
        };

        // React to molecule updates from AI
        useEffect(() => {
            if (ketcherRef.current && molecule?.mol_block) {
                ketcherRef.current.setMolecule(molecule.mol_block).then(() => {
                    if (ketcherRef.current) {
                        applySelections(ketcherRef.current, molecule);
                    }
                });
            }
        }, [molecule]);

        return (
            <div className="w-full h-full relative">
                <Editor
                    staticResourcesUrl={""}
                    structServiceProvider={structServiceProvider as any}
                    onInit={handleInit}
                    errorHandler={(err) => console.error(err)}
                />
            </div>
        );
    }
);

// Add window type
declare global {
    interface Window {
        ketcher: any;
    }
}
