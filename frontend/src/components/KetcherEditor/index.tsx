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
        setOptions: (opts: string) => void;
    };
}

export interface KetcherEditorRef {
    getSelectedAtoms: () => number[];
    getMolfile: () => Promise<string | null>;
    getSmiles: () => Promise<string | null>;
    setMolecule: (molBlock: string) => void;
    layout: () => void;
}

interface KetcherEditorProps {
    molecule: VisualMolecule | null;
    lastSelection?: number[];
    onInit?: (ketcher: Ketcher) => void;
}

const structServiceProvider = new StandaloneStructServiceProvider();

export const KetcherEditor = forwardRef<KetcherEditorRef, KetcherEditorProps>(
    ({ molecule, lastSelection = [], onInit }, ref) => {
        const ketcherRef = useRef<Ketcher | null>(null);
        const lastMoleculeId = useRef<string | null>(null);

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
            getSmiles: async () => {
                if (!ketcherRef.current) return null;
                try {
                    return await ketcherRef.current.getSmiles();
                } catch (e) {
                    console.error("Failed to get SMILES:", e);
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
                    const selectedAtomsFromMolecule = mol.atoms
                        .filter(a => a.ui_state === "SELECTED")
                        .map(a => a.id);

                    const addedAtoms = mol.atoms
                        .filter(a => a.diff_state === "ADDED")
                        .map(a => a.id);

                    // Merge AI-added atoms with the user's last selection to preserve context
                    const finalAtomSelection = [...new Set([
                        ...selectedAtomsFromMolecule,
                        ...addedAtoms,
                        ...lastSelection
                    ])];

                    // Filter selection to only include atoms that actually exist in current structure
                    const existingAtomIds = new Set(mol.atoms.map(a => a.id));
                    const validSelection = finalAtomSelection.filter(id => existingAtomIds.has(id));

                    // Also select bonds where both source and target atoms are selected
                    const atomSet = new Set(validSelection);
                    const finalBondSelection = mol.bonds
                        .map((b, idx) => ({ bond: b, idx }))
                        .filter(({ bond }) => atomSet.has(bond.source) && atomSet.has(bond.target))
                        .map(({ idx }) => idx);

                    if (validSelection.length > 0) {
                        ketcher.editor.selection({
                            atoms: validSelection,
                            bonds: finalBondSelection
                        });
                    } else {
                        ketcher.editor.selection(null);
                    }
                } catch (e) {
                    console.warn("Failed to apply selection to Ketcher:", e);
                }
            }, 100);
        };

        const handleInit = (ketcher: any) => {
            ketcherRef.current = ketcher;
            window.ketcher = ketcher; // For debugging

            // Hide atom IDs from the user - they're for internal use only
            try {
                ketcher.editor.setOptions(JSON.stringify({
                    showAtomIds: false
                }));
            } catch (e) {
                console.warn("Failed to set Ketcher options:", e);
            }

            if (onInit) {
                onInit(ketcher);
            }

            // Initial load - add small delay to ensure Ketcher engine is ready
            if (molecule?.mol_block) {
                setTimeout(() => {
                    try {
                        ketcher.setMolecule(molecule.mol_block).then(() => {
                            applySelections(ketcher, molecule);
                        });
                    } catch (e) {
                        console.error("Initial setMolecule failed:", e);
                    }
                }, 200);
            }
        };

        // Monitor Ketcher for atom map numbers and hide them
        useEffect(() => {
            // Add a style tag to hide map numbers via CSS as a first layer
            const style = document.createElement('style');
            style.innerHTML = `
                /* We can't target text content via CSS easily, so we use DOM approach below,
                   but we can hide the indices layer natively via CSS as well. */
                .Ketcher-root svg g[id^="indices"] { display: none !important; }
            `;
            document.head.appendChild(style);

            const hideMaps = () => {
                const editorEl = document.querySelector('.Ketcher-root') || document.body;

                // Ketcher renders map numbers as .N. (e.g. .33.)
                const texts = editorEl.querySelectorAll('text');
                texts.forEach(t => {
                    const content = t.textContent?.trim();
                    if (content && /^\.\d+\.$/.test(content)) {
                        t.setAttribute('style', 'display: none !important; visibility: hidden !important;');
                    }
                });
            };

            // Run frequently to ensure they stay hidden during interactions
            const interval = setInterval(hideMaps, 200);

            // Also observe mutations for smoother hiding
            const observer = new MutationObserver(hideMaps);
            observer.observe(document.body, {
                childList: true,
                subtree: true,
                characterData: true
            });

            return () => {
                clearInterval(interval);
                observer.disconnect();
                if (style.parentNode) {
                    style.parentNode.removeChild(style);
                }
            };
        }, [molecule?.molecule_id]); // Re-run when molecule changes

        // React to molecule updates from AI
        useEffect(() => {
            if (ketcherRef.current && molecule?.mol_block) {
                if (molecule.molecule_id !== lastMoleculeId.current) {
                    lastMoleculeId.current = molecule.molecule_id;
                    ketcherRef.current.setMolecule(molecule.mol_block).then(() => {
                        if (ketcherRef.current) {
                            applySelections(ketcherRef.current, molecule);
                        }
                    });
                }
            }
        }, [molecule, molecule?.mol_block]);

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
