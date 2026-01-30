import { useEffect, useRef, forwardRef, useImperativeHandle, useState } from 'react';
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
    getAtomsCoords: (ids: number[]) => { x: number, y: number }[];
}

interface KetcherEditorProps {
    molecule: VisualMolecule | null;
    lastSelection?: number[];
    onInit?: (ketcher: Ketcher) => void;
}

const structServiceProvider = new StandaloneStructServiceProvider();

export const KetcherEditor = forwardRef<KetcherEditorRef, KetcherEditorProps>(
    ({ molecule, lastSelection = [], onInit }, ref) => {
        const [ketcherInstance, setKetcherInstance] = useState<Ketcher | null>(null);
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
            },
            getAtomsCoords: (ids: number[]) => {
                if (!ketcherRef.current) return [];
                const coords: { x: number, y: number }[] = [];
                try {
                    // Access Ketcher's internal render objects
                    const ctab = (ketcherRef.current as any).editor.render.ctab;
                    for (const id of ids) {
                        const atom = ctab.atoms.get(id);
                        if (atom && atom.a && atom.a.pp) {
                            coords.push({ x: atom.a.pp.x, y: atom.a.pp.y });
                        }
                    }
                } catch (e) {
                    console.error("Failed to get atom coordinates:", e);
                }
                return coords;
            }
        }));

        const applySelections = (ketcher: Ketcher, mol: VisualMolecule) => {
            // Give Ketcher a moment to finish rendering the molecule
            setTimeout(() => {
                try {
                    // Prioritize atoms explicitly marked in the molecule data (Added/Selected).
                    // If the molecule has added or selected atoms (e.g. from an AI edit), 
                    // we only highlight those to focus on the changes.
                    // We only fallback to lastSelection if no structural highlights are provided.
                    const selectedAtomsFromMolecule = mol.atoms
                        .filter(a => a.ui_state === "SELECTED")
                        .map(a => a.id);

                    const addedAtoms = mol.atoms
                        .filter(a => a.diff_state === "ADDED")
                        .map(a => a.id);

                    let finalAtomSelection = [...new Set([
                        ...selectedAtomsFromMolecule,
                        ...addedAtoms
                    ])];

                    if (finalAtomSelection.length === 0 && lastSelection.length > 0) {
                        finalAtomSelection = [...lastSelection];
                    }

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
            setKetcherInstance(ketcher);
            window.ketcher = ketcher; // For debugging

            // Branding and UI customization for Ketcher
            try {
                ketcher.editor.setOptions(JSON.stringify({
                    showAtomIds: false
                }));

                // MONKEY-PATCH: Fix atom selection alignment bug within the chemist-ai project.
                // We delay the patch until the first molecule is loaded to safely access the ReAtom prototype.
                const originalSetMolecule = ketcher.setMolecule.bind(ketcher);
                ketcher.setMolecule = async (molString: string) => {
                    const result = await originalSetMolecule(molString);

                    try {
                        const ctab = ketcher.editor.render.ctab;
                        if (ctab && ctab.atoms && ctab.atoms.size > 0) {
                            const firstAtom = ctab.atoms.values().next().value;
                            const proto = Object.getPrototypeOf(firstAtom);

                            if (proto && !proto._patchedForApotheca) {
                                proto._patchedForApotheca = true;
                                proto.getLabeledSelectionContour = function (render: any, highlightPadding = 0) {
                                    const { paper, options } = render;
                                    const { fontszInPx, radiusScaleFactor, microModeScale } = options;
                                    const padding = fontszInPx * radiusScaleFactor + highlightPadding;
                                    const radius = fontszInPx * radiusScaleFactor * 2 + highlightPadding;

                                    // FIX: Use the atom's actual canvas position as the anchor
                                    const ps = this.a.pp.scaled(microModeScale);

                                    const box = this.getVBoxObj(render)!;
                                    const ps1 = box.p0.scaled(microModeScale);
                                    const ps2 = box.p1.scaled(microModeScale);

                                    const width = Math.max(ps2.x - ps1.x, fontszInPx * 0.8);
                                    const height = fontszInPx * 1.25;

                                    return paper.rect(
                                        ps.x - width / 2 - padding,
                                        ps.y - height / 2 - padding,
                                        width + padding * 2,
                                        height + padding * 2,
                                        radius,
                                    );
                                };
                            }
                        }
                    } catch (err) {
                        console.warn("Apotheca: Failed to apply rendering patch", err);
                    }
                    return result;
                };
            } catch (e) {
                console.warn("Failed to set Ketcher options or prepare patch:", e);
            }

            if (onInit) {
                onInit(ketcher);
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
            if (ketcherInstance && molecule?.mol_block) {
                if (molecule.molecule_id !== lastMoleculeId.current) {
                    lastMoleculeId.current = molecule.molecule_id;
                    ketcherInstance.setMolecule(molecule.mol_block).then(() => {
                        applySelections(ketcherInstance, molecule);
                    });
                }
            }
        }, [molecule, ketcherInstance]);

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
