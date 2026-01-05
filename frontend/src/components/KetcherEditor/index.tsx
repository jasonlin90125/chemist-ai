import { useEffect, useRef, forwardRef, useImperativeHandle } from 'react';
import { Editor } from 'ketcher-react';
import { StandaloneStructServiceProvider } from 'ketcher-standalone';
import 'ketcher-react/dist/index.css';
import { VisualMolecule } from '../../types/molecule';

// Define the Ketcher interface mostly for 'setMolecule'
interface Ketcher {
    setMolecule: (mol: string) => void;
    getMolfile: () => Promise<string>;
    getSmiles: () => Promise<string>;
    editor: {
        selection: () => {
            atoms?: number[];
            bonds?: number[];
        } | null;
    };
}

export interface KetcherEditorRef {
    getSelectedAtoms: () => number[];
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
            }
        }));

        const handleInit = (ketcher: any) => {
            ketcherRef.current = ketcher;
            window.ketcher = ketcher; // For debugging

            if (onInit) {
                onInit(ketcher);
            }

            // Initial load
            if (molecule?.mol_block) {
                ketcher.setMolecule(molecule.mol_block);
            }
        };

        // React to molecule updates from AI
        useEffect(() => {
            if (ketcherRef.current && molecule?.mol_block) {
                ketcherRef.current.setMolecule(molecule.mol_block);
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
