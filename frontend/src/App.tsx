import { useEffect, useRef } from 'react';
import { KetcherEditor, KetcherEditorRef } from './components/KetcherEditor';
import { ChatPanel } from './components/ChatPanel';
import { DiffActionBar } from './components/DiffActionBar';
import { useMolecule } from './hooks/useMolecule';

function App() {
    const {
        molecule,
        status,
        error,
        loadInitial,
        requestEdit,
        acceptChange,
        rejectChange
    } = useMolecule();

    const ketcherRef = useRef<KetcherEditorRef>(null);

    useEffect(() => {
        loadInitial();
    }, []);

    const handleSendPrompt = async (prompt: string) => {
        const selectedIds = ketcherRef.current?.getSelectedAtoms() || [];
        // Get live molecule from Ketcher to ensure we're editing what's displayed
        const liveMolfile = await ketcherRef.current?.getMolfile();
        requestEdit(prompt, selectedIds, liveMolfile);
    };

    return (
        <div className="flex w-full h-screen bg-white overflow-hidden">

            {/* Left Sidebar */}
            <ChatPanel
                onSendPrompt={handleSendPrompt}
                isLoading={status === 'LOADING'}
                error={error}
            />

            {/* Main Stage */}
            <div className="flex-1 relative h-full bg-gray-50">
                <KetcherEditor
                    ref={ketcherRef}
                    molecule={molecule}
                />

                {/* Diff Controls (Floating over editor) */}
                {status === 'DIFFING' && (
                    <DiffActionBar
                        onAccept={acceptChange}
                        onReject={rejectChange}
                    />
                )}

                {/* Loading Indicator (Centered over editor) */}
                {status === 'LOADING' && (
                    <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 bg-white/80 backdrop-blur shadow-xl rounded-xl px-8 py-6 flex flex-col items-center border border-gray-100">
                        <div className="w-10 h-10 border-3 border-blue-500 border-t-transparent rounded-full animate-spin mb-3"></div>
                        <span className="text-gray-600 text-sm font-medium tracking-wide">Synthesizing...</span>
                    </div>
                )}
            </div>

        </div>
    );
}

export default App;
