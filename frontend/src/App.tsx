import { useEffect } from 'react';
import { MoleculeCanvas } from './components/MoleculeCanvas';
import { ChatPanel } from './components/ChatPanel';
import { DiffActionBar } from './components/DiffActionBar';
import { useMolecule } from './hooks/useMolecule';

function App() {
    const {
        molecule,
        status,
        loadInitial,
        requestEdit,
        acceptChange,
        rejectChange,
        handleSelection
    } = useMolecule();

    useEffect(() => {
        loadInitial();
    }, []);

    return (
        <div className="w-full h-screen bg-chemist-bg relative overflow-hidden">

            {/* 3D/Graph Stage */}
            <MoleculeCanvas
                molecule={molecule}
                onSelectionChange={handleSelection}
            />

            {/* UI Overlay */}
            <ChatPanel
                onSendPrompt={requestEdit}
                isLoading={status === 'LOADING'}
            />

            {/* Diff Controls */}
            {status === 'DIFFING' && (
                <DiffActionBar
                    onAccept={acceptChange}
                    onReject={rejectChange}
                />
            )}

            {/* Loading Indicator */}
            {status === 'LOADING' && (
                <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 bg-black/50 backdrop-blur rounded-lg px-6 py-4 flex flex-col items-center">
                    <div className="w-8 h-8 border-2 border-chemist-accent border-t-transparent rounded-full animate-spin mb-2"></div>
                    <span className="text-white text-sm font-medium tracking-wide">Synthesizing...</span>
                </div>
            )}

        </div>
    );
}

export default App;
