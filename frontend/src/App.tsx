import { useEffect, useRef, useState } from 'react';
import { KetcherEditor, KetcherEditorRef } from './components/KetcherEditor';
import { ChatPanel } from './components/ChatPanel';
import { DiffActionBar } from './components/DiffActionBar';
import { useMolecule } from './hooks/useMolecule';
import { useClipboard } from './hooks/useClipboard';
import { ClipboardEntry } from './types/clipboard';
import { Message } from './types/chat';
import { v4 as uuidv4 } from 'uuid';

function App() {
    const {
        molecule,
        status,
        error,
        proposalsCount,
        currentProposalIdx,
        loadInitial,
        requestEdit,
        cycleVariant,
        acceptChange,
        rejectChange,
        revertTo
    } = useMolecule();

    const { entries, addToClipboard, clearClipboard } = useClipboard();
    const [activeTab, setActiveTab] = useState<'chat' | 'history'>('chat');
    const [messages, setMessages] = useState<Message[]>([]);
    const [lastSelection, setLastSelection] = useState<number[]>([]);

    const ketcherRef = useRef<KetcherEditorRef>(null);

    // Poll for selection to track last selection
    useEffect(() => {
        const interval = setInterval(() => {
            const current = ketcherRef.current?.getSelectedAtoms() || [];
            if (current.length > 0) {
                setLastSelection(current);
            }
        }, 200);

        return () => {
            clearInterval(interval);
        };
    }, []);

    useEffect(() => {
        const init = async () => {
            const data = await loadInitial();
            if (data) {
                // Log the starting material in history
                addToClipboard(data, 'accepted');
            }
        };
        init();
    }, []);

    const isFirstLoad = !molecule && status === 'LOADING';

    // Watch for initial molecule to log it specifically as "Starting Material"
    useEffect(() => {
        if (molecule && entries.length === 0 && status === 'IDLE') {
            addToClipboard(molecule, 'accepted');
        }
    }, [molecule, entries.length, status]);

    const handleSendPrompt = async (prompt: string) => {
        const userMsg: Message = {
            id: uuidv4(),
            role: 'user',
            content: prompt,
            timestamp: Date.now()
        };
        setMessages(prev => [...prev, userMsg]);

        // Get current selection from Ketcher
        let selectedIds = ketcherRef.current?.getSelectedAtoms() || [];

        if (selectedIds.length === 0 && lastSelection.length > 0) {
            selectedIds = lastSelection;
            console.log("Using tracked lastSelection:", lastSelection);
        }

        // Map indices to persistent map numbers
        const selectedMaps = selectedIds.map(idx => {
            const atom = molecule?.atoms?.find(a => a.id === idx);
            return atom?.atom_map;
        }).filter(m => m !== undefined && m !== null) as number[];

        console.log("Selection IDs:", selectedIds);
        console.log("Selection Maps:", selectedMaps);

        if (selectedMaps.length === 0 && selectedIds.length > 0) {
            console.log("WARNING: No atom maps found for selection. Backend will use raw indices.");
        }

        let liveMolfile = null;
        try {
            liveMolfile = await ketcherRef.current?.getMolfile();
            // liveSmiles is fetched but not needed for the request currently
            await ketcherRef.current?.getSmiles();
        } catch (e) {
            console.error("Error fetching live molecule data:", e);
        }

        const selectedCoords = ketcherRef.current?.getAtomsCoords(selectedIds) || [];
        console.log("Selection Coords:", selectedCoords);

        const result = await requestEdit(prompt, selectedIds, liveMolfile, selectedMaps, selectedCoords);

        // Clear selection after use
        setLastSelection([]);

        if (result) {
            const agentMsg: Message = {
                id: uuidv4(),
                role: 'assistant',
                content: `Generated ${result.proposalCount} variants for: ${result.action}`,
                timestamp: Date.now(),
                type: 'edit',
                metadata: {
                    action: result.action,
                    proposalCount: result.proposalCount
                }
            };
            setMessages(prev => [...prev, agentMsg]);
        } else if (error) {
            const errorMsg: Message = {
                id: uuidv4(),
                role: 'system',
                content: error,
                timestamp: Date.now(),
                type: 'error'
            };
            setMessages(prev => [...prev, errorMsg]);
        }
    };

    const handleAccept = () => {
        const molBlock = acceptChange();
        if (molBlock && ketcherRef.current) {
            if (molecule) {
                addToClipboard(molecule, 'accepted');
            }
            ketcherRef.current.setMolecule(molBlock);
            // Note: layout() call removed to preserve molecule orientation
        }
    };

    const handleReject = () => {
        const currentProposal = molecule;
        const revertTo = rejectChange();
        if (revertTo && ketcherRef.current) {
            if (currentProposal) {
                addToClipboard(currentProposal, 'rejected');
            }
            ketcherRef.current.setMolecule(revertTo);
        }
    };

    const handleRevert = (entry: ClipboardEntry) => {
        revertTo(entry);
        if (entry.mol_block && ketcherRef.current) {
            ketcherRef.current.setMolecule(entry.mol_block);
        }
    };

    return (
        <div className="flex w-full h-screen bg-white overflow-hidden font-sans">

            {/* Main Stage (Flex-1) */}
            <div className="flex-1 relative h-full bg-gray-50 flex flex-col">
                <div className="flex-1 relative">
                    {/* Full page loader for first visit */}
                    {isFirstLoad && (
                        <div className="absolute inset-0 bg-white z-[100] flex flex-col items-center justify-center p-8 transition-opacity duration-500">
                            <div className="w-16 h-16 border-4 border-chemist-primary border-t-transparent rounded-full animate-spin mb-6"></div>
                            <h2 className="text-2xl font-bold text-gray-800 mb-2">Chemist AI</h2>
                            <p className="text-gray-500 font-medium">Reconstructing lead compounds...</p>
                        </div>
                    )}

                    <KetcherEditor
                        ref={ketcherRef}
                        molecule={molecule}
                        lastSelection={lastSelection}
                        onInit={() => console.log("Ketcher Ready")}
                    />

                    {/* Diff Controls (Floating over editor) */}
                    {status === 'DIFFING' && (
                        <DiffActionBar
                            onAccept={handleAccept}
                            onReject={handleReject}
                            onCycle={cycleVariant}
                            proposalsCount={proposalsCount}
                            currentProposalIdx={currentProposalIdx}
                        />
                    )}

                    {/* Loading Indicator (Centered over editor) */}
                    {status === 'LOADING' && (
                        <div className="absolute top-1/2 left-1/2 -translate-x-1/2 -translate-y-1/2 bg-white/80 backdrop-blur shadow-xl rounded-xl px-8 py-6 flex flex-col items-center border border-gray-100 ring-1 ring-black/5 z-50">
                            <div className="w-10 h-10 border-3 border-chemist-primary border-t-transparent rounded-full animate-spin mb-3"></div>
                            <span className="text-gray-600 text-sm font-bold tracking-tight">Synthesizing...</span>
                        </div>
                    )}
                </div>
            </div>

            {/* Right Sidebar - Agent/Chat & History */}
            <ChatPanel
                onSendPrompt={handleSendPrompt}
                isLoading={status === 'LOADING'}
                error={error}
                activeTab={activeTab}
                onTabChange={setActiveTab}
                messages={messages}
                entries={entries}
                onRevert={handleRevert}
                onClear={clearClipboard}
            />

        </div>
    );
}

export default App;
