import { ClipboardEntry } from '../../types/clipboard';
import { History, Trash2, RotateCcw, CheckCircle2, XCircle } from 'lucide-react';

interface ClipboardPanelProps {
    entries: ClipboardEntry[];
    onRevert: (entry: ClipboardEntry) => void;
    onClear: () => void;
    isOpen: boolean;
    onClose: () => void;
}

export const ClipboardPanel = ({ entries, onRevert, onClear, isOpen, onClose }: ClipboardPanelProps) => {
    if (!isOpen) return null;

    return (
        <div className="w-80 h-full border-l border-white/20 bg-chemist-panel flex flex-col animate-in slide-in-from-right duration-300">
            {/* Header */}
            <div className="p-4 border-b border-white/10 flex items-center justify-between bg-white/5">
                <div className="flex items-center gap-2">
                    <History className="w-5 h-5 text-chemist-primary" />
                    <h2 className="text-white font-semibold text-sm uppercase tracking-wider">Checkpoint History</h2>
                </div>
                <button
                    onClick={onClose}
                    className="text-white/40 hover:text-white transition-colors text-xs"
                >
                    Close
                </button>
            </div>

            {/* List */}
            <div className="flex-1 overflow-y-auto p-3 space-y-3 custom-scrollbar">
                {entries.length === 0 ? (
                    <div className="h-full flex flex-col items-center justify-center text-white/20 py-10">
                        <History className="w-12 h-12 mb-3 opacity-20" />
                        <p className="text-sm">No checkpoints yet</p>
                    </div>
                ) : (
                    entries.map((entry) => (
                        <div
                            key={entry.id}
                            className="bg-white/5 border border-white/10 rounded-xl p-4 hover:border-white/30 transition-all group"
                        >
                            <div className="flex items-start justify-between mb-2">
                                <div className="flex items-center gap-2">
                                    {entry.status === 'accepted' ? (
                                        <CheckCircle2 className="w-4 h-4 text-chemist-success" />
                                    ) : (
                                        <XCircle className="w-4 h-4 text-chemist-danger" />
                                    )}
                                    <span className="text-white/80 text-[10px] uppercase font-bold tracking-tighter">
                                        {new Date(entry.timestamp).toLocaleTimeString()}
                                    </span>
                                </div>
                                <button
                                    onClick={() => onRevert(entry)}
                                    className="opacity-0 group-hover:opacity-100 p-1.5 rounded-lg bg-chemist-primary/20 text-chemist-primary hover:bg-chemist-primary/30 transition-all"
                                    title="Revert to this state"
                                >
                                    <RotateCcw className="w-3 h-3" />
                                </button>
                            </div>

                            <div className="bg-black/20 rounded-lg p-2 mb-2 font-mono text-[10px] text-white/40 truncate select-all">
                                {entry.smiles}
                            </div>

                            <div className="text-[10px] text-white/30 italic">
                                {entry.status === 'accepted' ? 'Merged structure' : 'Rejected proposal'}
                            </div>
                        </div>
                    ))
                )}
            </div>

            {/* Footer */}
            {entries.length > 0 && (
                <div className="p-4 border-t border-white/10">
                    <button
                        onClick={onClear}
                        className="w-full flex items-center justify-center gap-2 py-2 rounded-xl text-white/60 hover:text-white hover:bg-white/5 transition-all text-sm"
                    >
                        <Trash2 className="w-4 h-4" />
                        Clear History
                    </button>
                </div>
            )}
        </div>
    );
};
