import { useRef, useState } from 'react';
import { Send, Sparkles, AlertCircle, MessageSquare, History, Trash2, RotateCcw, CheckCircle2, XCircle, Bot, User, Zap } from 'lucide-react';
import { ClipboardEntry } from '../../types/clipboard';
import { Message } from '../../types/chat';
import clsx from 'clsx';
import { Virtuoso, VirtuosoHandle } from 'react-virtuoso';

interface ChatPanelProps {
    onSendPrompt: (prompt: string) => void;
    isLoading: boolean;
    error?: string | null;
    activeTab: 'chat' | 'history';
    onTabChange: (tab: 'chat' | 'history') => void;
    messages: Message[];
}

export const ChatPanel = ({
    onSendPrompt,
    isLoading,
    error,
    activeTab,
    onTabChange,
    messages,
    entries = [],
    onRevert,
    onClear
}: ChatPanelProps & {
    entries?: ClipboardEntry[],
    onRevert?: (entry: ClipboardEntry) => void,
    onClear?: () => void
}) => {
    const [input, setInput] = useState("");
    const virtuosoRef = useRef<VirtuosoHandle>(null);

    const handleSubmit = (e: React.FormEvent) => {
        e.preventDefault();
        if (input.trim() && !isLoading) {
            onSendPrompt(input);
            setInput("");
        }
    };

    return (
        <div className="w-[400px] h-full flex flex-col border-l border-gray-200 bg-white z-20 shadow-xl overflow-hidden">

            {/* Header & Tabs */}
            <div className="border-b border-gray-100 bg-gray-50/50">
                <div className="p-4 flex items-center gap-2">
                    <Sparkles className="text-blue-600 w-5 h-5" />
                    <h2 className="text-gray-800 font-semibold tracking-tight">Chemist.ai</h2>
                </div>
                <div className="flex px-4 gap-4">
                    <button
                        onClick={() => onTabChange('chat')}
                        className={clsx(
                            "pb-2 text-sm font-medium transition-all border-b-2 flex items-center gap-2",
                            activeTab === 'chat' ? "border-blue-500 text-blue-600" : "border-transparent text-gray-400 hover:text-gray-600"
                        )}
                    >
                        <MessageSquare className="w-4 h-4" />
                        Chat
                    </button>
                    <button
                        onClick={() => onTabChange('history')}
                        className={clsx(
                            "pb-2 text-sm font-medium transition-all border-b-2 flex items-center gap-2",
                            activeTab === 'history' ? "border-blue-500 text-blue-600" : "border-transparent text-gray-400 hover:text-gray-600"
                        )}
                    >
                        <History className="w-4 h-4" />
                        History
                    </button>
                </div>
            </div>

            {/* Content Area */}
            <div className="flex-1 overflow-hidden flex flex-col">
                {activeTab === 'chat' ? (
                    <div className="flex-1 flex flex-col overflow-hidden">
                        <div className="flex-1 bg-gray-50/30 overflow-hidden">
                            <Virtuoso
                                ref={virtuosoRef}
                                style={{ height: '100%' }}
                                className="custom-scrollbar"
                                data={messages}
                                initialTopMostItemIndex={messages.length - 1}
                                followOutput="smooth"
                                itemContent={(_, msg) => (
                                    <div className="px-4 py-3">
                                        <div
                                            key={msg.id}
                                            className={clsx(
                                                "flex flex-col gap-2 animate-in fade-in slide-in-from-bottom-2 duration-300",
                                                msg.role === 'user' ? "items-end" : "items-start"
                                            )}
                                        >
                                            <div className="flex items-center gap-2 text-[10px] font-bold text-gray-400 uppercase tracking-widest px-1">
                                                {msg.role === 'user' ? (
                                                    <><span>You</span><User className="w-3 h-3" /></>
                                                ) : (
                                                    <><Bot className="w-3 h-3" /><span>Chemist.ai</span></>
                                                )}
                                            </div>

                                            {msg.type === 'edit' ? (
                                                <div className="w-full flex items-center gap-3 bg-white border border-gray-100 p-4 rounded-2xl shadow-sm border-l-4 border-l-blue-500 animate-in zoom-in-95 duration-300">
                                                    <div className="bg-blue-50 p-2.5 rounded-xl">
                                                        <Zap className="w-5 h-5 text-blue-600" />
                                                    </div>
                                                    <div className="flex-1">
                                                        <div className="text-sm font-bold text-gray-800 leading-tight mb-0.5">{msg.metadata?.action || "Edit Applied"}</div>
                                                        <div className="text-[11px] text-gray-500 font-medium">
                                                            Explored {msg.metadata?.proposalCount || 0} unique structural variants
                                                        </div>
                                                    </div>
                                                </div>
                                            ) : (
                                                <div
                                                    className={clsx(
                                                        "max-w-[90%] px-4 py-3 rounded-2xl text-sm leading-relaxed shadow-sm transition-all",
                                                        msg.role === 'user'
                                                            ? "bg-gradient-to-br from-blue-600 to-blue-700 text-white rounded-tr-none font-medium"
                                                            : "bg-white border border-gray-100 text-gray-800 rounded-tl-none"
                                                    )}
                                                >
                                                    {msg.content}
                                                </div>
                                            )}
                                        </div>
                                    </div>
                                )}
                                components={{
                                    Header: () => (
                                        <>
                                            {/* Top padding equivalent */}
                                            <div className="pt-4" />
                                            {messages.length === 0 && (
                                                <div className="px-4 pb-6">
                                                    <div className="bg-blue-50 border border-blue-100 p-4 rounded-xl text-sm text-gray-700 shadow-sm">
                                                        <div className="flex items-center gap-2 text-blue-900 font-bold mb-2">
                                                            <Zap className="w-4 h-4" />
                                                            <span>Ready for Design</span>
                                                        </div>
                                                        <p className="leading-relaxed">Describe a modification or select atoms to begin. I'll pre-calculate symmetry variants for you instantly.</p>
                                                        <div className="mt-3 flex flex-wrap gap-2 text-[10px]">
                                                            <span className="px-2 py-1 bg-white border border-blue-200 rounded-lg text-blue-600 font-semibold">"Add furan"</span>
                                                            <span className="px-2 py-1 bg-white border border-blue-200 rounded-lg text-blue-600 font-semibold">"Replace with pyridine"</span>
                                                        </div>
                                                    </div>
                                                </div>
                                            )}
                                        </>
                                    ),
                                    Footer: () => (
                                        <div className="px-4 pb-4">
                                            {error && (
                                                <div className="bg-red-50 border border-red-100 p-3 rounded-xl text-xs text-red-700 flex items-start gap-2 animate-in pulse duration-1000 mb-2">
                                                    <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                                                    <p className="font-medium">{error}</p>
                                                </div>
                                            )}

                                            {isLoading && (
                                                <div className="flex flex-col gap-2 items-start animate-in fade-in duration-300">
                                                    <div className="flex items-center gap-2 text-[10px] font-bold text-gray-400 uppercase tracking-widest px-1">
                                                        <Bot className="w-3 h-3" /><span>Chemist.ai</span>
                                                    </div>
                                                    <div className="bg-white border border-gray-100 px-4 py-3 rounded-2xl shadow-sm rounded-tl-none flex items-center gap-2">
                                                        <div className="flex gap-1">
                                                            <div className="w-1.5 h-1.5 bg-gray-300 rounded-full animate-bounce"></div>
                                                            <div className="w-1.5 h-1.5 bg-gray-300 rounded-full animate-bounce [animation-delay:0.2s]"></div>
                                                            <div className="w-1.5 h-1.5 bg-gray-300 rounded-full animate-bounce [animation-delay:0.4s]"></div>
                                                        </div>
                                                    </div>
                                                </div>
                                            )}
                                        </div>
                                    )
                                }}
                            />
                        </div>

                        {/* Input (Only in Chat) */}
                        <div className="p-4 border-t border-gray-200 bg-white">
                            <form onSubmit={handleSubmit} className="relative">
                                <input
                                    type="text"
                                    className="w-full bg-gray-50 text-gray-900 rounded-lg pl-4 pr-10 py-3 focus:outline-none focus:ring-2 focus:ring-blue-500/50 border border-gray-200 transition-all placeholder:text-gray-400"
                                    placeholder="Describe your edit..."
                                    value={input}
                                    onChange={(e) => setInput(e.target.value)}
                                    disabled={isLoading}
                                />
                                <button
                                    type="submit"
                                    disabled={isLoading || !input.trim()}
                                    className={clsx(
                                        "absolute right-2 top-2 p-1.5 rounded-md transition-all",
                                        isLoading ? "opacity-50 cursor-wait" : "hover:bg-blue-50 text-blue-600"
                                    )}
                                >
                                    <Send className="w-4 h-4" />
                                </button>
                            </form>
                        </div>
                    </div>
                ) : (
                    <div className="flex-1 flex flex-col overflow-hidden bg-gray-50/30">
                        <div className="flex-1 overflow-y-auto p-4 space-y-3 custom-scrollbar">
                            {entries.length === 0 ? (
                                <div className="h-full flex flex-col items-center justify-center text-gray-400 py-10">
                                    <History className="w-12 h-12 mb-3 opacity-20" />
                                    <p className="text-sm">No checkpoints yet</p>
                                </div>
                            ) : (
                                entries.map((entry) => (
                                    <div
                                        key={entry.id}
                                        onClick={() => onRevert && onRevert(entry)}
                                        className="bg-white border border-gray-200 rounded-xl p-4 hover:border-blue-300 hover:shadow-md transition-all group shadow-sm cursor-pointer active:scale-[0.98]"
                                        title="Click to revert to this state"
                                    >
                                        <div className="flex items-start justify-between mb-2">
                                            <div className="flex items-center gap-2">
                                                {entry.status === 'accepted' ? (
                                                    <CheckCircle2 className="w-4 h-4 text-green-500" />
                                                ) : (
                                                    <XCircle className="w-4 h-4 text-red-400" />
                                                )}
                                                <span className="text-gray-500 text-[10px] uppercase font-bold tracking-tighter">
                                                    {new Date(entry.timestamp).toLocaleTimeString()}
                                                </span>
                                            </div>
                                            {onRevert && (
                                                <button
                                                    onClick={() => onRevert(entry)}
                                                    className="opacity-0 group-hover:opacity-100 p-1.5 rounded-lg bg-blue-50 text-blue-600 hover:bg-blue-100 transition-all shadow-sm"
                                                    title="Revert to this state"
                                                >
                                                    <RotateCcw className="w-3 h-3" />
                                                </button>
                                            )}
                                        </div>

                                        {entry.svg ? (
                                            <div
                                                className="bg-gray-50 rounded-lg p-2 mb-2 border border-gray-100 flex justify-center items-center h-32 overflow-hidden"
                                                dangerouslySetInnerHTML={{ __html: entry.svg }}
                                            />
                                        ) : (
                                            <div className="bg-gray-50 rounded-lg p-2 mb-2 font-mono text-[10px] text-gray-600 truncate select-all border border-gray-100">
                                                {entry.smiles}
                                            </div>
                                        )}

                                        <div className="text-[10px] text-gray-400 italic">
                                            {entry.status === 'accepted' ? 'Merged structure' : 'Rejected proposal'}
                                        </div>
                                    </div>
                                ))
                            )}
                        </div>
                        {entries.length > 0 && onClear && (
                            <div className="p-4 border-t border-gray-200 bg-white">
                                <button
                                    onClick={onClear}
                                    className="w-full flex items-center justify-center gap-2 py-2 rounded-lg text-gray-500 hover:text-red-500 hover:bg-red-50 transition-all text-sm font-medium border border-gray-100"
                                >
                                    <Trash2 className="w-4 h-4" />
                                    Clear History
                                </button>
                            </div>
                        )}
                    </div>
                )}
            </div>
        </div>
    );
};
