import { useState } from 'react';
import { Send, Sparkles, AlertCircle } from 'lucide-react';
import clsx from 'clsx';

interface ChatPanelProps {
    onSendPrompt: (prompt: string) => void;
    isLoading: boolean;
    error?: string | null;
}

export const ChatPanel = ({ onSendPrompt, isLoading, error }: ChatPanelProps) => {
    const [input, setInput] = useState("");

    const handleSubmit = (e: React.FormEvent) => {
        e.preventDefault();
        if (input.trim() && !isLoading) {
            onSendPrompt(input);
            setInput("");
        }
    };

    return (
        <div className="w-[400px] h-full flex flex-col border-r border-gray-200 bg-white z-20 shadow-xl">

            {/* Header */}
            <div className="p-4 border-b border-gray-100 flex items-center gap-2 bg-gray-50/50">
                <Sparkles className="text-blue-600 w-5 h-5" />
                <h2 className="text-gray-800 font-semibold tracking-tight">Chemist.ai</h2>
            </div>

            {/* History / Error / Info */}
            <div className="flex-1 p-4 overflow-y-auto space-y-4 bg-gray-50/30">
                {error && (
                    <div className="bg-red-50 border border-red-100 p-3 rounded-lg text-sm text-red-700 flex items-start gap-2 animate-in fade-in slide-in-from-top-2 duration-300">
                        <AlertCircle className="w-4 h-4 mt-0.5 shrink-0" />
                        <p>{error}</p>
                    </div>
                )}

                <div className="bg-blue-50 border border-blue-100 p-3 rounded-lg text-sm text-gray-700">
                    <p className="font-medium text-blue-900 mb-1">Welcome!</p>
                    <p>Select atoms on the canvas and describe your edit.</p>
                    <p className="mt-2 text-xs text-blue-600/80">Try: "Make this ring aromatic" or "Add a methyl group"</p>
                </div>
            </div>

            {/* Input */}
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
    );
};
