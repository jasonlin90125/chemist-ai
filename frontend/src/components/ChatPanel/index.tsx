import { useState } from 'react';
import { Send, Sparkles } from 'lucide-react';
import clsx from 'clsx';

interface ChatPanelProps {
    onSendPrompt: (prompt: string) => void;
    isLoading: boolean;
}

export const ChatPanel = ({ onSendPrompt, isLoading }: ChatPanelProps) => {
    const [input, setInput] = useState("");

    const handleSubmit = (e: React.FormEvent) => {
        e.preventDefault();
        if (input.trim() && !isLoading) {
            onSendPrompt(input);
            setInput("");
        }
    };

    return (
        <div className="absolute left-4 top-4 bottom-4 w-80 bg-chemist-panel/90 backdrop-blur-md rounded-xl border border-white/10 flex flex-col shadow-2xl z-10">

            {/* Header */}
            <div className="p-4 border-b border-white/10 flex items-center gap-2">
                <Sparkles className="text-chemist-accent w-5 h-5" />
                <h2 className="text-white font-semibold">Chemist.ai</h2>
            </div>

            {/* History (Placeholder for now) */}
            <div className="flex-1 p-4 overflow-y-auto space-y-4">
                <div className="bg-white/5 p-3 rounded-lg text-sm text-chemist-muted">
                    <p>Welcome! Select atoms on the canvas and describe your edit.</p>
                    <p className="mt-2 text-xs opacity-50">Try: "Make this ring aromatic" or "Add a methyl group"</p>
                </div>
            </div>

            {/* Input */}
            <div className="p-4 border-t border-white/10">
                <form onSubmit={handleSubmit} className="relative">
                    <input
                        type="text"
                        className="w-full bg-black/40 text-white rounded-lg pl-4 pr-10 py-3 focus:outline-none focus:ring-2 focus:ring-chemist-accent placeholder-chemist-muted border border-white/5 transition-all"
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
                            isLoading ? "opacity-50 cursor-wait" : "hover:bg-white/10 text-chemist-accent"
                        )}
                    >
                        <Send className="w-4 h-4" />
                    </button>
                </form>
            </div>
        </div>
    );
};
