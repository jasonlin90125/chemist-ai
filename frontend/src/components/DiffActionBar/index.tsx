import { Check, X, ChevronLeft, ChevronRight, Hash } from 'lucide-react';
import { useEffect } from 'react';

interface DiffActionBarProps {
    onAccept: () => void;
    onReject: () => void;
    onCycle: (direction: number) => void;
}

export const DiffActionBar = ({ onAccept, onReject, onCycle }: DiffActionBarProps) => {

    // Keyboard shortucts
    useEffect(() => {
        const handleKeyDown = (e: KeyboardEvent) => {
            if (e.key === 'Enter') onAccept();
            if (e.key === 'Escape') onReject();
            if (e.key === 'ArrowUp') onCycle(1);
            if (e.key === 'ArrowDown') onCycle(-1);
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [onAccept, onReject, onCycle]);

    return (
        <div className="absolute bottom-8 left-1/2 -translate-x-1/2 bg-chemist-panel border border-white/20 rounded-full px-6 py-2 flex items-center gap-4 shadow-xl z-50 animate-in fade-in slide-in-from-bottom-4">
            <span className="text-white/80 text-sm font-medium pr-2">Proposed Changes</span>

            <button
                onClick={onReject}
                className="flex items-center gap-2 px-3 py-1.5 rounded-full hover:bg-chemist-danger/20 text-chemist-danger transition-colors text-sm font-semibold"
            >
                <X className="w-4 h-4" />
                Reject
            </button>

            <div className="w-px h-6 bg-white/20"></div>

            <div className="flex items-center gap-1">
                <button
                    onClick={() => onCycle(-1)}
                    className="p-1.5 rounded-full hover:bg-white/10 text-white/80 transition-colors"
                    title="Previous Variant (Up Arrow)"
                >
                    <ChevronLeft className="w-4 h-4" />
                </button>
                <span className="text-white/40"><Hash className="w-3 h-3" /></span>
                <button
                    onClick={() => onCycle(1)}
                    className="p-1.5 rounded-full hover:bg-white/10 text-white/80 transition-colors"
                    title="Next Variant (Down Arrow)"
                >
                    <ChevronRight className="w-4 h-4" />
                </button>
            </div>

            <div className="w-px h-6 bg-white/20"></div>

            <button
                onClick={onAccept}
                className="flex items-center gap-2 px-3 py-1.5 rounded-full hover:bg-chemist-success/20 text-chemist-success transition-colors text-sm font-semibold"
            >
                <Check className="w-4 h-4" />
                Accept
            </button>
        </div>
    );
};
