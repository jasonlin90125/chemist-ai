import { Check, X } from 'lucide-react';
import { useEffect } from 'react';

interface DiffActionBarProps {
    onAccept: () => void;
    onReject: () => void;
}

export const DiffActionBar = ({ onAccept, onReject }: DiffActionBarProps) => {

    // Keyboard shortucts
    useEffect(() => {
        const handleKeyDown = (e: KeyboardEvent) => {
            if (e.key === 'Enter') onAccept();
            if (e.key === 'Escape') onReject();
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [onAccept, onReject]);

    return (
        <div className="absolute bottom-8 left-1/2 -translate-x-1/2 bg-chemist-panel border border-white/20 rounded-full px-6 py-2 flex items-center gap-4 shadow-xl z-50 animate-in fade-in slide-in-from-bottom-4">
            <span className="text-white/80 text-sm font-medium pr-2">Proposed Changes</span>

            <button
                onClick={onReject}
                className="flex items-center gap-2 px-3 py-1.5 rounded-full hover:bg-chemist-danger/20 text-chemist-danger transition-colors text-sm font-semibold"
            >
                <X className="w-4 h-4" />
                Reject
                <span className="text-[10px] opacity-60 ml-1 border border-current px-1 rounded">ESC</span>
            </button>

            <div className="w-px h-6 bg-white/20"></div>

            <button
                onClick={onAccept}
                className="flex items-center gap-2 px-3 py-1.5 rounded-full hover:bg-chemist-success/20 text-chemist-success transition-colors text-sm font-semibold"
            >
                <Check className="w-4 h-4" />
                Accept
                <span className="text-[10px] opacity-60 ml-1 border border-current px-1 rounded">ENT</span>
            </button>
        </div>
    );
};
