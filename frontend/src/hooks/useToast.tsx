import { useState } from 'react';
import { CheckCircle2, AlertCircle, X } from 'lucide-react';
import clsx from 'clsx';

export type ToastType = 'success' | 'info' | 'warning' | 'error';

interface Toast {
    id: string;
    message: string;
    type: ToastType;
}

export const useToast = () => {
    const [toasts, setToasts] = useState<Toast[]>([]);

    const showToast = (message: string, type: ToastType = 'success') => {
        const id = Math.random().toString(36).substring(2, 9);
        setToasts(prev => [...prev, { id, message, type }]);

        // Auto-remove after 3 seconds
        setTimeout(() => {
            setToasts(prev => prev.filter(t => t.id !== id));
        }, 3000);
    };

    return { toasts, showToast, removeToast: (id: string) => setToasts(prev => prev.filter(t => t.id !== id)) };
};

export const ToastContainer = ({ toasts, removeToast }: { toasts: Toast[], removeToast: (id: string) => void }) => {
    return (
        <div className="fixed top-6 left-1/2 -translate-x-1/2 z-[100] flex flex-col gap-3 pointer-events-none">
            {toasts.map(toast => (
                <div
                    key={toast.id}
                    className={clsx(
                        "pointer-events-auto flex items-center gap-3 px-4 py-3 rounded-xl shadow-2xl border min-w-[300px] animate-in slide-in-from-top-4 fade-in duration-300",
                        {
                            'bg-white border-gray-100 text-gray-800': toast.type === 'info',
                            'bg-green-50 border-green-100 text-green-800': toast.type === 'success',
                            'bg-amber-50 border-amber-100 text-amber-800': toast.type === 'warning',
                            'bg-red-50 border-red-100 text-red-800': toast.type === 'error'
                        }
                    )}
                >
                    {toast.type === 'success' ? (
                        <CheckCircle2 className="w-5 h-5 text-green-500" />
                    ) : (
                        <AlertCircle className={clsx("w-5 h-5", {
                            'text-blue-500': toast.type === 'info',
                            'text-amber-500': toast.type === 'warning',
                            'text-red-500': toast.type === 'error'
                        })} />
                    )}
                    <p className="text-sm font-semibold flex-1">{toast.message}</p>
                    <button
                        onClick={() => removeToast(toast.id)}
                        className="p-1 hover:bg-black/5 rounded-md transition-colors"
                    >
                        <X className="w-4 h-4 opacity-50" />
                    </button>
                </div>
            ))}
        </div>
    );
};
