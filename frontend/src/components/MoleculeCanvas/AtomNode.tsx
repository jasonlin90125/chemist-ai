import { Handle, Position } from '@xyflow/react';
import clsx from 'clsx';
import { Atom } from '../../types/molecule';

interface AtomNodeProps {
    data: Atom;
    selected?: boolean;
}

export const AtomNode = ({ data, selected }: AtomNodeProps) => {
    // Diff coloring
    let colorClass = 'text-white bg-chemist-panel';
    let borderClass = 'border-transparent';

    // Transparent style for Carbons (Skeletal look)
    const isCarbon = data.element === 'C';
    const isExisting = data.diff_state === 'EXISTING' || !data.diff_state;

    // Improved Logic:
    // If selected, show visual feedback regardless of element.
    // If not selected and Carbon, show skeletal (invisible bg) but maintain hover.

    if (data.diff_state === 'ADDED') {
        colorClass = 'text-chemist-success bg-green-900/20';
        borderClass = 'border-chemist-success';
    } else if (data.diff_state === 'REMOVED') {
        colorClass = 'text-chemist-danger bg-red-900/20';
        borderClass = 'border-chemist-danger';
    } else {
        // Default State
        if (isCarbon && isExisting) {
             colorClass = 'text-transparent bg-transparent hover:bg-white/10';
        } else {
             // Standard atom look
             colorClass = 'text-white bg-gray-800';
             // Add a subtle gradient or look
             // Actually, bg-gray-800 is fine, maybe a gradient?
             // let's keep it simple flat for now but cleaner
        }
    }

    if (selected) {
        borderClass = 'border-chemist-accent ring-2 ring-chemist-accent/50';
        if (isCarbon && isExisting) {
            // If selecting a carbon, make it visible
             colorClass = 'text-white bg-chemist-accent/20';
        }
    }

    // Element colors (cpk-ish)
    const elementColors: Record<string, string> = {
        O: 'text-red-400',
        N: 'text-blue-400',
        S: 'text-yellow-400',
        Cl: 'text-green-400',
        F: 'text-green-300',
        I: 'text-purple-400',
        Br: 'text-red-700',
        P: 'text-orange-400',
        C: 'text-gray-200'
    };

    const elementColor = elementColors[data.element] || 'text-gray-200';

    return (
        <div className={clsx(
            "w-10 h-10 rounded-full flex items-center justify-center border-2 transition-all shadow-md",
            colorClass,
            borderClass
        )}>
            <div className={clsx("font-bold text-lg select-none flex items-baseline", elementColor)}>
                {!isCarbon && (
                    <>
                        {data.element}
                        {data.implicit_h > 0 && (
                            <span className="text-xs ml-0.5 font-normal text-gray-400">
                                H{data.implicit_h > 1 ? <sub>{data.implicit_h}</sub> : null}
                            </span>
                        )}
                    </>
                )}
                {isCarbon && selected && (
                    <span className="text-xs text-white/50">C</span>
                )}
            </div>

            {/* Handles centered for bond radiation */}
            <Handle
                type="target"
                position={Position.Top}
                className="w-1 h-1 !bg-transparent !border-none"
                style={{ left: '50%', top: '50%', transform: 'translate(-50%, -50%)' }}
            />
            <Handle
                type="source"
                position={Position.Top}
                className="w-1 h-1 !bg-transparent !border-none"
                style={{ left: '50%', top: '50%', transform: 'translate(-50%, -50%)' }}
            />
        </div>
    );
};
