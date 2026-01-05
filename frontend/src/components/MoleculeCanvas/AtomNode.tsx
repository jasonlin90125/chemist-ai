import { Handle, Position } from '@xyflow/react';
import clsx from 'clsx';
import { Atom } from '../../types/molecule';

interface AtomNodeProps {
    data: Atom;
}

export const AtomNode = ({ data }: AtomNodeProps) => {
    const isSelected = data.ui_state === 'SELECTED';

    // Diff coloring
    let colorClass = 'text-white bg-chemist-panel';
    let borderClass = 'border-transparent';

    // Transparent style for Carbons (Skeletal look)
    if (data.element === 'C' && data.ui_state === 'DEFAULT' && data.diff_state === 'EXISTING') {
        colorClass = 'text-transparent bg-transparent hover:bg-white/10'; // Invisible but interactive
    }

    if (data.diff_state === 'ADDED') {
        colorClass = 'text-chemist-success bg-green-900/20';
        borderClass = 'border-chemist-success';
    } else if (data.diff_state === 'REMOVED') {
        colorClass = 'text-chemist-danger bg-red-900/20';
        borderClass = 'border-chemist-danger';
    } else if (isSelected) {
        borderClass = 'border-chemist-accent';
        colorClass = 'text-white bg-chemist-accent/20';
    }

    // Element colors (cpk-ish)
    const elementColors: Record<string, string> = {
        O: 'text-red-400',
        N: 'text-blue-400',
        S: 'text-yellow-400',
        Cl: 'text-green-400',
        F: 'text-green-300',
        C: 'text-gray-200'
    };

    const elementColor = elementColors[data.element] || 'text-gray-200';

    return (
        <div className={clsx(
            "w-10 h-10 rounded-full flex items-center justify-center border-2 transition-all shadow-lg",
            colorClass,
            borderClass
        )}>
            <div className={clsx("font-bold text-lg select-none flex items-baseline", elementColor)}>
                {data.element !== 'C' && (
                    <>
                        {data.element}
                        {data.implicit_h > 0 && (
                            <span className="text-xs ml-0.5">
                                H{data.implicit_h > 1 ? <sub>{data.implicit_h}</sub> : null}
                            </span>
                        )}
                    </>
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
