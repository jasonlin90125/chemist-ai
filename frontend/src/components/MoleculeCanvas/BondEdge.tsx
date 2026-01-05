import { BaseEdge, EdgeLabelRenderer, EdgeProps, getStraightPath } from '@xyflow/react';
import clsx from 'clsx';
import { Bond } from '../../types/molecule';

export const BondEdge = ({
    id,
    sourceX,
    sourceY,
    targetX,
    targetY,
    data,
    markerEnd,
}: EdgeProps<Bond>) => {
    const [edgePath, labelX, labelY] = getStraightPath({
        sourceX,
        sourceY,
        targetX,
        targetY,
    });

    const diffState = data?.diff_state || 'EXISTING';
    let strokeColor = '#64748B'; // Muted slate
    let strokeWidth = 2;

    if (diffState === 'ADDED') {
        strokeColor = '#22C55E';
        strokeWidth = 3;
    } else if (diffState === 'REMOVED') {
        strokeColor = '#EF4444';
        strokeWidth = 3;
        // Dashed for removed?
    }

    // Bond Visualization Logic
    // Bond Visualization Logic
    const isDouble = data?.order === 'DOUBLE';
    const isTriple = data?.order === 'TRIPLE';

    // Calculate normal vector for offsets
    const dx = targetX - sourceX;
    const dy = targetY - sourceY;
    const len = Math.sqrt(dx * dx + dy * dy);
    // Unit normal vector
    const nx = len ? -dy / len : 0;
    const ny = len ? dx / len : 0;

    const OFFSET = 4; // Increased from 2.5 to 4 for better visibility

    // Helper to generate offset path
    const getOffsetPath = (offset: number) => {
        return `M ${sourceX + nx * offset} ${sourceY + ny * offset} L ${targetX + nx * offset} ${targetY + ny * offset}`;
    };

    const baseStyle = {
        stroke: strokeColor,
        strokeWidth: 2,
        opacity: diffState === 'REMOVED' ? 0.3 : 1
    };

    if (isDouble) {
        return (
            <g>
                <path d={getOffsetPath(OFFSET)} style={baseStyle} fill="none" />
                <path d={getOffsetPath(-OFFSET)} style={baseStyle} fill="none" />
            </g>
        );
    }

    if (isTriple) {
        return (
            <g>
                <path d={edgePath} style={baseStyle} fill="none" />
                <path d={getOffsetPath(OFFSET * 1.5)} style={baseStyle} fill="none" />
                <path d={getOffsetPath(-OFFSET * 1.5)} style={baseStyle} fill="none" />
            </g>
        );
    }

    // Default Single
    return (
        <BaseEdge
            path={edgePath}
            markerEnd={markerEnd}
            style={{ ...baseStyle, strokeWidth: diffState === 'ADDED' ? 3 : 2 }}
        />
    );
};
