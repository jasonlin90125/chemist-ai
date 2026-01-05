import { BaseEdge, EdgeProps, getBezierPath } from '@xyflow/react';
import { Bond } from '../../types/molecule';

export const BondEdge = ({
    sourceX,
    sourceY,
    targetX,
    targetY,
    data,
    markerEnd,
    style
}: EdgeProps) => {
    // We removed 'id', 'labelX', 'labelY', 'EdgeLabelRenderer', 'clsx' to fix unused var errors
    // We added 'style' to satisfy EdgeProps requirement (though we override it)

    const [edgePath] = getBezierPath({
        sourceX,
        sourceY,
        targetX,
        targetY,
    });

    // Handle data type safely
    const bondData = data as unknown as Bond | undefined;
    const diff_state = bondData?.diff_state || 'EXISTING';
    const order = bondData?.order || 'SINGLE';

    let strokeColor = '#4B5563'; // gray-600
    let strokeWidth = 2;
    let strokeDasharray = undefined;

    if (diff_state === 'ADDED') {
        strokeColor = '#10B981'; // green-500
        strokeWidth = 3;
    } else if (diff_state === 'REMOVED') {
        strokeColor = '#EF4444'; // red-500
        strokeDasharray = '5,5';
    }

    // Double bond styling simulation
    if (order === 'DOUBLE') {
        strokeWidth = 4;
    } else if (order === 'TRIPLE') {
        strokeWidth = 6;
    }

    return (
        <BaseEdge
            path={edgePath}
            markerEnd={markerEnd}
            style={{
                ...style, // preserve incoming style (selection etc) but override colors
                stroke: strokeColor,
                strokeWidth,
                strokeDasharray,
                opacity: diff_state === 'REMOVED' ? 0.5 : 1
            }}
        />
    );
};
