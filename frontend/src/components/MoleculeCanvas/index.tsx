import { useState, useCallback, useEffect } from 'react';
import {
    ReactFlow,
    Background,
    Controls,
    Node,
    Edge,
    useNodesState,
    useEdgesState,
    NodeTypes,
    EdgeTypes,
    OnSelectionChangeParams,
    SelectionMode
} from '@xyflow/react';
import '@xyflow/react/dist/style.css';

import { AtomNode } from './AtomNode';
import { BondEdge } from './BondEdge';
import { VisualMolecule } from '../../types/molecule';

const nodeTypes: NodeTypes = {
    atom: AtomNode,
};

const edgeTypes: EdgeTypes = {
    bond: BondEdge,
};

interface MoleculeCanvasProps {
    molecule: VisualMolecule | null;
    onSelectionChange: (selectedIds: number[]) => void;
    // We might need a callback to trigger edits directly from shortcuts
    onManualEdit?: (action: string, selectedIds: number[]) => void;
}

export const MoleculeCanvas = ({ molecule, onSelectionChange, onManualEdit }: MoleculeCanvasProps) => {
    // Initialize with explicit types to avoid 'never[]' error
    const [nodes, setNodes, onNodesChange] = useNodesState<Node>([]);
    const [edges, setEdges, onEdgesChange] = useEdgesState<Edge>([]);

    // Explicitly unused variables removed or used
    // const [selectionMode, setSelectionMode] = useState<SelectionMode>(SelectionMode.Partial);
    // We just use isSelectMode boolean to toggle props.

    const [isSelectMode, setIsSelectMode] = useState(false);
    const [selectedIds, setSelectedIds] = useState<number[]>([]);

    // Convert VisualMolecule to React Flow Nodes/Edges
    useEffect(() => {
        if (!molecule) return;

        const SCALE = 60;
        const CENTER_X = window.innerWidth / 2;
        const CENTER_Y = window.innerHeight / 2;

        const newNodes: Node[] = molecule.atoms.map((atom) => ({
            id: atom.id.toString(),
            type: 'atom',
            position: { x: atom.x * SCALE + CENTER_X, y: atom.y * SCALE + CENTER_Y },
            data: { ...atom },
            draggable: false,
            connectable: false,
            selected: atom.ui_state === 'SELECTED',
            selectable: true,
        }));

        const newEdges: Edge[] = molecule.bonds.map((bond, idx) => ({
            id: `e${bond.source}-${bond.target}-${idx}`,
            source: bond.source.toString(),
            target: bond.target.toString(),
            type: 'bond',
            data: { ...bond },
            selectable: false,
        }));

        setNodes(newNodes);
        setEdges(newEdges);
    }, [molecule, setNodes, setEdges]);

    // Handle Selection Change from React Flow
    const onSelectionChangeHandler = useCallback(({ nodes }: OnSelectionChangeParams) => {
        const ids = nodes.map(n => parseInt(n.id));
        setSelectedIds(ids);
        onSelectionChange(ids);
    }, [onSelectionChange]);

    // Toggle Mode
    const toggleMode = () => {
        setIsSelectMode(!isSelectMode);
    };

    // Keyboard Shortcuts
    useEffect(() => {
        const handleKeyDown = (e: KeyboardEvent) => {
            if (selectedIds.length === 0) return;

            // Avoid triggering when typing in an input
            if (document.activeElement?.tagName === 'INPUT' || document.activeElement?.tagName === 'TEXTAREA') return;

            if (e.key === 'Delete' || e.key === 'Backspace') {
                onManualEdit?.('delete', selectedIds);
            } else {
                // Check for single letter elements
                const key = e.key.toUpperCase();
                const validElements = ['C', 'N', 'O', 'S', 'F', 'P', 'I'];
                if (validElements.includes(key)) {
                    onManualEdit?.(`replace:${key}`, selectedIds);
                }
            }
        };

        window.addEventListener('keydown', handleKeyDown);
        return () => window.removeEventListener('keydown', handleKeyDown);
    }, [selectedIds, onManualEdit]);

    return (
        <div className="w-full h-full bg-chemist-bg relative">
            {/* Toolbar */}
            <div className="absolute top-4 right-16 z-10 flex gap-2">
                <button
                    onClick={toggleMode}
                    className={`px-3 py-1.5 rounded text-sm font-medium shadow-md transition-colors ${
                        isSelectMode
                            ? 'bg-chemist-accent text-white'
                            : 'bg-chemist-panel text-gray-300 hover:bg-white/10'
                    }`}
                >
                    {isSelectMode ? 'Lasso Mode' : 'Pan Mode'}
                </button>
                <div className="bg-chemist-panel/80 text-xs text-gray-400 px-3 py-1.5 rounded backdrop-blur-sm">
                    {selectedIds.length} selected
                </div>
            </div>

            <ReactFlow
                nodes={nodes}
                edges={edges}
                onNodesChange={onNodesChange}
                onEdgesChange={onEdgesChange}
                nodeTypes={nodeTypes}
                edgeTypes={edgeTypes}
                onSelectionChange={onSelectionChangeHandler}
                fitView
                minZoom={0.5}
                maxZoom={2}
                nodesDraggable={false}
                panOnDrag={!isSelectMode}
                selectionOnDrag={isSelectMode}
                selectionMode={SelectionMode.Partial}
                zoomOnDoubleClick={false}
                multiSelectionKeyCode={null}
            >
                <Background gap={20} color="#333" />
                <Controls className="bg-chemist-panel border-chemist-muted text-white" />
            </ReactFlow>
        </div>
    );
};
