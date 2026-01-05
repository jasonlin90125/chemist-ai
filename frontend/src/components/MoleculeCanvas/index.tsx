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
    EdgeTypes
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
}

export const MoleculeCanvas = ({ molecule, onSelectionChange }: MoleculeCanvasProps) => {
    const [nodes, setNodes, onNodesChange] = useNodesState([]);
    const [edges, setEdges, onEdgesChange] = useEdgesState([]);

    // Convert VisualMolecule to React Flow Nodes/Edges
    useEffect(() => {
        if (!molecule) return;

        // Scale coordinates for better visibility (Angstroms -> Pixels)
        const SCALE = 60;

        // Center offset (approximate) - could be dynamic
        const CENTER_X = window.innerWidth / 2;
        const CENTER_Y = window.innerHeight / 2;

        const newNodes: Node[] = molecule.atoms.map((atom) => ({
            id: atom.id.toString(),
            type: 'atom',
            position: { x: atom.x * SCALE + CENTER_X, y: atom.y * SCALE + CENTER_Y }, // Simple projection
            data: { ...atom },
            draggable: false,
            connectable: false,
        }));

        const newEdges: Edge[] = molecule.bonds.map((bond, idx) => ({
            id: `e${bond.source}-${bond.target}-${idx}`,
            source: bond.source.toString(),
            target: bond.target.toString(),
            type: 'bond',
            data: { ...bond },
        }));

        setNodes(newNodes);
        setEdges(newEdges);
    }, [molecule, setNodes, setEdges]);

    const onNodeClick = useCallback((event: React.MouseEvent, node: Node) => {
        // Toggle selection logic could go here
        // For now, let's just propagate the click ID
        const atomId = parseInt(node.id);

        // We need to manage multi-selection state locally or upstream.
        // Let's assume onSelectionChange receives the NEW complete list.
        // For this MVP, let's just toggle the clicked one in the upstream state? 
        // Or we handle the "Blue Halo" here by updating local node data?

        // The requirement says: "Clicking an atom draws a blue halo (Active Selection)"
        // and "User selects atoms... + types prompt".

        // We'll trust the parent to handle the logic, we just notify click.
        // But to show visual feedback immediately, the parent needs to update the molecule prop
        // with `ui_state: SELECTED`.

        onSelectionChange([atomId]);
    }, [onSelectionChange]);

    return (
        <div className="w-full h-full bg-chemist-bg">
            <ReactFlow
                nodes={nodes}
                edges={edges}
                onNodesChange={onNodesChange}
                onEdgesChange={onEdgesChange}
                nodeTypes={nodeTypes}
                edgeTypes={edgeTypes}
                onNodeClick={onNodeClick}
                fitView
                minZoom={0.5}
                maxZoom={2}
                nodesDraggable={false}
                panOnDrag={true}
                zoomOnDoubleClick={false}
            >
                <Background gap={20} color="#333" />
                <Controls className="bg-chemist-panel border-chemist-muted text-white" />
            </ReactFlow>
        </div>
    );
};
