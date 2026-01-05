export interface Atom {
    id: number;
    element: string;
    x: number;
    y: number;
    charge: number;
    implicit_h: number;
    ui_state: "DEFAULT" | "SELECTED" | "DIMMED";
    diff_state: "EXISTING" | "ADDED" | "REMOVED";
}

export interface Bond {
    source: number;
    target: number;
    order: "SINGLE" | "DOUBLE" | "TRIPLE" | "AROMATIC";
    diff_state: "EXISTING" | "ADDED" | "REMOVED";
}

export interface VisualMolecule {
    molecule_id: string;
    atoms: Atom[];
    bonds: Bond[];
    mol_block?: string;
}

export interface EditRequest {
    current_molecule: VisualMolecule;
    user_prompt: string;
    selected_indices: number[];
}
