export interface ClipboardEntry {
    id: string;
    smiles: string;
    mol_block: string;
    svg?: string;
    status: 'accepted' | 'rejected' | 'pending';
    timestamp: number;
    parent_id?: string;
}
