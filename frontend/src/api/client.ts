import axios from 'axios';
import { VisualMolecule, EditRequest, SimpleEditRequest } from '../types/molecule';

const api = axios.create({
    baseURL: 'http://localhost:8000/api',
});

export const moleculeApi = {
    getInitial: async (): Promise<VisualMolecule> => {
        const res = await api.get('/molecule/ibrutinib');
        return res.data;
    },

    edit: async (req: EditRequest): Promise<VisualMolecule> => {
        const res = await api.post('/molecule/edit', req);
        return res.data;
    },

    simpleEdit: async (req: SimpleEditRequest): Promise<VisualMolecule> => {
        const res = await api.post('/molecule/simple-edit', req);
        return res.data;
    },

    multiEdit: async (req: SimpleEditRequest): Promise<VisualMolecule[]> => {
        const res = await api.post('/molecule/multi-edit', req);
        return res.data;
    },

    visualize: async (mol_block: string): Promise<VisualMolecule> => {
        const res = await api.post('/molecule/visualize', { mol_block });
        return res.data;
    },

    exportSDF: async (molecules: { mol_block: string, id: string }[]): Promise<Blob> => {
        const res = await api.post('/export/sdf', molecules, { responseType: 'blob' });
        return res.data;
    }
};
