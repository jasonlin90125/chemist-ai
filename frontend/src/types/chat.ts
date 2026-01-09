export interface Message {
    id: string;
    role: 'user' | 'assistant' | 'system';
    content: string;
    timestamp: number;
    type?: 'text' | 'edit' | 'error';
    metadata?: {
        action?: string;
        proposalCount?: number;
    };
}
