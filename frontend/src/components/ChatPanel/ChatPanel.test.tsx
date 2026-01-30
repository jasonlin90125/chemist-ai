import { describe, it, expect, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import { ChatPanel } from './index';
import { Message } from '../../types/chat';

// Mock ResizeObserver properly
global.ResizeObserver = class ResizeObserver {
  observe() {}
  unobserve() {}
  disconnect() {}
};

// Mock HTMLElement.offsetParent
Object.defineProperty(HTMLElement.prototype, 'offsetParent', {
  get() {
    return this.parentNode;
  },
});

describe('ChatPanel Performance', () => {
  const generateMessages = (count: number): Message[] => {
    return Array.from({ length: count }, (_, i) => ({
      id: `msg-${i}`,
      role: i % 2 === 0 ? 'user' : 'assistant',
      content: `Message content ${i}`,
      timestamp: Date.now(),
      type: 'text',
    }));
  };

  const defaultProps = {
    onSendPrompt: vi.fn(),
    isLoading: false,
    activeTab: 'chat' as const,
    onTabChange: vi.fn(),
    messages: [],
    entries: [],
  };

  it('renders only visible messages in the DOM (Virtualization)', () => {
    const messageCount = 1000;
    const messages = generateMessages(messageCount);

    const { container } = render(
      <ChatPanel
        {...defaultProps}
        messages={messages}
      />
    );

    // Look for the message containers
    const messageElements = container.querySelectorAll('.animate-in.fade-in.slide-in-from-bottom-2');

    // With virtualization, we should NOT see all 1000 messages.
    // In JSDOM without layout, Virtuoso might render 0 or a few items.
    console.log(`Rendered ${messageElements.length} messages out of ${messageCount}`);

    expect(messageElements.length).toBeLessThan(messageCount);
    expect(messageElements.length).toBeLessThan(100);
  });

  it('renders welcome message when chat is empty', () => {
     render(
      <ChatPanel
        {...defaultProps}
        messages={[]}
      />
    );

    // Check for "Ready for Design" text
    expect(screen.getByText('Ready for Design')).toBeInTheDocument();
  });
});
