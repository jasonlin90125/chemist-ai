/** @type {import('tailwindcss').Config} */
export default {
    content: [
        "./index.html",
        "./src/**/*.{js,ts,jsx,tsx}",
    ],
    theme: {
        extend: {
            colors: {
                chemist: {
                    bg: '#0F1115',
                    panel: '#1A1D23',
                    accent: '#3B82F6',
                    text: '#E2E8F0',
                    muted: '#64748B',
                    success: '#22C55E', // Green for Added
                    danger: '#EF4444',  // Red for Removed
                }
            }
        },
    },
    plugins: [],
}
