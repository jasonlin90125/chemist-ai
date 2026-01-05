# Chemist.ai

The "Cursor" for Medicinal Chemistry. A natural language-driven molecular editor with "diff-based" visualizations.

## Project Status

- **Backend:** Fully operational (Python/FastAPI/RDKit).
- **Frontend:** Code complete (React/TS), but requires `node` and `npm` to run.
- **Verification:** Core logic verified via `verify_backend.py`.

## Quick Start (Backend)

The backend handles the core logic: loading molecules, AI tool execution, and diff generation.

1. **Install Dependencies:**
   ```bash
   cd backend
   python3 -m venv venv
   source venv/bin/activate
   pip install -r requirements.txt
   ```

2. **Configure Environment:**
   Add your OpenRouter API Key to `.env`:
   ```bash
   echo "OPENROUTER_API_KEY=sk-or-v1-..." > .env
   ```
   *(Without a key, the AI editor will likely return the original molecule or mock response).*

3. **Run the Server:**
   ```bash
   uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload
   ```
   Open `http://localhost:8000/docs` to see the API Swagger.

4. **Verify Logic:**
   Run the verification script to test "Ibrutinib" loading and "Phenyl -> Pyridine" transformation logic:
   ```bash
   python verify_backend.py
   ```
   *Expected Output:*
   ```text
   SUCCESS: Found 1 added atoms (Green).
   Sample Added Atom: N at ...
   ```

## Frontend Setup

**Prerequisite:** Ensure `node` (v18+) and `npm` are installed.

1. **Install & Run:**
   ```bash
   cd frontend
   npm install
   npm run dev
   ```

2. **Usage:**
   - Select atoms on the canvas.
   - Type "Convert this ring to pyridine".
   - See the Green (Added) / Red (Removed) diff.
   - Press Enter to accept.

## Architecture

- **Backend:** FastAPI + RDKit (Cheminformatics) + OpenAI/OpenRouter (AI Agent).
- **Frontend:** React + React Flow (Graph Visualization).
- **Key Feature:** `align_and_diff` ensures the molecule doesn't rotate disjointedly when modified, preserving user mental model.
