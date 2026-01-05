from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from app.models import VisualMolecule, EditRequest
from app.chemistry.molecule import get_ibrutinib
from app.ai.agent import process_molecule_edit
import os
from dotenv import load_dotenv

load_dotenv()

app = FastAPI(title="Chemist.ai Backend")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

@app.get("/")
def health_check():
    return {"status": "ok", "version": "0.1.0"}

@app.get("/api/molecule/ibrutinib", response_model=VisualMolecule)
def get_initial_molecule():
    """Returns the Ibrutinib molecule to start the session."""
    return get_ibrutinib()

@app.post("/api/molecule/edit", response_model=VisualMolecule)
async def edit_molecule(request: EditRequest):
    """
    Takes a molecule context + prompt, runs the AI agent, 
    and returns a diffed proposal molecule.
    """
    return await process_molecule_edit(request)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.main:app", host="0.0.0.0", port=8000, reload=True)
