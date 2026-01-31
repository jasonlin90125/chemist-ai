from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from app.models import VisualMolecule, EditRequest, SimpleEditRequest
from app.chemistry.molecule import get_ibrutinib
from app.ai.agent import process_molecule_edit, process_simple_edit, process_multi_edit
import os
import io
from fastapi.responses import StreamingResponse
from rdkit import Chem
from dotenv import load_dotenv

load_dotenv()

app = FastAPI(title="Chemist.ai Backend")

@app.post("/api/molecule/multi-edit", response_model=list[VisualMolecule])
async def multi_edit_molecule(request: SimpleEditRequest):
    """
    Pre-calculate all unique symmetry variants for a simple edit.
    """
    return await process_multi_edit(request)

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

@app.post("/api/molecule/visualize", response_model=VisualMolecule)
async def visualize_molecule(request: dict):
    """
    Standardizes a mol_block and returns full metadata (SVG, SMILES, etc.)
    """
    from app.chemistry.molecule import VisualMoleculeBuilder
    from rdkit import Chem
    
    mol_block = request.get("mol_block")
    if not mol_block:
        raise HTTPException(status_code=400, detail="mol_block is required")
    
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
    if not mol:
        # Try SMILES if mol_block fails (unlikely from Ketcher but good fallback)
        mol = Chem.MolFromSmiles(mol_block)
        
    if not mol:
        raise HTTPException(status_code=400, detail="Invalid molecule data")
        
    return VisualMoleculeBuilder.mol_to_visual_json(mol)

@app.post("/api/molecule/edit", response_model=VisualMolecule)
async def edit_molecule(request: EditRequest):
    """
    Takes a molecule context + prompt, runs the AI agent, 
    and returns a diffed proposal molecule.
    """
    return await process_molecule_edit(request)

@app.post("/api/molecule/simple-edit", response_model=VisualMolecule)
async def simple_edit_molecule(request: SimpleEditRequest):
    """
    Direct tool execution for simple tasks (add, remove, search)
    without involving the LLM or images.
    """
    return await process_simple_edit(request)

@app.post("/api/export/sdf")
async def export_sdf(molecules: list[dict]):
    """
    Combines a list of molecules into a single SDF file.
    Expects [{"mol_block": "...", "id": "..."}]
    """
    try:
        output = io.StringIO()
        writer = Chem.SDWriter(output)
        
        for mol_data in molecules:
            mol_block = mol_data.get("mol_block")
            mol_id = mol_data.get("id")
            
            if mol_block:
                # Ketcher often provides V3000, ensure RDKit handles it
                mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
                if mol:
                    mol.SetProp("_Name", str(mol_id))
                    writer.write(mol)
        
        writer.close()
        
        # Convert string output to bytes for the response
        sdf_data = output.getvalue().encode('utf-8')
        return StreamingResponse(
            io.BytesIO(sdf_data),
            media_type="chemical/x-mdl-sdfile",
            headers={"Content-Disposition": "attachment; filename=apothecary.sdf"}
        )
    except Exception as e:
        print(f"SDF Export Error: {e}")
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run("app.main:app", host="0.0.0.0", port=8000, reload=True)
