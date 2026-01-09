from typing import List, Optional, Literal
from pydantic import BaseModel

class Atom(BaseModel):
    id: int
    element: str
    x: float
    y: float
    charge: int = 0
    implicit_h: int = 0
    atom_map: Optional[int] = None
    ui_state: Literal["DEFAULT", "SELECTED", "DIMMED"] = "DEFAULT"
    diff_state: Literal["EXISTING", "ADDED", "REMOVED"] = "EXISTING"

class Bond(BaseModel):
    source: int
    target: int
    order: Literal["SINGLE", "DOUBLE", "TRIPLE", "AROMATIC"]
    diff_state: Literal["EXISTING", "ADDED", "REMOVED"] = "EXISTING"

class VisualMolecule(BaseModel):
    molecule_id: str
    atoms: List[Atom]
    bonds: List[Bond]
    mol_block: Optional[str] = None
    smiles: Optional[str] = None
    svg: Optional[str] = None

class EditRequest(BaseModel):
    current_molecule: VisualMolecule
    user_prompt: str
    selected_indices: List[int]
    selected_maps: List[int] = []

class SimpleEditRequest(BaseModel):
    action: str  # e.g., "add_substructure", "remove_atoms"
    current_molecule: VisualMolecule
    selected_indices: List[int]
    selected_maps: List[int] = []
    parameters: dict
