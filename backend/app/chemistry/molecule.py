from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from app.models import VisualMolecule, Atom, Bond
import uuid

class VisualMoleculeBuilder:
    @staticmethod
    def mol_to_visual_json(mol: Chem.Mol, mol_id: str = None) -> VisualMolecule:
        if not mol_id:
            mol_id = f"mol_{uuid.uuid4().hex[:8]}"
            
        atoms = []
        conf = mol.GetConformer()
        
        for atom in mol.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append(Atom(
                id=atom.GetIdx(),
                element=atom.GetSymbol(),
                x=pos.x,
                y=pos.y,
                charge=atom.GetFormalCharge(),
                implicit_h=atom.GetTotalNumHs(),
                ui_state="DEFAULT",
                diff_state="EXISTING"
            ))
            
        bonds = []
        for bond in mol.GetBonds():
            order_map = {
                Chem.BondType.SINGLE: "SINGLE",
                Chem.BondType.DOUBLE: "DOUBLE",
                Chem.BondType.TRIPLE: "TRIPLE",
                Chem.BondType.AROMATIC: "AROMATIC"
            }
            bonds.append(Bond(
                source=bond.GetBeginAtomIdx(),
                target=bond.GetEndAtomIdx(),
                order=order_map.get(bond.GetBondType(), "SINGLE"),
                diff_state="EXISTING"
            ))
            
        return VisualMolecule(
            molecule_id=mol_id,
            atoms=atoms,
            bonds=bonds
        )

    @staticmethod
    def visual_json_to_mol(visual_mol: VisualMolecule) -> Chem.Mol:
        # Create an editable molecule
        rw_mol = Chem.RWMol()
        
        # Add atoms
        atom_map = {} # visual_id -> rdkit_idx
        # We need to sort by ID to ensure consistent ordering if IDs are sequential 
        # but RDKit indices are 0-based.
        # Ideally, we trust the input order or re-index. 
        # For this simple version, we'll assume IDs map to indices if we insert in order,
        # but to be safe we track the mapping.
        
        sorted_atoms = sorted(visual_mol.atoms, key=lambda a: a.id)
        
        for atom_data in sorted_atoms:
            a = Chem.Atom(atom_data.element)
            a.SetFormalCharge(atom_data.charge)
            idx = rw_mol.AddAtom(a)
            atom_map[atom_data.id] = idx
            
        # Add bonds
        for bond in visual_mol.bonds:
            if bond.source in atom_map and bond.target in atom_map:
                order_map = {
                    "SINGLE": Chem.BondType.SINGLE,
                    "DOUBLE": Chem.BondType.DOUBLE,
                    "TRIPLE": Chem.BondType.TRIPLE,
                    "AROMATIC": Chem.BondType.AROMATIC
                }
                rw_mol.AddBond(
                    atom_map[bond.source],
                    atom_map[bond.target],
                    order_map.get(bond.order, Chem.BondType.SINGLE)
                )
                
        mol = rw_mol.GetMol()
        
        # Set coordinates
        conf = Chem.Conformer(mol.GetNumAtoms())
        for atom_data in sorted_atoms:
            idx = atom_map[atom_data.id]
            conf.SetAtomPosition(idx, (atom_data.x, atom_data.y, 0.0))
            
        mol.AddConformer(conf)
        return mol

def get_ibrutinib_smiles():
    return "C=CC(=O)N1CCC[C@H](C1)N2C3=NC=NC(=C3C(=N2)C4=CC=C(C=C4)OC5=CC=CC=C5)N"

def get_ibrutinib() -> VisualMolecule:
    smiles = get_ibrutinib_smiles()
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    Chem.Kekulize(mol)
    return VisualMoleculeBuilder.mol_to_visual_json(mol, "mol_ibrutinib")

def align_and_diff(original_mol: Chem.Mol, new_mol: Chem.Mol) -> VisualMolecule:
    # 1. Find MCS with strict comparison to prevent flipping/wrong alignment
    mcs = rdFMCS.FindMCS([original_mol, new_mol], timeout=5)
    
    if mcs.numAtoms > 0:
        try:
            # 2. Try partial structure alignment
            AllChem.GenerateDepictionMatching2DStructure(new_mol, original_mol)
        except Exception:
            # 3. Fallback: Index-based coordinate mapping
            coordMap = {}
            conf = original_mol.GetConformer()
            from rdkit.Geometry import Point2D
            for i in range(min(original_mol.GetNumAtoms(), new_mol.GetNumAtoms())):
                pos = conf.GetAtomPosition(i)
                coordMap[i] = Point2D(pos.x, pos.y)
            AllChem.Compute2DCoords(new_mol, coordMap=coordMap)
    else:
        # Fallback if no common substructure
        AllChem.Compute2DCoords(new_mol)

    # 3. Compute Diff States
    
    # Needs a robust graph isomorphism or substructure match to map atoms
    # For this MVP, we will try to assume RWMol operations might preserve some indices, 
    # but strictly speaking, RDKit re-indexes.
    # We will use the RDKit GetSubstructMatch to find the "Old" part in the "New" part.
    
    match = new_mol.GetSubstructMatch(original_mol)
    
    # This is a simplification. A real diff needs a complex graph edit distance or 
    # ID tracking if we maintained persistent IDs through RDKit edits.
    # Since we rebuild RDKit mols from JSON, IDs might be lost or shifted.
    # STRATEGY: 
    #   - Mark all atoms in new_mol as "ADDED" initially.
    #   - Atoms that map to original_mol via MCS/SubstructMatch are "EXISTING".
    #   - Atoms in original_mol that are NOT in the match are "REMOVED" (we need to return them too to visualize red).
    
    # To visualize "REMOVED" atoms, we actually need to return a "Super Graph" containing both.
    # OR, the client handles the "REMOVED" by rendering the old molecule's ghosts.
    # The requirement says: "Atom 8 (Carbon) turns Red... A new Nitrogen appears... in Green".
    # This implies the returned JSON might need to contain the union of atoms, OR we rely on IDs.
    
    # Let's try to map IDs effectively.
    
    try:
        Chem.Kekulize(new_mol)
    except Exception:
        pass # If kekulization fails (e.g. radicals), continue

    vis_mol = VisualMoleculeBuilder.mol_to_visual_json(new_mol)
    
    # In a real heavy implementation, we'd use the unique atom IDs to track life-cycle.
    # For now, let's just return the new molecule with everything 'EXISTING' to ensure basic aligned rendering first.
    # The 'Diff' logic with Green/Red is refined in tools.py where we explicitly know what we changed.
    
    return vis_mol
