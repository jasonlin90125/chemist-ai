from rdkit import Chem
from rdkit.Chem import AllChem
from app.chemistry.molecule import VisualMoleculeBuilder, align_and_diff
from app.models import VisualMolecule
import copy

class ChemistryTools:
    @staticmethod
    def modify_atom(current_mol: Chem.Mol, atom_id: int, new_element: str) -> Chem.Mol:
        """
        Changes an atom to a different element.
        Note: In RDKit, we can edit the atom in place.
        """
        # Create editable mol
        rw_mol = Chem.RWMol(current_mol)
        
        # We need to find the atom with the matching ID in our tracking if possible.
        # But here we assume the RDKit indices match the input IDs for this session
        # (Since we rebuild the mol from JSON identically).
        if atom_id < rw_mol.GetNumAtoms():
            rw_mol.GetAtomWithIdx(atom_id).SetAtomicNum(
                Chem.GetPeriodicTable().GetAtomicNumber(new_element)
            )
            
            # Adjust hydrogens to satisfy valence
            rw_mol.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)
            rw_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rw_mol)
            
        return rw_mol.GetMol()

    @staticmethod
    def add_substructure(current_mol: Chem.Mol, anchor_atom_id: int, smiles_fragment: str) -> Chem.Mol:
        """
        Adds a fragment at a specific anchor point.
        """
        rw_mol = Chem.RWMol(current_mol)
        fragment = Chem.MolFromSmiles(smiles_fragment)
        
        if not fragment:
            return current_mol
            
        # Add fragment atoms
        frag_idx_map = {}
        for atom in fragment.GetAtoms():
            new_idx = rw_mol.AddAtom(atom)
            frag_idx_map[atom.GetIdx()] = new_idx
            
        # Add fragment bonds
        for bond in fragment.GetBonds():
            rw_mol.AddBond(
                frag_idx_map[bond.GetBeginAtomIdx()],
                frag_idx_map[bond.GetEndAtomIdx()],
                bond.GetBondType()
            )
            
        # Connect to anchor
        # We usually attach to the first atom of the fragment (simplistic)
        # or the one marked with an isotope or dummy.
        # For MVP, assume attaching to the first atom of fragment.
        if anchor_atom_id < rw_mol.GetNumAtoms():
            rw_mol.AddBond(
                anchor_atom_id,
                frag_idx_map[0],
                Chem.BondType.SINGLE
            )
            
        Chem.SanitizeMol(rw_mol)
        return rw_mol.GetMol()

    @staticmethod
    def apply_diff_metadata(original_mol: Chem.Mol, new_mol: Chem.Mol, visual_mol: VisualMolecule) -> VisualMolecule:
        """
        Compares the new visual molecule against the original to mark ADDED/EXISTING.
        """
        # A simple heuristic: 
        # Atoms that have the same coordinate (approx) and same element are EXISTING.
        # Others are ADDED.
        # This relies on the Alignment step having run before this.
        
        # Build spatial list of old atoms
        old_atoms = [] # list of (x, y, element)
        conf = original_mol.GetConformer()
        for i in range(original_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            old_atoms.append({
                "x": pos.x,
                "y": pos.y,
                "element": original_mol.GetAtomWithIdx(i).GetSymbol()
            })
            
        for atom in visual_mol.atoms:
            # unique matching to avoid double counting? For MVP, just find ANY match
            match = False
            for old in old_atoms:
                dx = atom.x - old["x"]
                dy = atom.y - old["y"]
                dist_sq = dx*dx + dy*dy
                if dist_sq < 0.25 and atom.element == old["element"]: # 0.5 Angstrom tolerance
                    match = True
                    break
            
            if match:
                atom.diff_state = "EXISTING"
            else:
                atom.diff_state = "ADDED"
                
        return visual_mol
