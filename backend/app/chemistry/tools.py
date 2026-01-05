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
        rw_mol = Chem.RWMol(current_mol)
        
        if atom_id < rw_mol.GetNumAtoms():
            rw_mol.GetAtomWithIdx(atom_id).SetAtomicNum(
                Chem.GetPeriodicTable().GetAtomicNumber(new_element)
            )
            rw_mol.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)
            rw_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rw_mol)
            
        return rw_mol.GetMol()

    @staticmethod
    def replace_atoms(current_mol: Chem.Mol, atom_ids: list[int], new_element: str) -> Chem.Mol:
        """
        Replaces multiple atoms with a new element.
        """
        rw_mol = Chem.RWMol(current_mol)

        for atom_id in atom_ids:
            if atom_id < rw_mol.GetNumAtoms():
                rw_mol.GetAtomWithIdx(atom_id).SetAtomicNum(
                    Chem.GetPeriodicTable().GetAtomicNumber(new_element)
                )
                rw_mol.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)

        rw_mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(rw_mol)
        except Exception:
            pass

        return rw_mol.GetMol()

    @staticmethod
    def remove_atoms(current_mol: Chem.Mol, atom_ids: list[int]) -> Chem.Mol:
        """
        Removes atoms from the molecule.
        Must remove in descending order of indices to avoid shifting problems.
        """
        rw_mol = Chem.RWMol(current_mol)

        # Sort descending
        sorted_ids = sorted(atom_ids, reverse=True)

        for atom_id in sorted_ids:
            if atom_id < rw_mol.GetNumAtoms():
                rw_mol.RemoveAtom(atom_id)

        try:
            Chem.SanitizeMol(rw_mol)
        except Exception:
            pass # Sanitize might fail on fragments, but we return what we have

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
            
        frag_idx_map = {}
        for atom in fragment.GetAtoms():
            new_idx = rw_mol.AddAtom(atom)
            frag_idx_map[atom.GetIdx()] = new_idx
            
        for bond in fragment.GetBonds():
            rw_mol.AddBond(
                frag_idx_map[bond.GetBeginAtomIdx()],
                frag_idx_map[bond.GetEndAtomIdx()],
                bond.GetBondType()
            )
            
        if anchor_atom_id < rw_mol.GetNumAtoms():
            # Attach to first atom of fragment
            rw_mol.AddBond(
                anchor_atom_id,
                frag_idx_map[0],
                Chem.BondType.SINGLE
            )
            
        try:
            Chem.SanitizeMol(rw_mol)
        except Exception:
            pass

        return rw_mol.GetMol()

    @staticmethod
    def replace_substructure(current_mol: Chem.Mol, atom_ids: list[int], smiles_fragment: str) -> Chem.Mol:
        """
        Removes the specified atoms and attaches the new fragment to the neighbor of the first removed atom (simplistic).
        A better implementation would identify attachment points.
        """
        # 1. Identify neighbors of the substructure that are NOT in the substructure.
        # These will be the attachment points.
        neighbors = set()
        atoms_to_remove = set(atom_ids)

        # We need a map of which removed atom connects to which outside atom
        # But for MVP, let's just find the first neighbor and attach the fragment there.
        # Then remove the atoms.

        # This is complex because removing atoms destroys bonds.
        # Strategy:
        # 1. Find attachment point(s).
        # 2. Remove atoms.
        # 3. Add fragment.
        # 4. Bond fragment to attachment point.

        rw_mol = Chem.RWMol(current_mol)

        attachment_point = -1

        # Find an attachment point (neighbor of the substructure)
        for atom_id in atom_ids:
            if atom_id >= rw_mol.GetNumAtoms(): continue
            atom = rw_mol.GetAtomWithIdx(atom_id)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in atoms_to_remove:
                    attachment_point = n_idx
                    break
            if attachment_point != -1:
                break

        # Remove atoms
        sorted_ids = sorted(atom_ids, reverse=True)
        for atom_id in sorted_ids:
            if atom_id < rw_mol.GetNumAtoms():
                rw_mol.RemoveAtom(atom_id)

        # Indices shifted. We need to track where attachment_point went?
        # RemoveAtom shifts indices > atom_id down by 1.
        # So we need to adjust attachment_point.
        # Re-calculating attachment point is hard after removal if we don't track it.

        # Actually, simpler way: use ReplaceSubstructs if we can pattern match.
        # But we have IDs.

        # Let's adjust attachment_point index.
        # If we remove an atom with index < attachment_point, attachment_point decreases.
        current_attachment = attachment_point
        if current_attachment != -1:
            for atom_id in sorted_ids:
                if atom_id < current_attachment:
                    current_attachment -= 1
                elif atom_id == current_attachment:
                    current_attachment = -1 # Should not happen based on logic above
                    break

        # Add fragment
        if current_attachment != -1:
            fragment = Chem.MolFromSmiles(smiles_fragment)
            if fragment:
                frag_idx_map = {}
                for atom in fragment.GetAtoms():
                    new_idx = rw_mol.AddAtom(atom)
                    frag_idx_map[atom.GetIdx()] = new_idx

                for bond in fragment.GetBonds():
                    rw_mol.AddBond(
                        frag_idx_map[bond.GetBeginAtomIdx()],
                        frag_idx_map[bond.GetEndAtomIdx()],
                        bond.GetBondType()
                    )

                # Attach
                rw_mol.AddBond(
                    current_attachment,
                    frag_idx_map[0],
                    Chem.BondType.SINGLE
                )

        try:
            Chem.SanitizeMol(rw_mol)
        except Exception:
            pass

        return rw_mol.GetMol()

    @staticmethod
    def get_mapped_smiles(mol: Chem.Mol) -> str:
        """
        Returns SMILES with atom indices mapped.
        """
        # Create a copy to not modify original properties
        m = Chem.Mol(mol)
        for atom in m.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx())
        return Chem.MolToSmiles(m)

    @staticmethod
    def apply_diff_metadata(original_mol: Chem.Mol, new_mol: Chem.Mol, visual_mol: VisualMolecule) -> VisualMolecule:
        """
        Compares the new visual molecule against the original to mark ADDED/EXISTING.
        """
        # A simple heuristic: 
        # Atoms that have the same coordinate (approx) and same element are EXISTING.
        # Others are ADDED.
        
        old_atoms = []
        conf = original_mol.GetConformer()
        for i in range(original_mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            old_atoms.append({
                "x": pos.x,
                "y": pos.y,
                "element": original_mol.GetAtomWithIdx(i).GetSymbol()
            })
            
        for atom in visual_mol.atoms:
            match = False
            for old in old_atoms:
                dx = atom.x - old["x"]
                dy = atom.y - old["y"]
                dist_sq = dx*dx + dy*dy
                if dist_sq < 0.25 and atom.element == old["element"]:
                    match = True
                    break
            
            if match:
                atom.diff_state = "EXISTING"
            else:
                atom.diff_state = "ADDED"
                
        return visual_mol
