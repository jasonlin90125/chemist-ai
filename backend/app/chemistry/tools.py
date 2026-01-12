from rdkit import Chem
from rdkit.Chem import AllChem
from app.chemistry.molecule import VisualMoleculeBuilder, align_and_diff
from app.chemistry.constants import FUNCTIONAL_GROUPS
from app.chemistry.substructures import SUBSTRUCTURE_PATTERNS
from app.models import VisualMolecule
import copy
from typing import Optional, Union, List, Tuple

class ChemistryTools:
    @staticmethod
    def _resolve_atom_idx(mol: Chem.Mol, provided_id: int, look_for_map: bool = True) -> Optional[int]:
        """
        Resolves an atom ID to its current RDKit index.
        If look_for_map is True, treats provided_id as a Persistent Map Number.
        Otherwise, treats it as a direct index and validates range.
        """
        if look_for_map:
            for atom in mol.GetAtoms():
                if atom.GetAtomMapNum() == provided_id:
                    return atom.GetIdx()
            return None
        else:
            if provided_id is not None and provided_id >= 0 and provided_id < mol.GetNumAtoms():
                return provided_id
            return None

    @staticmethod
    def _get_fragment_and_attachment(smiles_or_name: str, variant_idx: int = 0) -> tuple[Chem.Mol, int]:
        """
        Helper to get RDKit mol and the index of the atom to attach.
        Identifies all symmetry-unique attachment possibilities if variant_idx > 0.
        Returns (mol_without_dummy, attachment_idx).
        """
        # 1. Lookup in dictionaries
        query = smiles_or_name.lower()
        smiles = FUNCTIONAL_GROUPS.get(query)
        if not smiles:
            smiles = SUBSTRUCTURE_PATTERNS.get(query, smiles_or_name)
        
        # 2. Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            mol = Chem.MolFromSmiles(smiles.capitalize())
            if not mol:
                return None, -1
            
        # 3. Handle explicit dummy atom '*' (the specific intended attachment)
        dummy_idx = -1
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummy_idx = atom.GetIdx()
                break
        
        # If we have a dummy AND variant_idx is 0, we use the intended attachment
        if dummy_idx != -1 and variant_idx == 0:
            neighbors = mol.GetAtomWithIdx(dummy_idx).GetNeighbors()
            attachment_idx = neighbors[0].GetIdx() if neighbors else 0
            
            nm = Chem.RWMol(mol)
            nm.RemoveAtom(dummy_idx)
            if dummy_idx < attachment_idx:
                attachment_idx -= 1
            return nm.GetMol(), attachment_idx

        # 4. Otherwise (no dummy or cycling requested), find all unique atoms
        # Remove dummy first if it exists to get clean symmetry ranks
        if dummy_idx != -1:
            nm = Chem.RWMol(mol)
            nm.RemoveAtom(dummy_idx)
            mol = nm.GetMol()

        # Calculate symmetry ranks to find unique attachment points
        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        unique_ranks = sorted(list(set(ranks)))
        
        # Map each unique rank to the first atom index that has it
        unique_atom_indices = []
        for r in unique_ranks:
            for i, rank in enumerate(ranks):
                if rank == r:
                    unique_atom_indices.append(i)
                    break
        
        # Select the attachment index based on variant_idx
        attachment_idx = unique_atom_indices[variant_idx % len(unique_atom_indices)]
        
        return mol, attachment_idx

    @staticmethod
    def get_variants_count(smiles_or_name: str) -> int:
        """
        Returns the number of symmetry-unique attachment points for a fragment.
        """
        query = smiles_or_name.lower()
        smiles = FUNCTIONAL_GROUPS.get(query)
        if not smiles:
            smiles = SUBSTRUCTURE_PATTERNS.get(query, smiles_or_name)
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            mol = Chem.MolFromSmiles(smiles.capitalize())
            if not mol: return 1
            
        dummy_idx = -1
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummy_idx = atom.GetIdx()
                break
        
        if dummy_idx != -1:
            nm = Chem.RWMol(mol)
            nm.RemoveAtom(dummy_idx)
            mol = nm.GetMol()

        ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
        return len(set(ranks))

    @staticmethod
    def modify_atom(current_mol: Chem.Mol, atom_id: int, new_element: str) -> Chem.Mol:
        """
        Changes an atom to a different element.
        """
        rw_mol = Chem.RWMol(current_mol)
        
        if atom_id >= rw_mol.GetNumAtoms():
            raise ValueError(f"Atom ID {atom_id} is out of range. Molecule has {rw_mol.GetNumAtoms()} atoms.")

        rw_mol.GetAtomWithIdx(atom_id).SetAtomicNum(
            Chem.GetPeriodicTable().GetAtomicNumber(new_element)
        )
        rw_mol.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)
        rw_mol.UpdatePropertyCache(strict=False)
        
        try:
            Chem.SanitizeMol(rw_mol)
        except Exception as e:
            raise ValueError(f"Modification failed during sanitization: {str(e)}")
            
        return rw_mol.GetMol()

    @staticmethod
    def replace_atoms(current_mol: Chem.Mol, atom_ids: list[int], new_element: str) -> Chem.Mol:
        """
        Replaces multiple atoms with a new element.
        """
        rw_mol = Chem.RWMol(current_mol)

        for atom_id in atom_ids:
            if atom_id >= rw_mol.GetNumAtoms():
                raise ValueError(f"Atom ID {atom_id} is out of range.")
            rw_mol.GetAtomWithIdx(atom_id).SetAtomicNum(
                Chem.GetPeriodicTable().GetAtomicNumber(new_element)
            )
            rw_mol.GetAtomWithIdx(atom_id).SetNumExplicitHs(0)

        rw_mol.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(rw_mol)
        except Exception as e:
            raise ValueError(f"Batch replacement failed during sanitization: {str(e)}")

        return rw_mol.GetMol()

    @staticmethod
    def remove_atoms(current_mol: Chem.Mol, atom_ids: list[int]) -> Chem.Mol:
        """
        Removes atoms from the molecule.
        """
        rw_mol = Chem.RWMol(current_mol)
        sorted_ids = sorted(atom_ids, reverse=True)

        for atom_id in sorted_ids:
            if atom_id >= rw_mol.GetNumAtoms():
                continue # or raise error? consistent with others let's raise
            rw_mol.RemoveAtom(atom_id)

        try:
            Chem.SanitizeMol(rw_mol)
        except Exception as e:
            # For removal, we might occasionally have fragments. 
            # If sanitization fails, we still try to return the result but with a warning?
            # Actually, let's just return if it contains valid atoms.
            pass

        return rw_mol.GetMol()

    @staticmethod
    def add_substructure(current_mol: Chem.Mol, anchor_atom_id: int, smiles_fragment: str, variant_idx: int = 0) -> Chem.Mol:
        """
        Adds a fragment at a specific anchor point.
        """
        # Explicate Hs to ensure we can substitute them if needed
        mol_with_hs = Chem.AddHs(current_mol)
        rw_mol = Chem.RWMol(mol_with_hs)
        
        if anchor_atom_id >= rw_mol.GetNumAtoms():
             raise ValueError(f"Anchor atom ID {anchor_atom_id} is out of range.")

        fragment, frag_attach_idx = ChemistryTools._get_fragment_and_attachment(smiles_fragment, variant_idx)
        
        if not fragment:
            raise ValueError(f"Failed to parse or find fragment: {smiles_fragment}")
        
        # Check if anchor atom has implicit hydrogens we can replace
        anchor_atom = rw_mol.GetAtomWithIdx(anchor_atom_id)
        
        # Update property cache to get accurate H count
        # Update property cache to get accurate H count
        rw_mol.UpdatePropertyCache(strict=False)
        
        # Since we added Hs, we only check explicit neighbors
        explicit_h_neighbors = [n for n in anchor_atom.GetNeighbors() if n.GetAtomicNum() == 1]
        
        # If we need to make room (valence check would happen here in reality, 
        # but for now we opportunistically remove H if one exists and we are attaching something)
        if explicit_h_neighbors:
            # Remove one explicit H to make room
            h_to_remove = explicit_h_neighbors[0].GetIdx()
            rw_mol.RemoveAtom(h_to_remove)
            
            # Since we removed an atom, indices might shift if the removed atom had a lower index
            # But AddHs usually adds Hs at the end, so anchor_atom_id (original atom) should be safe.
            # However, if we remove an atom < anchor_atom_id, we must adjust.
            if h_to_remove < anchor_atom_id:
                anchor_atom_id -= 1
            
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
            
        # Connect fragment to anchor
        if frag_attach_idx != -1:
            rw_mol.AddBond(
                anchor_atom_id,
                frag_idx_map[frag_attach_idx],
                Chem.BondType.SINGLE
            )
            
        try:
            rw_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rw_mol)
            # Remove explicit Hs (clean up)
            final_mol = Chem.RemoveHs(rw_mol.GetMol())
            AllChem.Compute2DCoords(final_mol)
        except Exception as e:
            raise ValueError(f"Adding substructure failed: {str(e)}")

        return final_mol

    @staticmethod
    def replace_substructure(current_mol: Chem.Mol, atom_ids: list[int], smiles_fragment: str, variant_idx: int = 0) -> Chem.Mol:
        """
        Removes the specified atoms and attaches the new fragment.
        """
        rw_mol = Chem.RWMol(current_mol)
        atoms_to_remove = set(atom_ids)
        
        # 1. Identify all boundary bonds (one atom in remove_set, one not)
        # Store as (neighbor_index, bond_type)
        boundary_connections = []
        for atom_id in atom_ids:
            if atom_id >= rw_mol.GetNumAtoms(): continue
            atom = rw_mol.GetAtomWithIdx(atom_id)
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if n_idx not in atoms_to_remove:
                    bond = rw_mol.GetBondBetweenAtoms(atom_id, n_idx)
                    boundary_connections.append((n_idx, bond.GetBondType()))

        # 2. Remove atoms
        sorted_ids = sorted(atom_ids, reverse=True)
        for atom_id in sorted_ids:
            if atom_id >= rw_mol.GetNumAtoms():
                 raise ValueError(f"Atom ID {atom_id} is out of range.")
            rw_mol.RemoveAtom(atom_id)

        # 3. Adjust boundary connection indices (since atoms were removed)
        adjusted_connections = []
        for n_idx, b_type in boundary_connections:
            new_n_idx = n_idx
            for removed_id in sorted_ids:
                if removed_id < n_idx:
                    new_n_idx -= 1
                elif removed_id == n_idx:
                    # Neighbor was also removed? Should be impossible due to check above
                    new_n_idx = -1
                    break
            if new_n_idx != -1:
                adjusted_connections.append((new_n_idx, b_type))

        # 4. Add new fragment
        fragment, frag_attach_idx = ChemistryTools._get_fragment_and_attachment(smiles_fragment, variant_idx)
        if not fragment:
            raise ValueError(f"Failed to parse or find fragment: {smiles_fragment}")
            
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

        # 5. Reconnect adjusted boundary connections to the new fragment's attachment index
        if frag_attach_idx != -1:
            new_frag_attach_idx = frag_idx_map[frag_attach_idx]
            for n_idx, b_type in adjusted_connections:
                rw_mol.AddBond(n_idx, new_frag_attach_idx, b_type)

        try:
            rw_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rw_mol)
            AllChem.Compute2DCoords(rw_mol)
        except Exception as e:
            raise ValueError(f"Replacing substructure failed: {str(e)}")

        return rw_mol.GetMol()


    @staticmethod
    def get_mapped_smiles(mol: Chem.Mol) -> str:
        """
        Returns SMILES with atom indices mapped using persistent map numbers.
        """
        # Create a copy to not modify original properties
        m = Chem.Mol(mol)
        for atom in m.GetAtoms():
            # If for some reason map number is missing, fallback to 1-based index
            if atom.GetAtomMapNum() == 0:
                atom.SetAtomMapNum(atom.GetIdx() + 1)
        return Chem.MolToSmiles(m)

    @staticmethod
    def find_substructure(mol: Chem.Mol, pattern_name: str) -> list[list[int]]:
        """
        Find all occurrences of a named substructure in a molecule.
        
        Args:
            mol: RDKit molecule object
            pattern_name: Name of the substructure (e.g., "phenyl", "benzimidazole")
        
        Returns:
            List of atom ID lists, one for each match found.
        """
        from app.chemistry.substructures import find_substructure as _find
        return _find(mol, pattern_name)

    @staticmethod
    def aromatize(mol: Chem.Mol) -> Chem.Mol:
        """
        Convert molecule to aromatic form.
        """
        new_mol = Chem.Mol(mol)
        Chem.SetAromaticity(new_mol)
        return new_mol

    @staticmethod
    def dearomatize(mol: Chem.Mol) -> Chem.Mol:
        """
        Convert molecule to Kekulé form (alternating single/double bonds).
        """
        new_mol = Chem.RWMol(mol)
        Chem.Kekulize(new_mol, clearAromaticFlags=True)
        return new_mol.GetMol()

    @staticmethod
    def calculate_properties(mol: Chem.Mol) -> dict:
        """
        Calculate molecular properties.
        """
        from rdkit.Chem import Descriptors, rdMolDescriptors
        
        return {
            "molecular_weight": round(Descriptors.MolWt(mol), 2),
            "exact_mass": round(Descriptors.ExactMolWt(mol), 4),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "num_rings": rdMolDescriptors.CalcNumRings(mol),
            "num_aromatic_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "num_rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "num_hbd": rdMolDescriptors.CalcNumHBD(mol),
            "num_hba": rdMolDescriptors.CalcNumHBA(mol),
            "tpsa": round(rdMolDescriptors.CalcTPSA(mol), 2),
            "logp": round(Descriptors.MolLogP(mol), 2),
        }


    @staticmethod
    def check_structure(mol: Chem.Mol) -> dict:
        """
        Validate structure for common issues.
        """
        issues = []
        
        # Check valence
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            issues.append(f"Sanitization error: {str(e)}")
        
        # Check for radicals
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                issues.append(f"Atom {atom.GetIdx()} ({atom.GetSymbol()}) has radical electrons")
        
        return {
            "valid": len(issues) == 0,
            "issues": issues
        }
