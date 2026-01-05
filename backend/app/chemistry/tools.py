from rdkit import Chem
from rdkit.Chem import AllChem
from app.chemistry.molecule import VisualMoleculeBuilder, align_and_diff
from app.chemistry.constants import FUNCTIONAL_GROUPS
from app.models import VisualMolecule
import copy

class ChemistryTools:
    @staticmethod
    def _get_fragment_and_attachment(smiles_or_name: str) -> tuple[Chem.Mol, int]:
        """
        Helper to get RDKit mol and the index of the atom to attach.
        If SMILES has '*', it uses the neighbor of '*' as attachment point.
        Returns (mol_without_dummy, attachment_idx).
        """
        # 1. Lookup in dictionary
        smiles = FUNCTIONAL_GROUPS.get(smiles_or_name.lower(), smiles_or_name)
        
        # 2. Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return None, -1
            
        # 3. Find dummy atom '*' (atomic number 0)
        attachment_idx = 0 # Default
        dummy_idx = -1
        
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                dummy_idx = atom.GetIdx()
                # Find the neighbor of this dummy atom
                neighbors = atom.GetNeighbors()
                if neighbors:
                    attachment_idx = neighbors[0].GetIdx()
                break
        
        # 4. Remove dummy atom if found
        if dummy_idx != -1:
            nm = Chem.RWMol(mol)
            nm.RemoveAtom(dummy_idx)
            # Adjust attachment_idx if it was shifted
            if dummy_idx < attachment_idx:
                attachment_idx -= 1
            return nm.GetMol(), attachment_idx
            
        return mol, attachment_idx

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
    def add_substructure(current_mol: Chem.Mol, anchor_atom_id: int, smiles_fragment: str) -> Chem.Mol:
        """
        Adds a fragment at a specific anchor point.
        Automatically replaces an implicit hydrogen if needed.
        """
        rw_mol = Chem.RWMol(current_mol)
        
        if anchor_atom_id >= rw_mol.GetNumAtoms():
             raise ValueError(f"Anchor atom ID {anchor_atom_id} is out of range.")

        fragment, frag_attach_idx = ChemistryTools._get_fragment_and_attachment(smiles_fragment)
        
        if not fragment:
            raise ValueError(f"Failed to parse or find fragment: {smiles_fragment}")
        
        # Check if anchor atom has implicit hydrogens we can replace
        anchor_atom = rw_mol.GetAtomWithIdx(anchor_atom_id)
        
        # Update property cache to get accurate H count
        rw_mol.UpdatePropertyCache(strict=False)
        implicit_hs = anchor_atom.GetNumImplicitHs()
        
        if implicit_hs == 0:
            # Check explicit Hs
            explicit_h_neighbors = [n for n in anchor_atom.GetNeighbors() if n.GetAtomicNum() == 1]
            if explicit_h_neighbors:
                # Remove one explicit H to make room
                h_to_remove = explicit_h_neighbors[0].GetIdx()
                rw_mol.RemoveAtom(h_to_remove)
                # Adjust anchor_atom_id if needed
                if h_to_remove < anchor_atom_id:
                    anchor_atom_id -= 1
            # If no Hs available, we'll try anyway and let sanitization catch errors
            
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
            
            # Attempt to generate coords ONLY for the new atoms, preserving the old ones.
            # But 'add_substructure' might have messed up indices or connected things weirdly.
            # The standard Compute2DCoords resets everything.
            # Instead, we will try to use GenerateDepictionMatching2DStructure if possible,
            # BUT we need a reference. 'current_mol' has the original coords.

            try:
                # We want to preserve the layout of 'current_mol' as much as possible.
                # GenerateDepictionMatching2DStructure works best when the template has standard bond lengths (~1.5A).
                # If 'current_mol' is scaled (e.g. from UI), the new atoms will be generated at standard size,
                # creating a distortion. We need to normalize scale first.

                if rw_mol.GetNumConformers() == 0:
                     AllChem.Compute2DCoords(rw_mol)
                else:
                    if current_mol.GetNumConformers() > 0:
                        # 1. Calculate scale factor
                        total_len = 0.0
                        count = 0
                        conf = current_mol.GetConformer()
                        for bond in current_mol.GetBonds():
                            u = bond.GetBeginAtomIdx()
                            v = bond.GetEndAtomIdx()
                            p1 = conf.GetAtomPosition(u)
                            p2 = conf.GetAtomPosition(v)
                            dist = ((p1.x - p2.x)**2 + (p1.y - p2.y)**2)**0.5
                            total_len += dist
                            count += 1

                        avg_len = total_len / count if count > 0 else 1.5
                        scale_factor = avg_len / 1.5

                        # 2. Normalize template if needed
                        # Threshold: if scale is off by more than 10%
                        if abs(scale_factor - 1.0) > 0.1:
                            template_mol = Chem.Mol(current_mol)
                            template_conf = template_mol.GetConformer()
                            norm_scale = 1.0 / scale_factor
                            for i in range(template_mol.GetNumAtoms()):
                                pos = template_conf.GetAtomPosition(i)
                                template_conf.SetAtomPosition(i, (pos.x * norm_scale, pos.y * norm_scale, pos.z * norm_scale))

                            # Generate coords using normalized template.
                            # The new atoms will be generated at standard size (~1.5), which matches the template now.
                            AllChem.GenerateDepictionMatching2DStructure(rw_mol, template_mol, acceptFailure=True)

                            # 3. Scale result back to original size
                            rw_conf = rw_mol.GetConformer()
                            for i in range(rw_mol.GetNumAtoms()):
                                pos = rw_conf.GetAtomPosition(i)
                                rw_conf.SetAtomPosition(i, (pos.x * scale_factor, pos.y * scale_factor, pos.z * scale_factor))
                        else:
                            # Standard scale, proceed as usual
                            AllChem.GenerateDepictionMatching2DStructure(rw_mol, current_mol, acceptFailure=True)
                    else:
                         AllChem.Compute2DCoords(rw_mol)

            except Exception:
                # Fallback
                AllChem.Compute2DCoords(rw_mol)
            
        except Exception as e:
            raise ValueError(f"Adding substructure failed: {str(e)}")

        return rw_mol.GetMol()

    @staticmethod
    def replace_substructure(current_mol: Chem.Mol, atom_ids: list[int], smiles_fragment: str) -> Chem.Mol:
        """
        Removes the specified atoms and attaches the new fragment to the neighbor of the first removed atom.
        """
        rw_mol = Chem.RWMol(current_mol)
        atoms_to_remove = set(atom_ids)
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
            if atom_id >= rw_mol.GetNumAtoms():
                 raise ValueError(f"Atom ID {atom_id} is out of range.")
            rw_mol.RemoveAtom(atom_id)

        # Adjust attachment_point index
        current_attachment = attachment_point
        if current_attachment != -1:
            for atom_id in sorted_ids:
                if atom_id < current_attachment:
                    current_attachment -= 1
                elif atom_id == current_attachment:
                    current_attachment = -1
                    break

        # Add fragment
        if current_attachment != -1:
            fragment, frag_attach_idx = ChemistryTools._get_fragment_and_attachment(smiles_fragment)
            if fragment and frag_attach_idx != -1:
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
                    frag_idx_map[frag_attach_idx],
                    Chem.BondType.SINGLE
                )
            elif not fragment:
                raise ValueError(f"Failed to parse or find fragment: {smiles_fragment}")

        try:
            rw_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rw_mol)
        except Exception as e:
            raise ValueError(f"Replacing substructure failed: {str(e)}")

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

    # ==================== NEW TOOLS ====================

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
