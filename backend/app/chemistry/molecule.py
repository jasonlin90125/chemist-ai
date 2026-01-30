from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from app.models import VisualMolecule, Atom, Bond
import uuid

class VisualMoleculeBuilder:
    @staticmethod
    def mol_to_visual_json(mol: Chem.Mol, mol_id: str = None, diff_map: dict[int, str] = None) -> VisualMolecule:
        if not mol_id:
            mol_id = f"mol_{uuid.uuid4().hex[:8]}"
            
        atoms = []
        conf = mol.GetConformer()
        
        # 1. Collect existing map numbers to avoid collisions
        existing_maps = {atom.GetAtomMapNum() for atom in mol.GetAtoms() if atom.GetAtomMapNum() > 0}
        next_map = max(existing_maps) + 1 if existing_maps else 1

        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            
            # Ensure we have a persistent atom map number
            map_num = atom.GetAtomMapNum()
            if map_num == 0:
                while next_map in existing_maps:
                    next_map += 1
                map_num = next_map
                atom.SetAtomMapNum(map_num)
                existing_maps.add(map_num)
                next_map += 1

            # Use diff_map if provided, else default to EXISTING
            diff_state = diff_map.get(idx, "EXISTING") if diff_map else "EXISTING"
            
            atoms.append(Atom(
                id=idx,
                element=atom.GetSymbol(),
                x=pos.x,
                y=pos.y,
                charge=atom.GetFormalCharge(),
                implicit_h=atom.GetTotalNumHs(),
                atom_map=map_num,
                ui_state="DEFAULT",
                diff_state=diff_state
            ))
            
        bonds = []
        for bond in mol.GetBonds():
            b1, b2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            order_map = {
                Chem.BondType.SINGLE: "SINGLE",
                Chem.BondType.DOUBLE: "DOUBLE",
                Chem.BondType.TRIPLE: "TRIPLE",
                Chem.BondType.AROMATIC: "AROMATIC"
            }
            
            # A bond is EXISTING if both atoms are EXISTING
            bond_diff = "EXISTING"
            if diff_map:
                if diff_map.get(b1) == "ADDED" or diff_map.get(b2) == "ADDED":
                    bond_diff = "ADDED"

            bonds.append(Bond(
                source=b1,
                target=b2,
                order=order_map.get(bond.GetBondType(), "SINGLE"),
                diff_state=bond_diff
            ))
            
        # 3. Create a version for the UI that keeps map numbers (needed for targeting reliability)
        try:
            # Use V3000 to ensure we don't truncate high indices and keep coords exact.
            mol_block = Chem.MolToMolBlock(mol, forceV3000=True)
        except Exception as e:
            print(f"Error generating mol_block: {e}")
            mol_block = None

        # Generate SVG for visual preview
        try:
            # Create a copy and STRIP map numbers for user display
            mc = Chem.Mol(mol)
            for a in mc.GetAtoms():
                a.SetAtomMapNum(0)
                
            AllChem.Compute2DCoords(mc)
            drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)
            options = drawer.drawOptions()
            options.padding = 0.1
            options.addAtomIndices = False
            drawer.DrawMolecule(mc)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
        except:
            svg = None
            
        return VisualMolecule(
            molecule_id=mol_id,
            atoms=atoms,
            bonds=bonds,
            mol_block=mol_block,
            smiles=Chem.MolToSmiles(mol),
            svg=svg
        )

    @staticmethod
    def reconstruct_molecule(visual_mol: VisualMolecule) -> Chem.Mol:
        """
        Reconstructs an RDKit molecule from the editor's MolBlock,
        then uses geometric proximity to restore persistent map numbers from the 
        VisualMolecule data. This is robust against index shifts caused by manual edits.
        """
        if not visual_mol.mol_block:
             return VisualMoleculeBuilder.visual_json_to_mol(visual_mol)
             
        mol = Chem.MolFromMolBlock(visual_mol.mol_block, removeHs=False)
        if not mol:
             return VisualMoleculeBuilder.visual_json_to_mol(visual_mol)

        # IMPORTANT: Check for degenerate coordinates (common after paste)
        # If all atoms are at (0,0), we need to compute a layout
        conf = mol.GetConformer()
        is_degenerate = True
        for i in range(mol.GetNumAtoms()):
            pos = conf.GetAtomPosition(i)
            if abs(pos.x) > 0.01 or abs(pos.y) > 0.01:
                is_degenerate = False
                break
        
        if is_degenerate:
            # Clear existing degenerate conformer and compute fresh ones
            AllChem.Compute2DCoords(mol)
            # Fresh conformer after layout
            conf = mol.GetConformer()

        # Build a lookup of (rounded_x, rounded_y) -> map_num from the persistent data
        # We use a small epsilon for coordinate matching
        coord_map = {}
        for atom_data in visual_mol.atoms:
            if atom_data.atom_map:
                # Round to 3 decimal places to handle floating point noise
                key = (round(atom_data.x, 3), round(atom_data.y, 3))
                coord_map[key] = atom_data.atom_map

        # Apply maps to the new RDKit molecule founded on geometry
        # Re-fetch conf in case it was recomputed
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            # If map is already there (e.g. from V3000 AAM), keep it but still check for consistency
            if atom.GetAtomMapNum() > 0:
                 continue
                 
            # Fallback: geometric match for maps lost during manual edits or reconstruction
            pos = conf.GetAtomPosition(atom.GetIdx())
            key = (round(pos.x, 3), round(pos.y, 3))
            
            if key in coord_map:
                atom.SetAtomMapNum(coord_map[key])
            else:
                best_map = None
                best_dist = 0.015
                for (cx, cy), m_num in coord_map.items():
                    dist = ((pos.x - cx)**2 + (pos.y - cy)**2)**0.5
                    if dist < best_dist:
                        best_dist = dist
                        best_map = m_num
                if best_map:
                    atom.SetAtomMapNum(best_map)

        return mol

    @staticmethod
    def visual_json_to_mol(visual_mol: VisualMolecule) -> Chem.Mol:
        # Create an editable molecule
        rw_mol = Chem.RWMol()
        
        # Add atoms
        atom_map = {} # visual_id -> rdkit_idx
        
        sorted_atoms = sorted(visual_mol.atoms, key=lambda a: a.id)
        
        for atom_data in sorted_atoms:
            a = Chem.Atom(atom_data.element)
            a.SetFormalCharge(atom_data.charge)
            if atom_data.atom_map:
                a.SetAtomMapNum(atom_data.atom_map)
            idx = rw_mol.AddAtom(a)
            atom_map[atom_data.id] = idx
            
        # Add coordinates to conformer
        conf = Chem.Conformer(len(sorted_atoms))
        for atom_data in sorted_atoms:
            conf.SetAtomPosition(atom_map[atom_data.id], (atom_data.x, atom_data.y, 0.0))
        rw_mol.AddConformer(conf)
            
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

def rigid_align(original_mol: Chem.Mol, new_mol: Chem.Mol) -> bool:
    """
    Manually aligns new_mol to original_mol by measuring scale and centroid 
    of their Maximum Common Substructure (MCS).
    Also attempts basic 2D rotation alignment.
    """
    mcs = rdFMCS.FindMCS([original_mol, new_mol], timeout=2)
    if mcs.numAtoms < 2:
        return False

    pattern = Chem.MolFromSmarts(mcs.smartsString)
    match_orig = original_mol.GetSubstructMatch(pattern)
    match_new = new_mol.GetSubstructMatch(pattern)

    if not match_orig or not match_new:
        return False

    conf_orig = original_mol.GetConformer()
    conf_new = new_mol.GetConformer()

    # 1. Centered coordinates
    def get_centroid(mol, indices):
        c = mol.GetConformer()
        x, y = 0.0, 0.0
        for i in indices:
            p = c.GetAtomPosition(i)
            x += p.x ; y += p.y
        return x/len(indices), y/len(indices)

    cx_o, cy_o = get_centroid(original_mol, match_orig)
    cx_n, cy_n = get_centroid(new_mol, match_new)

    # 2. Scale & Rotation (Simplified Procrustes/Kabsch in 2D)
    # Solve for s, theta:  r_orig = s * R(theta) * r_new
    # We use complex numbers approach: (x+iy)_orig = (s*e^itheta) * (x+iy)_new
    
    num = 0j
    den = 0j
    for i in range(len(match_orig)):
        po = conf_orig.GetAtomPosition(match_orig[i])
        pn = conf_new.GetAtomPosition(match_new[i])
        
        zo = complex(po.x - cx_o, po.y - cy_o)
        zn = complex(pn.x - cx_n, pn.y - cy_n)
        
        num += zo * zn.conjugate()
        den += zn * zn.conjugate()

    if abs(den) < 1e-6:
        return False

    # Transformation T = num / den
    t_complex = num / den
    scale = abs(t_complex)
    angle = 0 # Default if scale is 0
    import math
    if scale > 1e-6:
        angle = math.atan2(t_complex.imag, t_complex.real)

    print(f"DEBUG: RigidAlign Scale={scale:.2f}, Angle={math.degrees(angle):.1f}°")

    # 3. Apply Transformation
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    
    for i in range(new_mol.GetNumAtoms()):
        p = conf_new.GetAtomPosition(i)
        # 1. Translate to origin (relative to new centroid)
        nx, ny = p.x - cx_n, p.y - cy_n
        # 2. Scale and Rotate
        rx = scale * (nx * cos_a - ny * sin_a)
        ry = scale * (nx * sin_a + ny * cos_a)
        # 3. Translate to original centroid
        conf_new.SetAtomPosition(i, (rx + cx_o, ry + cy_o, 0.0))

    return True

def align_and_diff(original_mol: Chem.Mol, new_mol: Chem.Mol) -> VisualMolecule:
    # 1. Alignment
    success = False
    try:
        success = rigid_align(original_mol, new_mol)
    except Exception as e:
        print(f"DEBUG: Rigid alignment failed: {e}")

    # 2. Diffing using MCS (Graph-based)
    # We find the mapping of new_mol atoms back to original_mol
    diff_map = {}
    for atom in new_mol.GetAtoms():
        diff_map[atom.GetIdx()] = "ADDED"

    mcs = rdFMCS.FindMCS([original_mol, new_mol], timeout=2)
    if mcs.numAtoms > 0:
        pattern = Chem.MolFromSmarts(mcs.smartsString)
        # Find which atoms in new_mol correspond to the original skeleton
        match_new = new_mol.GetSubstructMatch(pattern)
        for idx in match_new:
            diff_map[idx] = "EXISTING"

    # Final visual molecule with correct diff_state
    vis_mol = VisualMoleculeBuilder.mol_to_visual_json(new_mol, diff_map=diff_map)
    return vis_mol
