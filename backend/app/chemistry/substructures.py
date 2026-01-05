"""
SMARTS patterns for common substructure detection.
Used by find_substructure() to locate named groups in molecules.
"""

from rdkit import Chem

# SMARTS patterns for common ring systems and functional groups
SUBSTRUCTURE_PATTERNS = {
    # --- Aromatic Rings ---
    "phenyl": "c1ccccc1",
    "benzene": "c1ccccc1",
    "pyridine": "c1ccncc1",
    "pyrimidine": "c1ncncc1",
    "pyrazine": "c1cnccn1",
    "pyridazine": "c1ccnnc1",
    "furan": "c1ccoc1",
    "thiophene": "c1ccsc1",
    "pyrrole": "c1cc[nH]c1",
    "imidazole": "c1cnc[nH]1",
    "pyrazole": "c1cc[nH]n1",
    "oxazole": "o1cncc1",
    "isoxazole": "o1nccc1",
    "thiazole": "s1cncc1",
    "isothiazole": "s1nccc1",
    "triazole": "c1ncnn1",
    "1,2,3-triazole": "c1nnn[nH]1",
    "1,2,4-triazole": "c1nc[nH]n1",
    "tetrazole": "c1nnn[nH]1",
    "1,2,4-oxadiazole": "n1conn1",
    "1,3,4-oxadiazole": "n1cnno1",
    "1,3,5-triazine": "c1ncncn1",
    "1,2,4-triazine": "c1nnccn1",

    # --- Fused Rings ---
    "naphthalene": "c1ccc2ccccc2c1",
    "indole": "c1ccc2[nH]ccc2c1",
    "indazole": "c1ccc2[nH]ncc2c1",
    "benzimidazole": "c1ccc2[nH]cnc2c1",
    "benzofuran": "c1ccc2occc2c1",
    "benzothiophene": "c1ccc2sccc2c1",
    "benzothiazole": "c1ccc2scnc2c1",
    "benzoxazole": "c1ccc2ocnc2c1",
    "quinoline": "c1ccc2ncccc2c1",
    "isoquinoline": "c1ccc2cnnc2c1",
    "quinazoline": "c1ccc2ncncc2c1",
    "phthalazine": "c1ccc2cnnc2c1",
    "coumarin": "O=C1C=Cc2ccccc2O1",
    "chromone": "O=C1C=Cc2ccccc2O1",
    "purine": "c1nc2[nH]cnc2n1",
    "adenine": "Nc1ncnc2[nH]cnc12",

    # --- Saturated Heterocycles & Cages ---
    "epoxide": "C1OC1",
    "oxetane": "C1COC1",
    "azetidine": "C1CNC1",
    "tetrahydrofuran": "C1CCOC1",
    "pyrrolidine": "C1CCNC1",
    "piperidine": "C1CCNCC1",
    "piperazine": "C1CNCCN1",
    "morpholine": "C1COCCN1",
    "tetrahydropyran": "C1CCOCC1",
    "1,4-dioxane": "C1COCCO1",
    "adamantane": "C1C2CC3CC1CC(C2)C3",
    "norbornane": "C1CC2CCC1C2",

    # --- Functional Groups ---
    "carboxylic_acid": "C(=O)O",
    "ester": "C(=O)OC",
    "amide": "C(=O)N",
    "urea": "NC(=O)N",
    "carbamate": "NC(=O)O",
    "ketone": "CC(=O)C",
    "aldehyde": "[CH]=O",
    "alcohol": "[OH]",
    "ether": "COC",
    "amine": "[NH2]",
    "imine": "C=N",
    "hydrazine": "NN",
    "guanidine": "NC(=N)N",
    "amidine": "C(=N)N",
    "nitro": "[N+](=O)[O-]",
    
    # --- Sulfur/Phosphorus Groups ---
    "thiol": "[SH]",
    "thioether": "CSC",
    "disulfide": "SS",
    "sulfonamide": "S(=O)(=O)N",
    "sulfone": "S(=O)(=O)",
    "sulfonyl_chloride": "S(=O)(=O)Cl",
    "phosphate": "P(=O)(O)(O)",
    "phosphonate": "P(=O)(O)C",
    "thiourea": "NC(=S)N",

    # --- Common Substituents ---
    "trifluoromethyl": "C(F)(F)F",
    "methoxy": "CO",
    "acetyl": "CC(=O)",
    "cyano": "C#N",
    "isopropyl": "C(C)C",
    "tert_butyl": "C(C)(C)C",
    "cyclopropyl": "C1CC1",
    "cyclobutyl": "C1CCC1",
    "cyclopentyl": "C1CCCC1",
    "cyclohexyl": "C1CCCCC1",
}


def find_substructure(mol: Chem.Mol, pattern_name: str) -> list[list[int]]:
    """
    Find all occurrences of a named substructure in a molecule.
    
    Args:
        mol: RDKit molecule object
        pattern_name: Name of the substructure (e.g., "phenyl", "benzimidazole")
    
    Returns:
        List of atom ID lists, one for each match found.
        Empty list if pattern not found or not recognized.
    """
    pattern_name = pattern_name.lower().replace(" ", "_").replace("-", "_")
    
    smarts = SUBSTRUCTURE_PATTERNS.get(pattern_name)
    if not smarts:
        return []
    
    pattern = Chem.MolFromSmarts(smarts)
    if not pattern:
        return []
    
    matches = mol.GetSubstructMatches(pattern)
    return [list(match) for match in matches]


def get_available_patterns() -> list[str]:
    """Return list of all available substructure pattern names."""
    return list(SUBSTRUCTURE_PATTERNS.keys())


def search_patterns(query: str) -> list[str]:
    """Search for pattern names matching a query string."""
    query = query.lower()
    return [name for name in SUBSTRUCTURE_PATTERNS.keys() if query in name]
