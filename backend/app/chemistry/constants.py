FUNCTIONAL_GROUPS = {
    # --- Common Carbon Frameworks ---
    "phenyl": "c1ccccc1*",
    "benzyl": "*Cc1ccccc1",
    "isopropyl": "*C(C)C",
    "tert-butyl": "*C(C)(C)C",
    "cyclopropyl": "*C1CC1",
    "cyclobutyl": "*C1CCC1",
    "cyclopentyl": "*C1CCCC1",
    "cyclohexyl": "*C1CCCCC1",
    
    # --- Common Heterocycles (Aromatic) ---
    "pyridine_2_yl": "*c1ccccn1",
    "pyridine_3_yl": "*c1cnccc1",
    "pyridine_4_yl": "*c1ccncc1",
    "pyridine": "*c1cnccc1",
    "pyrimidine": "*c1ncccn1",
    "pyrazine": "*c1cnccn1",
    "pyridazine": "*c1ccnnc1",
    "furan": "*c1occc1",
    "thiophene": "*c1sccc1",
    "oxazole": "*c1ncno1",
    "isoxazole": "*c1nocc1",
    "thiazole": "*c1ncns1",
    "isothiazole": "*c1nscc1",
    "pyrazole": "*c1n[nH]cc1",
    "imidazole": "*c1nc[nH]c1",
    "triazole_123": "*c1nn[nH]c1",
    "triazole_124": "*c1ncn[nH]1",
    "tetrazole": "*c1nnn[nH]1",
    "indole": "*c1ccc2[nH]ccc2c1",
    "benzimidazole": "*c1ccc2[nH]cnc2c1",
    "quinoline": "*c1cccc2ncccc12",
    "isoquinoline": "*c1cccc2cnccc12",
    
    # --- Saturated Heterocycles ---
    "piperidine": "*N1CCCCC1",
    "morpholine": "*N1CCOCC1",
    "piperazine": "*N1CCNCC1",
    "pyrrolidine": "*N1CCCC1",
    "tetrahydrofuran": "*C1CCCO1",
    "tetrahydropyran": "*C1CCCCO1",
    "azetidine": "*N1CCC1",
    "oxetane": "*C1COC1",  # Common gem-dimethyl bioisostere
    
    # --- Polar / Reactive Groups ---
    "hydroxyl": "*O",
    "ether": "*OC",
    "amine_primary": "*N",
    "amine_dimethyl": "*N(C)C",
    "carboxylic_acid": "*C(=O)O",
    "ester_methyl": "*C(=O)OC",
    "amide_primary": "*C(=O)N",
    "amide_methyl": "*C(=O)NC",
    "cyano": "*C#N",
    "nitro": "*[N+](=O)[O-]",
    "fluorine": "*F",
    "chlorine": "*Cl",
    "bromine": "*Br",
    "trifluoromethyl": "*C(F)(F)F",
    "trifluoromethoxy": "*OC(F)(F)F",
    
    # --- Sulfur / Phosphorus ---
    "thiol": "*S",
    "thioether": "*SC",
    "sulfone": "*S(=O)(=O)C",
    "sulfonamide": "*S(=O)(=O)N",
    "phosphate": "*OP(=O)(O)O",
    "phosphonate": "*P(=O)(O)O",
    
    # --- Modern Bioisosteres (Advanced) ---
    "bicyclo_111_pentane": "*C12CC(C1)C2",  # Phenyl bioisostere
    "spiro_33_heptane": "*C12CCC(C1)C2",
    "cubane": "*C12C3C4C1C5C2C3C45",
}
