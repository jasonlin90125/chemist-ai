import os
import json
import re
import ast
import base64
from io import BytesIO
from openai import OpenAI
from app.models import EditRequest, VisualMolecule
from app.chemistry.molecule import VisualMoleculeBuilder, align_and_diff
from app.chemistry.tools import ChemistryTools
from rdkit import Chem
from rdkit.Chem import Draw

# Configure OpenRouter inside function to pick up env vars dynamically

TOOLS_SCHEMA = [
    {
        "type": "function",
        "function": {
            "name": "modify_atom",
            "description": "Change an atom's element (e.g. Carbon to Nitrogen).",
            "parameters": {
                "type": "object",
                "properties": {
                    "atom_id": { "type": "integer", "description": "The ID of the atom to change." },
                    "new_element_symbol": { "type": "string", "description": "Periodic table symbol (e.g., 'N', 'O', 'Cl')." }
                },
                "required": ["atom_id", "new_element_symbol"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "replace_atoms",
            "description": "Replace multiple atoms with a new element.",
            "parameters": {
                "type": "object",
                "properties": {
                    "atom_ids": { "type": "array", "items": { "type": "integer" }, "description": "List of atom IDs to replace." },
                    "new_element_symbol": { "type": "string", "description": "Periodic table symbol." }
                },
                "required": ["atom_ids", "new_element_symbol"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "remove_atoms",
            "description": "Remove atoms from the molecule.",
            "parameters": {
                "type": "object",
                "properties": {
                    "atom_ids": { "type": "array", "items": { "type": "integer" }, "description": "List of atom IDs to remove." }
                },
                "required": ["atom_ids"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "add_substructure",
            "description": "Add a functional group or fragment to a specific anchor atom.",
            "parameters": {
                "type": "object",
                "properties": {
                    "anchor_atom_id": { "type": "integer", "description": "The atom to attach to." },
                    "smiles_fragment": { "type": "string", "description": "Valid RDKit SMILES fragment." }
                },
                "required": ["anchor_atom_id", "smiles_fragment"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "replace_substructure",
            "description": "Replace a set of atoms with a new substructure.",
            "parameters": {
                "type": "object",
                "properties": {
                    "atom_ids": { "type": "array", "items": { "type": "integer" }, "description": "Atoms to be replaced." },
                    "smiles_fragment": { "type": "string", "description": "Valid RDKit SMILES fragment." }
                },
                "required": ["atom_ids", "smiles_fragment"]
            }
        }
    }
]

SYSTEM_PROMPT = """You are Chemist.ai, an expert computational chemist.
You do not generate SMILES directly. You operate on the provided chemical graph by calling tools.
Always prioritize stability and validity.

GUIDELINES:
1. When asked to "make this a pyridine", specific atomic modifications are preferred over rebuilding.
2. For `add_substructure`: Provide valid RDKit SMILES. 
   - Methyl: 'C'
   - Trifluoromethyl: 'C(F)(F)F' (NOT 'CF3' - numbers indicate rings)
   - Do not use abbreviations like 'Ph', 'Me'. Use explicit SMILES.
3. Refer to the provided Mapped SMILES to identify atom indices. The indices in the mapped SMILES `[Element:Index]` correspond to the atom_ids you should use in tool calls.
"""

def generate_molecule_image(mol: Chem.Mol) -> str:
    """Generates a base64 encoded PNG of the molecule with indices."""
    try:
        mol = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol)

        # Add atom indices to labels
        for atom in mol.GetAtoms():
            atom.SetProp('atomLabel', f"{atom.GetSymbol()}{atom.GetIdx()}")

        img = Draw.MolToImage(mol, size=(400, 400))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode("utf-8")
    except Exception as e:
        print(f"Error generating image: {e}")
        return ""

async def process_molecule_edit(request: EditRequest) -> VisualMolecule:
    """
    Orchestrates the AI editing process.
    """
    api_key = os.getenv("OPENROUTER_API_KEY")
    if not api_key:
        print("Warning: OPENROUTER_API_KEY not found.")
    
    client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key=api_key or "sk-or-placeholder",
    )

    # 1. Reconstruct RDKit Mol from JSON
    current_mol = VisualMoleculeBuilder.visual_json_to_mol(request.current_molecule)
    
    # Generate Context
    mapped_smiles = ChemistryTools.get_mapped_smiles(current_mol)
    # image_b64 = generate_molecule_image(current_mol) # Ready for vision models

    # 2. Call AI
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": f"""
        Current Molecule (Mapped SMILES): {mapped_smiles}
        Selection Indices: {request.selected_indices}
        
        User Request: "{request.user_prompt}"
        """}
    ]
    
    try:
        response = client.chat.completions.create(
            model="nvidia/nemotron-3-nano-30b-a3b:free",
            messages=messages,
            tools=TOOLS_SCHEMA,
            tool_choice="auto"
        )
        
        message = response.choices[0].message
        print(f"AI Response: {message}")
        
        # 3. Handle Tool Calls
        new_mol = current_mol
        
        tool_calls = message.tool_calls or []
        
        # Fallback for text-based tool calls (Gemma/Free models)
        if not tool_calls and message.content:
            matches = re.finditer(r"\[(\w+)\((.*?)\)\]", message.content)
            for m in matches:
                fn_name = m.group(1)
                args_str = m.group(2)
                try:
                    dict_str = re.sub(r"(\w+)\s*=", r"'\1':", args_str)
                    args = ast.literal_eval(f"{{{dict_str}}}")

                    class PseudoToolCall:
                        class Function:
                            def __init__(self, n, a):
                                self.name = n
                                self.arguments = json.dumps(a)
                        def __init__(self, n, a):
                            self.function = self.Function(n, a)
                    tool_calls.append(PseudoToolCall(fn_name, args))
                except Exception as e:
                    print(f"Fallback parse error: {e}")

        if tool_calls:
            print(f"Tool calls found: {len(tool_calls)}")
            for tool_call in tool_calls:
                fn_name = tool_call.function.name
                args = json.loads(tool_call.function.arguments)
                
                print(f"Executing {fn_name} with {args}")
                
                if fn_name == "modify_atom":
                    new_mol = ChemistryTools.modify_atom(
                        new_mol, 
                        args["atom_id"], 
                        args["new_element_symbol"]
                    )
                elif fn_name == "replace_atoms":
                    new_mol = ChemistryTools.replace_atoms(
                        new_mol,
                        args["atom_ids"],
                        args["new_element_symbol"]
                    )
                elif fn_name == "remove_atoms":
                    new_mol = ChemistryTools.remove_atoms(
                        new_mol,
                        args["atom_ids"]
                    )
                elif fn_name == "add_substructure":
                    new_mol = ChemistryTools.add_substructure(
                        new_mol, 
                        args["anchor_atom_id"], 
                        args["smiles_fragment"]
                    )
                elif fn_name == "replace_substructure":
                    new_mol = ChemistryTools.replace_substructure(
                        new_mol,
                        args["atom_ids"],
                        args["smiles_fragment"]
                    )
        
        # 4. Alignment & Diff
        aligned_vis = align_and_diff(current_mol, new_mol)
        
        # 5. Apply Diff Metadata
        final_vis = ChemistryTools.apply_diff_metadata(current_mol, new_mol, aligned_vis)
        
        return final_vis

    except Exception as e:
        print(f"AI Error: {e}")
        return request.current_molecule
