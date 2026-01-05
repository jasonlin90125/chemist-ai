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
from rdkit.Chem import Draw, AllChem

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

SYSTEM_PROMPT = """
You are a master organic chemist and molecular editor AI. 
Your goal is to modify chemical structures based on natural language requests.

Available tools:
- modify_atom(atom_id, new_element_symbol)
- replace_atoms(atom_ids, new_element_symbol)
- remove_atoms(atom_ids)
- add_substructure(anchor_atom_id, smiles_fragment)
- replace_substructure(atom_ids, smiles_fragment)

Guidelines:
1. You will be provided with:
    - A Mapped SMILES representation of the current molecule where atoms are noted as [Element:Index].
    - A set of 'Selected Indices' from the user.
    - An image of the molecule with atom indices labeled as 'SymbolINDEX' (e.g., C4, N10).
2. Use the image and mapped SMILES together to confirm which atom_ids correspond to the user's request.
3. Only use the provided tools.
4. Ensure SMILES fragments are chemically valid and use explicit parentheses for branching:
   - Methyl: 'C'
   - Trifluoromethyl: 'C(F)(F)F' (NOT 'CF3' - numbers indicate rings)
   - Do not use abbreviations like 'Ph', 'Me'. Use explicit SMILES.
5. The indices in the mapped SMILES and image correspond to the atom_ids you should use in tool calls.
"""

def generate_molecule_image(mol: Chem.Mol) -> str:
    """Generates a base64 encoded PNG of the molecule with indices."""
    try:
        mol = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol)

        # Add atom indices to labels
        for atom in mol.GetAtoms():
            atom.SetProp('atomLabel', f"{atom.GetSymbol()}{atom.GetIdx()}")

        img = Draw.MolToImage(mol, size=(600, 600))
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
    image_b64 = generate_molecule_image(current_mol)

    # 2. Call AI
    user_content = [
        {
            "type": "text", 
            "text": f"Current Molecule (Mapped SMILES): {mapped_smiles}\nSelection Indices: {request.selected_indices}\n\nUser Request: \"{request.user_prompt}\""
        }
    ]
    
    if image_b64:
        user_content.append({
            "type": "image_url",
            "image_url": {
                "url": f"data:image/png;base64,{image_b64}"
            }
        })

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": user_content}
    ]
    
    try:
        response = client.chat.completions.create(
            model="google/gemma-3-27b-it:free",
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
            # Clean up content for parsing
            clean_content = message.content.strip()
            # Remove markdown code blocks if present
            clean_content = re.sub(r"^```python\s*|\s*```$", "", clean_content, flags=re.MULTILINE)
            # Ensure it looks like a list or a series of calls
            if not clean_content.startswith("["):
                clean_content = f"[{clean_content}]"
            
            try:
                tree = ast.parse(clean_content)
                if tree.body and isinstance(tree.body[0], ast.Expr) and isinstance(tree.body[0].value, ast.List):
                    for node in tree.body[0].value.elts:
                        if isinstance(node, ast.Call):
                            fn_name = node.func.id
                            # Extract keyword arguments
                            args = {kw.arg: ast.literal_eval(kw.value) for kw in node.keywords}
                            
                            class PseudoToolCall:
                                class Function:
                                    def __init__(self, n, a):
                                        self.name = n
                                        self.arguments = json.dumps(a)
                                def __init__(self, n, a):
                                    self.function = self.Function(n, a)
                            tool_calls.append(PseudoToolCall(fn_name, args))
            except Exception as e:
                print(f"Fallback parse error: {e}. Content was: {clean_content}")

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
        
        # Check if any changes were made
        has_changes = any(a.diff_state != "EXISTING" for a in final_vis.atoms) or \
                      any(b.diff_state != "EXISTING" for b in final_vis.bonds)
        
        if not has_changes:
            raise Exception("No changes proposed by AI.")

        return final_vis

    except Exception as e:
        print(f"AI Error: {e}")
        # If it's our direct "No changes" error, re-raise to bubble up
        if "No changes proposed" in str(e):
            raise e
        return request.current_molecule
