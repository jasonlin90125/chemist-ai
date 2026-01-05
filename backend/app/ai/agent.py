import os
import json
import re
import ast
from openai import OpenAI
from app.models import EditRequest, VisualMolecule
from app.chemistry.molecule import VisualMoleculeBuilder, align_and_diff
from app.chemistry.tools import ChemistryTools

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
            "name": "add_substructure",
            "description": "Add a functional group or fragment to a specific anchor atom.",
            "parameters": {
                "type": "object",
                "properties": {
                    "anchor_atom_id": { "type": "integer", "description": "The atom to attach to." },
                    "smiles_fragment": { "type": "string", "description": "Valid RDKit SMILES fragment (e.g., 'C', 'C(F)(F)F', 'N'). Do NOT use 'CF3'." }
                },
                "required": ["anchor_atom_id", "smiles_fragment"]
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
"""

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
    
    # 2. Call AI
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": f"""
        Current Molecule: {len(request.current_molecule.atoms)} atoms.
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
        new_mol = current_mol # Default if no change
        
        tool_calls = message.tool_calls or []
        
        # Fallback for text-based tool calls (Gemma/Free models)
        if not tool_calls and message.content:
            matches = re.finditer(r"\[(\w+)\((.*?)\)\]", message.content)
            for m in matches:
                fn_name = m.group(1)
                args_str = m.group(2)
                try:
                    # Parse kwargs string: "a=1, b='c'" -> dict
                    # Transform keys to quoted keys for JSON-like dictionary: a=1 -> 'a': 1
                    # Regex finds identifiers followed by = and quotes them.
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
                elif fn_name == "add_substructure":
                    new_mol = ChemistryTools.add_substructure(
                        new_mol, 
                        args["anchor_atom_id"], 
                        args["smiles_fragment"]
                    )
        
        # 4. Alignment & Diff
        # We align new_mol to current_mol
        aligned_vis = align_and_diff(current_mol, new_mol)
        
        # 5. Apply Diff Metadata (Logic to mark ADDED vs EXISTING)
        final_vis = ChemistryTools.apply_diff_metadata(current_mol, new_mol, aligned_vis)
        
        return final_vis

    except Exception as e:
        print(f"AI Error: {e}")
        # Return original on error
        return request.current_molecule
