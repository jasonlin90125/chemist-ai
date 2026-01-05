"""
AI Agent for molecule editing.
Simplified for capable LLMs (Gemini 3.0 Flash, Claude 3.5, etc.)
"""

import os
import json
import base64
from io import BytesIO
from openai import OpenAI
from app.models import EditRequest, VisualMolecule
from app.chemistry.molecule import VisualMoleculeBuilder, align_and_diff
from app.chemistry.tools import ChemistryTools
from app.chemistry.substructures import get_available_patterns, search_patterns
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

# ==================== TOOL SCHEMA ====================

TOOLS_SCHEMA = [
    {
        "type": "function",
        "function": {
            "name": "find_substructure",
            "description": "Find all occurrences of a named substructure in the molecule. Use this before replacing a group.",
            "parameters": {
                "type": "object",
                "properties": {
                    "pattern_name": {
                        "type": "string",
                        "description": "Name of the substructure (e.g., 'phenyl', 'benzimidazole', 'pyridine')."
                    }
                },
                "required": ["pattern_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "search_patterns",
            "description": "Search for available substructure pattern names. Use when unsure of exact name.",
            "parameters": {
                "type": "object",
                "properties": {
                    "query": {
                        "type": "string",
                        "description": "Search term (e.g., 'pyrid' to find pyridine, pyrimidine, etc.)."
                    }
                },
                "required": ["query"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "modify_atom",
            "description": "Change an atom's element.",
            "parameters": {
                "type": "object",
                "properties": {
                    "atom_id": {"type": "integer", "description": "Atom ID to modify."},
                    "new_element": {"type": "string", "description": "New element symbol (e.g., 'N', 'O', 'S')."}
                },
                "required": ["atom_id", "new_element"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "add_substructure",
            "description": "Add a fragment to an anchor atom.",
            "parameters": {
                "type": "object",
                "properties": {
                    "anchor_atom_id": {"type": "integer", "description": "Atom ID to attach to."},
                    "fragment": {"type": "string", "description": "Fragment name (e.g., 'phenyl', 'methyl') or SMILES."}
                },
                "required": ["anchor_atom_id", "fragment"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "replace_substructure",
            "description": "Replace atoms with a new fragment.",
            "parameters": {
                "type": "object",
                "properties": {
                    "atom_ids": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "description": "Atom IDs to replace."
                    },
                    "fragment": {"type": "string", "description": "New fragment name or SMILES."}
                },
                "required": ["atom_ids", "fragment"]
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
                    "atom_ids": {
                        "type": "array",
                        "items": {"type": "integer"},
                        "description": "Atom IDs to remove."
                    }
                },
                "required": ["atom_ids"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "aromatize",
            "description": "Convert molecule to aromatic representation.",
            "parameters": {"type": "object", "properties": {}}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "dearomatize",
            "description": "Convert molecule to Kekulé form (explicit single/double bonds).",
            "parameters": {"type": "object", "properties": {}}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "calculate_properties",
            "description": "Calculate molecular properties (MW, formula, LogP, etc.).",
            "parameters": {"type": "object", "properties": {}}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "check_structure",
            "description": "Validate the molecule for chemical errors.",
            "parameters": {"type": "object", "properties": {}}
        }
    }
]

# ==================== SYSTEM PROMPT ====================

SYSTEM_PROMPT = """You are a molecular editor AI. Modify molecules using the provided tools.

## Workflow
1. **To replace a substructure**: First use `find_substructure` to get atom IDs, then use `replace_substructure`.
2. **To add a group**: Use `add_substructure` with the anchor atom ID (from user selection or image).
3. **To modify atoms**: Use `modify_atom` or `remove_atoms`.

## Context Provided
- **Mapped SMILES**: Atoms labeled as [Element:ID] (e.g., [C:5], [N:12]).
- **Image**: Atoms labeled as "SymbolID" (e.g., C5, N12).
- **Selected Indices**: Atoms the user has selected in the editor.

## Fragment Names
Common fragments: phenyl, pyridine, pyrimidine, furan, thiophene, imidazole, benzimidazole, 
piperidine, piperazine, morpholine, trifluoromethyl, methoxy, cyano, isopropyl, tert_butyl.

Use `search_patterns` if unsure of exact name.

## Rules
- Use atom IDs from the SMILES mapping or image labels.
- For raw SMILES fragments, use proper notation: C(F)(F)F for trifluoromethyl, c1ccccc1 for benzene.
"""


def generate_molecule_image(mol: Chem.Mol) -> str:
    """Generates a base64 encoded PNG of the molecule with atom indices."""
    try:
        mol = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol)

        for atom in mol.GetAtoms():
            label = f"{atom.GetSymbol()}{atom.GetIdx()}"
            atom.SetProp("atomNote", label)

        drawer = Draw.MolDraw2DCairo(600, 500)
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        img_bytes = drawer.GetDrawingText()
        return base64.b64encode(img_bytes).decode("utf-8")
    except Exception as e:
        print(f"Image generation failed: {e}")
        return ""


async def process_molecule_edit(request: EditRequest) -> VisualMolecule:
    """
    Main entry point for AI-powered molecule editing.
    """
    client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key=os.getenv("OPENROUTER_API_KEY")
    )

    # Get current molecule from request
    mol_block = request.current_molecule.mol_block
    if not mol_block:
        raise ValueError("No mol_block provided in request")

    current_mol = Chem.MolFromMolBlock(mol_block)
    if not current_mol:
        raise ValueError("Failed to parse mol_block")

    # Generate context
    mapped_smiles = ChemistryTools.get_mapped_smiles(current_mol)
    mol_image = generate_molecule_image(current_mol)

    # Build messages
    user_content = [
        {"type": "text", "text": f"**Mapped SMILES**: {mapped_smiles}"},
        {"type": "text", "text": f"**Selected Indices**: {request.selected_indices}"},
        {"type": "text", "text": f"**User Request**: {request.user_prompt}"},
    ]

    if mol_image:
        user_content.append({
            "type": "image_url",
            "image_url": {"url": f"data:image/png;base64,{mol_image}"}
        })

    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": user_content}
    ]

    # Run agent loop
    new_mol = current_mol
    has_modified = False
    max_turns = 5

    for turn in range(max_turns):
        response = client.chat.completions.create(
            model="google/gemini-3-flash-preview",
            messages=messages,
            tools=TOOLS_SCHEMA,
            tool_choice="auto",
            max_tokens=2048  # Limit to prevent credit exhaustion
        )

        ai_msg = response.choices[0].message
        messages.append(ai_msg)

        tool_calls = ai_msg.tool_calls or []
        
        if not tool_calls:
            # No more tool calls, exit loop
            break

        print(f"Turn {turn}: {len(tool_calls)} tool call(s)")

        for tool_call in tool_calls:
            fn_name = tool_call.function.name
            args = json.loads(tool_call.function.arguments)
            result = ""

            try:
                if fn_name == "find_substructure":
                    matches = ChemistryTools.find_substructure(new_mol, args["pattern_name"])
                    result = json.dumps({"matches": matches, "count": len(matches)})

                elif fn_name == "search_patterns":
                    patterns = search_patterns(args["query"])
                    result = json.dumps({"patterns": patterns})

                elif fn_name == "modify_atom":
                    new_mol = ChemistryTools.modify_atom(new_mol, args["atom_id"], args["new_element"])
                    has_modified = True
                    result = "Success: Atom modified."

                elif fn_name == "add_substructure":
                    new_mol = ChemistryTools.add_substructure(new_mol, args["anchor_atom_id"], args["fragment"])
                    has_modified = True
                    result = "Success: Fragment added."

                elif fn_name == "replace_substructure":
                    new_mol = ChemistryTools.replace_substructure(new_mol, args["atom_ids"], args["fragment"])
                    has_modified = True
                    result = "Success: Substructure replaced."

                elif fn_name == "remove_atoms":
                    new_mol = ChemistryTools.remove_atoms(new_mol, args["atom_ids"])
                    has_modified = True
                    result = "Success: Atoms removed."

                elif fn_name == "aromatize":
                    new_mol = ChemistryTools.aromatize(new_mol)
                    has_modified = True
                    result = "Success: Aromatized."

                elif fn_name == "dearomatize":
                    new_mol = ChemistryTools.dearomatize(new_mol)
                    has_modified = True
                    result = "Success: Dearomatized."

                elif fn_name == "calculate_properties":
                    props = ChemistryTools.calculate_properties(new_mol)
                    result = json.dumps(props)

                elif fn_name == "check_structure":
                    check = ChemistryTools.check_structure(new_mol)
                    result = json.dumps(check)

                else:
                    result = f"Error: Unknown function {fn_name}"

            except Exception as e:
                result = f"Error: {str(e)}"
                print(f"Tool {fn_name} failed: {e}")

            messages.append({
                "role": "tool",
                "tool_call_id": tool_call.id,
                "name": fn_name,
                "content": result
            })

    # Finalize
    if not has_modified:
        raise Exception("No changes proposed by AI.")

    aligned_vis = align_and_diff(current_mol, new_mol)
    final_vis = ChemistryTools.apply_diff_metadata(current_mol, new_mol, aligned_vis)

    # Verify changes
    has_changes = any(a.diff_state != "EXISTING" for a in final_vis.atoms) or \
                  any(b.diff_state != "EXISTING" for b in final_vis.bonds)

    if not has_changes:
        raise Exception("No visible changes in molecule.")

    return final_vis
