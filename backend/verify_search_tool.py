import asyncio
import os
from dotenv import load_dotenv
from app.ai.agent import process_molecule_edit
from app.models import EditRequest, VisualMolecule
from app.chemistry.molecule import get_ibrutinib

load_dotenv()

async def test_search_and_edit():
    print("Testing Multi-turn Tool calling (Search -> Edit)...")
    mol = get_ibrutinib()
    request = EditRequest(
        current_molecule=mol,
        user_prompt="I need to add a 'phenyl' group to atom 32. Use the search tool to find the name 'phenyl' first to confirm it exists.",
        selected_indices=[32]
    )
    
    try:
        # We can't see the internal prints easily, so we rely on the final result
        result = await process_molecule_edit(request)
        print(f"Final atom count: {len(result.atoms)} (Original: {mol.GetNumAtoms()})")
        if len(result.atoms) > mol.GetNumAtoms():
            print("SUCCESS: Molecule was modified after tool search.")
    except Exception as e:
        print(f"Caught error: {e}")

if __name__ == "__main__":
    asyncio.run(test_search_and_edit())
