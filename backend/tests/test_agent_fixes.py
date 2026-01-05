import pytest
from rdkit import Chem
from app.ai.agent import process_molecule_edit
from app.models import EditRequest, VisualMolecule, Atom
from app.chemistry.molecule import VisualMoleculeBuilder
from unittest.mock import MagicMock, patch

# Helper to create a basic molecule
def create_benzene_visual():
    mol = Chem.MolFromSmiles("c1ccccc1")
    from rdkit.Chem import AllChem
    AllChem.Compute2DCoords(mol)
    Chem.Kekulize(mol)
    return VisualMoleculeBuilder.mol_to_visual_json(mol)

def create_mock_tool_call(name, args):
    # Correctly mocking the structure of tool_call.function.name
    # When agent access tool_call.function.name, it should get the string
    function = MagicMock()
    function.name = name
    function.arguments = args

    tool_call = MagicMock()
    tool_call.function = function
    return tool_call

@pytest.mark.asyncio
async def test_process_molecule_edit_selection():
    # Setup
    vis_mol = create_benzene_visual()
    request = EditRequest(
        current_molecule=vis_mol,
        user_prompt="Select the benzene ring",
        selected_indices=[]
    )

    with patch("app.ai.agent.OpenAI") as mock_openai, \
         patch("app.ai.agent.ChemistryTools") as mock_tools:

        # Mock OpenAI response
        tool_call = create_mock_tool_call("find_substructure", '{"pattern_name": "benzene"}')

        mock_client = mock_openai.return_value
        mock_client.chat.completions.create.side_effect = [
            MagicMock(choices=[MagicMock(message=MagicMock(tool_calls=[tool_call]))]),
            MagicMock(choices=[MagicMock(message=MagicMock(tool_calls=[]))])
        ]

        # Mock ChemistryTools.find_substructure to return indices
        # Benzene indices 0-5
        mock_tools.find_substructure.return_value = [[0, 1, 2, 3, 4, 5]]

        # Mock other needed methods
        mock_tools.get_mapped_smiles.return_value = "c1ccccc1"
        # We need apply_diff_metadata to just return the mol
        mock_tools.apply_diff_metadata.side_effect = lambda o, n, v: v

        # Execute
        result_mol = await process_molecule_edit(request)

        # Verify
        assert isinstance(result_mol, VisualMolecule)
        # Should have atoms 0-5 selected
        selected_count = sum(1 for a in result_mol.atoms if a.ui_state == "SELECTED")
        assert selected_count == 6

@pytest.mark.asyncio
async def test_process_molecule_edit_add_substructure():
    # Setup
    vis_mol = create_benzene_visual()
    request = EditRequest(
        current_molecule=vis_mol,
        user_prompt="Add a methyl group",
        selected_indices=[]
    )

    with patch("app.ai.agent.OpenAI") as mock_openai, \
         patch("app.ai.agent.ChemistryTools") as mock_tools, \
         patch("app.ai.agent.align_and_diff") as mock_align:

        # Mock OpenAI
        tool_call = create_mock_tool_call("add_substructure", '{"anchor_atom_id": 0, "fragment": "methyl"}')

        mock_client = mock_openai.return_value
        mock_client.chat.completions.create.side_effect = [
            MagicMock(choices=[MagicMock(message=MagicMock(tool_calls=[tool_call]))]),
            MagicMock(choices=[MagicMock(message=MagicMock(tool_calls=[]))])
        ]

        # Mock ChemistryTools.add_substructure
        # Return a new molecule object (mocked)
        mock_tools.add_substructure.return_value = MagicMock()
        mock_tools.get_mapped_smiles.return_value = "c1ccccc1"

        # Mock align_and_diff to return a visual molecule with 7 atoms
        mock_vis_mol = VisualMolecule(
            molecule_id="test",
            atoms=[Atom(id=i, element="C", x=0, y=0) for i in range(7)],
            bonds=[]
        )
        mock_align.return_value = mock_vis_mol
        mock_tools.apply_diff_metadata.return_value = mock_vis_mol

        # Execute
        result_mol = await process_molecule_edit(request)

        # Verify
        assert len(result_mol.atoms) == 7
        assert mock_tools.add_substructure.called
