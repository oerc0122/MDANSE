
import os

import pytest

from MDANSE.Framework.AtomSelector.selector import ReusableSelection
from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData


short_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "Converted",
    "short_trajectory_after_changes.mdt",
)
mdmc_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "Converted",
    "Ar_mdmc_h5md.h5",
)
com_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "Converted",
    "com_trajectory.mdt",
)
traj_2vb1 = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "Converted",
    "2vb1.mdt"
)


@pytest.mark.parametrize("trajectory", [short_traj, mdmc_traj, com_traj])
def test_select_all(trajectory):
    traj_object = HDFTrajectoryInputData(trajectory)
    n_atoms = len(traj_object.chemical_system.atom_list)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_all'})
    selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(selection) == n_atoms


@pytest.mark.parametrize("trajectory", [short_traj, mdmc_traj, com_traj])
def test_select_none(trajectory):
    traj_object = HDFTrajectoryInputData(trajectory)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_none'})
    selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(selection) == 0


@pytest.mark.parametrize("trajectory", [short_traj, mdmc_traj, com_traj])
def test_inverted_all_is_none(trajectory):
    traj_object = HDFTrajectoryInputData(trajectory)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_all'})
    reusable_selection.set_selection(None, {'function_name': 'invert_selection'})
    selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(selection) == 0


def test_json_saving_is_reversible():
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_all'})
    reusable_selection.set_selection(None, {'function_name': 'invert_selection'})
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    json_string_2 = another_selection.convert_to_json()
    assert json_string == json_string_2


@pytest.mark.parametrize("trajectory", [short_traj, mdmc_traj, com_traj])
def test_selection_from_json_is_the_same_as_from_runtime(trajectory):
    traj_object = HDFTrajectoryInputData(trajectory)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_all'})
    reusable_selection.set_selection(None, {'function_name': 'invert_selection'})
    selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    selection2 = another_selection.select_in_trajectory(traj_object.trajectory)
    print(f"original: {reusable_selection.operations}")
    print(f"another: {another_selection.operations}")
    assert selection == selection2

def test_select_atoms_selects_by_element():
    traj_object = HDFTrajectoryInputData(short_traj)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'atom_types': ['S']})
    s_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(s_selection) == 208
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'atom_types': ['Cu']})
    cu_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(cu_selection) == 208
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'atom_types': ['Sb']})
    sb_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(sb_selection) == 64
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'atom_types': ['S', 'Sb', 'Cu']})
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    cusbs_selection = another_selection.select_in_trajectory(traj_object.trajectory)
    assert len(cusbs_selection) == 480


def test_select_atoms_selects_by_range():
    traj_object = HDFTrajectoryInputData(short_traj)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'index_range': [15,35]})
    range_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(range_selection) == 20
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'index_range': [470,510]})
    overshoot_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(overshoot_selection) == 10


def test_select_atoms_selects_by_slice():
    traj_object = HDFTrajectoryInputData(short_traj)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'index_slice': [150,350,10]})
    range_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(range_selection) == 20
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_atoms', 'index_slice': [470,510,5]})
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    overshoot_selection = another_selection.select_in_trajectory(traj_object.trajectory)
    assert len(overshoot_selection) == 2


def test_select_molecules_selects_water():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_molecules', 'molecule_names': ['H2 O1']})
    water_selection = reusable_selection.select_in_trajectory(traj_object.trajectory)
    assert len(water_selection) == 28746


def test_select_molecules_selects_protein():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_molecules', 'molecule_names': ['C613 H959 N193 O185 S10']})
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    protein_selection = another_selection.select_in_trajectory(traj_object.trajectory)
    assert len(protein_selection) == 613 + 959 + 193 + 185 + 10


def test_select_molecules_inverted_selects_ions():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_molecules', 'molecule_names': ['C613 H959 N193 O185 S10', 'H2 O1']})
    reusable_selection.set_selection(None, {'function_name': 'invert_selection'})
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    non_molecules_selection = another_selection.select_in_trajectory(traj_object.trajectory)
    assert all([traj_object.chemical_system.atom_list[index] in ['Na', 'Cl'] for index in non_molecules_selection])


def test_select_pattern_selects_water():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    reusable_selection = ReusableSelection()
    reusable_selection.set_selection(None, {'function_name': 'select_pattern', 'rdkit_pattern': "[#8X2;H2](~[H])~[H]"})
    json_string = reusable_selection.convert_to_json()
    another_selection = ReusableSelection()
    another_selection.read_from_json(json_string)
    water_selection = another_selection.select_in_trajectory(traj_object.trajectory)
    assert len(water_selection) == 28746
