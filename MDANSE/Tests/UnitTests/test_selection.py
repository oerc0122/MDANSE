
import os

import pytest

from MDANSE.Framework.AtomSelector.general_selection import select_all, select_none, invert_selection
from MDANSE.Framework.AtomSelector.atom_selection import select_atoms
from MDANSE.Framework.AtomSelector.molecule_selection import select_molecules
from MDANSE.Framework.AtomSelector.group_selection import select_labels, select_pattern
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
    selection = select_all(traj_object.trajectory)
    assert len(selection) == n_atoms


@pytest.mark.parametrize("trajectory", [short_traj, mdmc_traj, com_traj])
def test_select_none(trajectory):
    traj_object = HDFTrajectoryInputData(trajectory)
    selection = select_none(traj_object.trajectory)
    assert len(selection) == 0


@pytest.mark.parametrize("trajectory", [short_traj, mdmc_traj, com_traj])
def test_inverted_none_is_all(trajectory):
    traj_object = HDFTrajectoryInputData(trajectory)
    none_selection = select_none(traj_object.trajectory)
    all_selection = select_all(traj_object.trajectory)
    inverted_none = invert_selection(traj_object.trajectory, none_selection)
    assert all_selection == inverted_none


def test_select_atoms_selects_by_element():
    traj_object = HDFTrajectoryInputData(short_traj)
    s_selection = select_atoms(traj_object.trajectory, atom_types=['S'])
    assert len(s_selection) == 208
    cu_selection = select_atoms(traj_object.trajectory, atom_types=['Cu'])
    assert len(cu_selection) == 208
    sb_selection = select_atoms(traj_object.trajectory, atom_types=['Sb'])
    assert len(sb_selection) == 64
    cusbs_selection = select_atoms(traj_object.trajectory, atom_types=['Cu','Sb', 'S'])
    assert len(cusbs_selection) == 480


def test_select_atoms_selects_by_range():
    traj_object = HDFTrajectoryInputData(short_traj)
    range_selection = select_atoms(traj_object.trajectory, index_range=[15,35])
    assert len(range_selection) == 20
    overshoot_selection = select_atoms(traj_object.trajectory, index_range=[470,510])
    assert len(overshoot_selection) == 10


def test_select_atoms_selects_by_slice():
    traj_object = HDFTrajectoryInputData(short_traj)
    range_selection = select_atoms(traj_object.trajectory, index_slice=[150,350, 10])
    assert len(range_selection) == 20
    overshoot_selection = select_atoms(traj_object.trajectory, index_slice=[470,510,5])
    assert len(overshoot_selection) == 2


def test_select_molecules_selects_water():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    water_selection = select_molecules(traj_object.trajectory, molecule_names = ['H2 O1'])
    assert len(water_selection) == 28746


def test_select_molecules_selects_protein():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    protein_selecton = select_molecules(traj_object.trajectory, molecule_names = ['C613 H959 N193 O185 S10'])
    assert len(protein_selecton) == 613 + 959 + 193 + 185 + 10


def test_select_molecules_inverted_selects_ions():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    all_molecules_selection = select_molecules(traj_object.trajectory, molecule_names = ['C613 H959 N193 O185 S10', 'H2 O1'])
    non_molecules_selection = invert_selection(traj_object.trajectory, all_molecules_selection)
    assert all([traj_object.chemical_system.atom_list[index] in ['Na', 'Cl'] for index in non_molecules_selection])


def test_select_pattern_selects_water():
    traj_object = HDFTrajectoryInputData(traj_2vb1)
    water_selection = select_pattern(traj_object.trajectory, rdkit_pattern="[#8X2;H2](~[H])~[H]")
    assert len(water_selection) == 28746
