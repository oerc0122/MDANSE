import tempfile
import os
from os import path

import h5py
import numpy as np
import pytest

from MDANSE.Framework.AtomSelector.general_selection import select_all, select_none
from MDANSE.Framework.AtomSelector.atom_selection import select_atoms
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
