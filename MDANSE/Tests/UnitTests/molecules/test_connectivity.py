import os
import pytest
from MDANSE.MolecularDynamics.Connectivity import Connectivity
from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from MDANSE.MolecularDynamics.Trajectory import Trajectory


short_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "Data", "co2gas_md3.mdt"
)


@pytest.fixture
def trajectory() -> HDFTrajectoryInputData:
    trajectory = HDFTrajectoryInputData(short_traj)
    yield trajectory.trajectory


def test_create_connectivity(trajectory: Trajectory):
    conn = Connectivity(trajectory=trajectory)
    print(conn._unique_elements)
    assert len(conn._unique_elements) == 2


def test_find_bonds(trajectory: Trajectory):
    conn = Connectivity(trajectory=trajectory)
    conn.find_bonds()
    assert len(conn._unique_bonds) == 40


def test_find_molecules(trajectory: Trajectory):
    conn = Connectivity(trajectory=trajectory)
    conn.find_molecules()
    assert len(conn._molecules) == 20


def test_rebuild_molecules(trajectory: Trajectory):
    print(trajectory.chemical_system.atom_list)
    conn = Connectivity(trajectory=trajectory)
    conn.find_molecules()
    atoms_before = int(trajectory.chemical_system.number_of_atoms)
    chemical_system = trajectory.chemical_system
    print(conn._molecules)
    chemical_system.rebuild(conn._molecules)
    atoms_after = int(trajectory.chemical_system.number_of_atoms)
    assert atoms_before == atoms_after


def test_identify_molecules(trajectory: Trajectory):
    conn = Connectivity(trajectory=trajectory)
    conn.find_bonds()
    chemical_system = trajectory.chemical_system
    conn.add_bond_information(chemical_system)
    molstrings = []
    for molname, mollist in chemical_system._clusters.items():
        assert len(mollist) == 20
    result = True
    for ms in molstrings[1:]:
        result = result and ms == molstrings[0]
    assert result
    print(chemical_system.unique_molecules())
    assert len(chemical_system.unique_molecules()) == 1
