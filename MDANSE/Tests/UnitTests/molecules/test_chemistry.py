import os
import pytest
import numpy as np
from rdkit.Chem.rdmolops import SanitizeMol
from rdkit.Chem.rdmolops import GetMolFrags
from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from MDANSE.Chemistry.Structrures import Topology


short_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "..", "Data", "co2gas_md3.mdt"
)


@pytest.fixture(scope="module")
def trajectory():
    trajectory = HDFTrajectoryInputData(short_traj)
    yield trajectory


def test_unit_cell(trajectory: HDFTrajectoryInputData):
    chem_system = trajectory.chemical_system
    configuration = chem_system.configuration
    unit_cell = configuration.unit_cell
    print(unit_cell.abc_and_angles)


def test_molecule_finder(trajectory: HDFTrajectoryInputData):
    chem_system = trajectory.chemical_system
    configuration = chem_system.configuration
    coordinates = configuration._variables["coordinates"]
    print(coordinates.shape)
