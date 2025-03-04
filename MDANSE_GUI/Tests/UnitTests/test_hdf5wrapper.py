#    This file is part of MDANSE_GUI.
#
#    MDANSE_GUI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
import pytest

import numpy as np

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.MolecularDynamics.Configuration import RealConfiguration
from MDANSE.MolecularDynamics.Trajectory import Trajectory, TrajectoryWriter
from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData

from MDANSE_GUI.MolecularViewer.readers.hdf5wrapper import HDF5Wrapper

N_ATOMS = 4
N_TIMESTEPS = 150


@pytest.fixture(scope="module")
def chemical_system():
    chem_sys = ChemicalSystem("Dummy test system")
    chem_sys.initialise_atoms(N_ATOMS * ["H"])
    return chem_sys


@pytest.fixture(scope="module")
def sample_configuration(chemical_system):
    unit_cell = 15.0 * np.eye(3)
    coords = np.array([[1.0, 1.0, 1.0],
                       [1.0, 2.0, 1.0],
                       [10.0, 1.0, 5.11],
                       [10.0, 2.0, 5.09]], dtype=float)

    temp = RealConfiguration(
        chemical_system,
        coords,
        # unit_cell
    )
    return temp


@pytest.fixture(scope="module")
def sample_trajectory(tmp_path, chemical_system, sample_configuration):
    # here we write to a file
    fname = tmp_path / "output"
    writer = TrajectoryWriter(fname, chemical_system, n_steps=N_TIMESTEPS)
    for ts in range(N_TIMESTEPS):
        writer.dump_configuration(sample_configuration, ts)
    return fname


@pytest.fixture(scope="module")
def wrapper(sample_trajectory):
    traj_object = HDFTrajectoryInputData(sample_trajectory)
    wrapper = HDF5Wrapper(
        sample_trajectory, traj_object.trajectory, traj_object.chemical_system
    )
    return wrapper


def test_frame_read(wrapper: HDF5Wrapper):
    coords = wrapper.read_frame(0)
    assert len(coords.shape) == 2
