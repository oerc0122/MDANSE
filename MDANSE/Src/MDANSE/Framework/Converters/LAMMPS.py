#    This file is part of MDANSE.
#
#    MDANSE is free software: you can redistribute it and/or modify
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
from __future__ import annotations

from abc import ABC, abstractmethod
from collections import defaultdict
from contextlib import suppress
from enum import Enum, auto
from itertools import count
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import h5py
import numpy as np
from more_itertools import consume as drop
from more_itertools import ilen, take
from numpy.typing import NDArray

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.Core.Error import Error
from MDANSE.Framework.AtomMapping import get_element_from_mapping
from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Parameters import (
    AtomMapping,
    Boolean,
    Float,
    Integer,
    OutputTrajectory,
    PathParam,
    SingleChoice,
)
from MDANSE.Framework.Parsers import (
    LAMMPSConfigFile,
    LAMMPScustom,
    LAMMPSh5md,
    LAMMPSReader,
    LAMMPSxyz,
)
from MDANSE.Framework.Units import measure
from MDANSE.MLogging import LOG
from MDANSE.MolecularDynamics.Configuration import (
    PeriodicBoxConfiguration,
    RealConfiguration,
)
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter
from MDANSE.MolecularDynamics.UnitCell import UnitCell

if TYPE_CHECKING:
    from MDANSE.Framework.Configurators.ConfigFileConfigurator import (
        ConfigFileConfigurator,
    )

ELECTRON_CHARGE = 1.6021765e-19
DIMS = ("x", "y", "z")

LAMMPS_UNITS = {
    "real": {
        "energy": "kcal_per_mole",
        "time": "fs",
        "length": "ang",
        "velocity": "ang/fs",
        "mass": "Da",
        "charge_conv": 1.0,
    },
    "metal": {
        "energy": "eV",
        "time": "ps",
        "length": "ang",
        "velocity": "ang/ps",
        "mass": "Da",
        "charge_conv": 1.0,
    },
    "si": {
        "energy": "J",
        "time": "s",
        "length": "m",
        "velocity": "m/s",
        "mass": "kg",
        "charge_conv": 1.0 / ELECTRON_CHARGE,
    },
    "cgs": {
        "energy": "erg",  # this will fail
        "time": "s",
        "length": "cm",
        "velocity": "cm/s",
        "mass": "g",
        "charge_conv": 1.0 / 4.8032044e-10,
    },
    "electron": {
        "energy": "Ha",
        "time": "fs",
        "length": "Bohr",
        "velocity": "ang/fs",
        "mass": "Da",
        "charge_conv": 1.0,
    },
    "micro": {
        "energy": "pg*m/s",
        "time": "us",
        "length": "um",
        "velocity": "m/s",
        "mass": "pg",
        "charge_conv": 1.0 / (ELECTRON_CHARGE * 1e12),
    },
    "nano": {
        "energy": "ag*m/s",
        "time": "ns",
        "length": "nm",
        "velocity": "m/s",
        "mass": "ag",
        "charge_conv": 1.0,
    },
}

UnitSchemes = Literal["real", "metal", "si", "cgs", "electron", "micro", "nano"]


class LAMMPSTrajectoryFileError(Error):
    pass


class BoxStyle(Enum):
    """Different styles of "box" provided by LAMMPS.

    Handles conversion to standard 3x3 lattice vectors.
    """

    ORTHOGONAL = auto()
    NONORTHOGONAL = auto()
    TRICLINIC = auto()

    def to_cell(
        self, value: NDArray[float], *, bounds: bool = False
    ) -> tuple[NDArray[float], NDArray[float]]:
        """Convert from LAMMPS box definition to unit cell, origin.

        Parameters
        ----------
        value : NDArray[float]
            LAMMPS unit cell specification.
        bounds : bool
            For LAMMPS' non-orthogonal cells dumps provide *bounds*
            rather than ``[xyz](lo|hi)``. If this is true, convert them.

        Returns
        -------
        NDArray[float]
            Unit cell as 3x3 matrix.
        NDArray[float]
            Origin as 3-vector.

        Examples
        --------
        >>> BoxStyle.ORTHOGONAL.to_cell([2, 5, 2, 5, 2, 5])
        (array([[3., 0., 0.],
               [0., 3., 0.],
               [0., 0., 3.]]), array([2., 2., 2.]))
        >>> BoxStyle.ORTHOGONAL.to_cell([[2, 5], [2, 5], [2, 5]])
        (array([[3., 0., 0.],
               [0., 3., 0.],
               [0., 0., 3.]]), array([2., 2., 2.]))
        >>> BoxStyle.NONORTHOGONAL.to_cell([[2, 5, 1], [2, 5, 2], [2, 5, 4]])
        (array([[3., 0., 0.],
               [1., 3., 0.],
               [2., 4., 3.]]), array([2., 2., 2.]))
        >>> a = BoxStyle.NONORTHOGONAL.to_cell([[-7.6875, 20.0293, 0.],
        ...                                     [0., 11.5962, -7.6875],
        ...                                     [0., 19.6804, 0.]], bounds=True)
        >>> b = BoxStyle.NONORTHOGONAL.to_cell([[0., 20.0293, 0.],
        ...                                     [0., 11.5962, -7.6875],
        ...                                     [0., 19.6804, 0.]], bounds=False)
        >>> np.allclose(a[0], b[0]), np.allclose(a[1], b[1])
        (True, True)
        >>> BoxStyle.TRICLINIC.to_cell([[1, 2, 3, 2], [4, 5, 6, 2], [7, 8, 9, 2]])
        (array([[1., 2., 3.],
               [4., 5., 6.],
               [7., 8., 9.]]), array([2., 2., 2.]))
        """
        value = np.array(value, dtype=float)
        if self is BoxStyle.ORTHOGONAL:
            value = value.reshape((3, 2))
            unit_cell, origin = np.diag(value[:, 1] - value[:, 0]), value[:, 0]
        elif self is BoxStyle.NONORTHOGONAL:
            value = value.reshape((3, 3))
            xy, xz, yz = value[:, 2]

            if bounds:
                value[0, 0] -= min(0.0, xy, xz, xy + xz)
                value[0, 1] -= max(0.0, xy, xz, xy + xz)
                value[1, 0] -= min(0.0, yz)
                value[1, 1] -= max(0.0, yz)

            unit_cell, origin = np.diag(value[:, 1] - value[:, 0]), value[:, 0]
            unit_cell[1, 0] = xy
            unit_cell[2, 0] = xz
            unit_cell[2, 1] = yz
        elif self is BoxStyle.TRICLINIC:
            value = value.reshape((3, 4))
            unit_cell, origin = value[:3, :3], value[:, 3]

        return unit_cell, origin


class LAMMPS(Converter):
    """Converts a LAMMPS trajectory to an MDT trajectory."""

    label = "LAMMPS"

    _READERS = {
        "custom": LAMMPScustom,
        "xyz": LAMMPSxyz,
        "h5md": LAMMPSh5md,
    }

    config_file = PathParam(
        mode="r",
        label="LAMMPS configuration file.",
        extensions={"LAMMPS config": "*.config"},
        default="INPUT_FILENAME.config",
    )
    trajectory_file = PathParam(
        mode="r",
        label="LAMMPS trajectory file.",
        extensions={
            "XYZ files": "*.xyz",
            "Dump files": "*.dump",
            "HDF5 files": "*.h5",
            "LAMMPS files": "*.lammps",
        },
        default="INPUT_FILENAME.lammps",
    )
    trajectory_format = SingleChoice(
        label="LAMMPS trajectory format",
        choices=("custom", "xyz", "h5md"),
        default="custom",
    )
    lammps_units = SingleChoice(
        label="LAMMPS unit system",
        default="real",
        choices=(
            "From config",
            "real",
            "metal",
            "si",
            "cgs",
            "electron",
            "micro",
            "nano",
        ),
    )
    atom_type = SingleChoice(
        label="LAMMPS atom type",
        choices=[
            "From config",
            "angle",
            "atomic",
            "body",
            "bond",
            "bpm/sphere",
            "charge",
            "dielectric",
            "dipole",
            "dpd",
            "edpd",
            "electron",
            "ellipsoid",
            "full",
            "hybrid",
            "line",
            "mdpd",
            "molecular",
            "peri",
            "rheo",
            "rheo/thermal",
            "smd",
            "sph",
            "sphere",
            "spin",
            "tdpd",
            "template",
            "tri",
            "wavepacket",
        ],
        default="From config",
    )
    time_step = Float(
        label="Time step (lammps units, depends on unit system)",
        default=1.0,
        minimum=1e-9,
    )
    n_steps = Integer(
        label="Number of time steps (0 for automatic detection)",
        default=0,
        minimum=0,
    )
    atom_aliases = AtomMapping(
        depends={"trajectory": "trajectory_file"},
        label="Atom mapping",
        default={},
    )
    fold = Boolean(label="Fold coordinates into box")
    output_files = OutputTrajectory()

    def initialize(self):
        """
        Initialize the job.
        """
        super().initialize()

        # The number of steps of the analysis.
        self.numberOfSteps = self.n_steps

        self._lammps_config = LAMMPSConfigFile(self.config_file)
        self._lammps_config.parse()

        if self.lammps_units == "From config":
            if "units" not in self.config_file:
                raise KeyError("Units not found in lammps config")
            self.lammps_units = self._lammps_config["units"]

        self._reader = self.create_reader(self.trajectory_format)

        self._reader.set_units(self.lammps_units)
        self._chemical_system = self.parse_first_step(self.atom_aliases)
        self._chemical_system.find_clusters_from_bonds()

        if self.numberOfSteps == 0:
            self.numberOfSteps = self._reader.get_time_steps(self.trajectory_file)

        charges_single_cell = self._lammps_config["charges"]
        charges_single_cell *= self._reader.units["charge_conv"]

        # Assume replicate
        if len(charges_single_cell) < self._chemical_system.number_of_atoms:
            charges = list(charges_single_cell) * int(
                self._chemical_system.number_of_atoms // len(charges_single_cell)
            )
        else:
            charges = list(charges_single_cell)

        # A trajectory is opened for writing.
        self._trajectory = TrajectoryWriter(
            self.output_files.path,
            self._chemical_system,
            self.numberOfSteps,
            positions_dtype=self.output_files.dtype,
            chunking_limit=self.output_files.chunk_size,
            compression=self.output_files.compression,
            initial_charges=charges,
        )

        self._start = 0
        self._reader.open_file(self.trajectory_file)
        self._reader.set_output(self._trajectory)

    def create_reader(
        self, trajectory_type: Literal["custom", "xyz", "h5md"]
    ) -> LAMMPSReader:
        """Create the required reader.

        Parameters
        ----------
        trajectory_type : Literal["custom", "xyz", "h5md"]
            Type of file to load.

        Returns
        -------
        LAMMPSReader
            Configured reader.
        """
        return self._READERS[trajectory_type](
            lammps_units=self.lammps_units,
            atom_type=self.atom_type,
            timestep=self.time_step,
            fold_coordinates=self.fold,
        )

    def run_step(self, index) -> tuple[int, None]:
        """Runs a single step of the job.

        Parameters
        ----------
        index : int.
            The index of the step.

        Returns
        -------
        int
            Index of job step.
        """

        val1, val2 = self._reader.run_step(index)

        self._start = val2

        return val1, None

    def combine(self, index: int, x: Any) -> None:
        """No work to be done.

        Parameters
        ----------
        index : int.
            The index of the step.
        x : Any
            Data.
        """

    def finalize(self) -> None:
        """
        Finalize the job.
        """

        self._reader.close()

        # Close the output trajectory.
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()

        super().finalize()

    def parse_first_step(self, aliases: dict[str, dict[str, str]]) -> ChemicalSystem:
        """Wrapper for `parse_first_step` on readers.

        Parameters
        ----------
        aliases : dict[str, dict[str, str]]
            Mapping of atomic aliases to elements.

        Returns
        -------
        ChemicalSystem
            Initialised structure.

        See Also
        --------
        LAMMPScustom.parse_first_step
        LAMMPSxyz.parse_first_step
        LAMMPSh5md.parse_first_step
        """

        self._reader.open_file(self.trajectory_file)
        chemical_system = self._reader.parse_first_step(aliases, self._lammps_config)
        self._reader.close()
        return chemical_system
