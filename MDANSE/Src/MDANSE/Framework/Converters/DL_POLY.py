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

import collections
from typing import Any

import numpy as np
from more_itertools import consume as drop

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.Core.Error import Error
from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Parameters import (
    AtomMapping,
    Boolean,
    OutputTrajectory,
    PathParam,
)
from MDANSE.Framework.Parsers import DLPField, DLPHistory
from MDANSE.MolecularDynamics.Configuration import (
    PeriodicRealConfiguration,
    RealConfiguration,
)
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter
from MDANSE.MolecularDynamics.UnitCell import UnitCell


class HistoryFileError(Error):
    pass


class DL_POLYConverterError(Error):
    pass


class DL_POLY(Converter):
    """Converts a DL_POLY trajectory to an MDT trajectory."""

    label = "DL-POLY"

    field_file = PathParam(
        mode="r",
        extensions={"FIELD files": "FIELD*"},
        default="INPUT_FILENAME",
        label="Input FIELD file",
    )
    history_file = PathParam(
        mode="r",
        extensions={"HISTORY files": "HISTORY*"},
        default="INPUT_FILENAME",
        label="Input HISTORY file",
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

        self._atomic_aliases = self.atom_aliases
        self._field_file = DLPField(self.field_file)
        self._history_file = DLPHistory(self.history_file)

        self.frames = self._history_file.frames

        # The number of steps of the analysis.
        self.numberOfSteps = self._history_file.n_frames
        self._chemical_system = ChemicalSystem()

        self._field_file.build_chemical_system(
            self._chemical_system, self._atomic_aliases
        )

        self._trajectory = TrajectoryWriter(
            self.output_files.path,
            self._chemical_system,
            self.numberOfSteps,
            positions_dtype=self.output_files.dtype,
            chunking_limit=self.output_files.chunk_size,
            compression=self.output_files.compression,
            initial_charges=self._field_file.get_atom_charges(),
        )

    def run_step(self, index: int) -> tuple[int, None]:
        """Runs a single step of the job.

        Parameters
        ----------
        index : int
            The index of the loop.

        Notes
        -----
        The argument index is note the index of the frame.
        the index of the step.

        Returns
        -------
        int
            Index.
        """

        frame = next(self.frames)

        if self._history_file.imcon:
            conf = PeriodicRealConfiguration(
                self._trajectory.chemical_system, frame["positions"], frame["cell"]
            )
        else:
            conf = RealConfiguration(
                self._trajectory.chemical_system, frame["positions"]
            )

        if self.fold:
            conf.fold_coordinates()

        if "velocities" in frame:
            conf["velocities"] = frame["velocities"]
        if "gradients" in frame:
            conf["gradients"] = frame["gradients"]

        self._trajectory.dump_configuration(
            conf,
            frame["step"],
            units={
                "time": "ps",
                "unit_cell": "nm",
                "coordinates": "nm",
                "velocities": "nm/ps",
                "gradients": "Da nm/ps2",
            },
        )

        self._trajectory.write_charges(frame["charge"], index)

        return index, None

    def combine(self, index: int, x: Any):
        """Join job steps.

        Parameters
        ----------
        index : int
            Current index.
        x : Any
            Misc data.
        """
        pass

    def finalize(self):
        """
        Finalize the job.
        """

        # Close the output trajectory.
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()

        super().finalize()
