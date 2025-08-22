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

import numpy as np
from more_itertools import ilen

from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Parameters import (
    Boolean,
    Float,
    OutputTrajectory,
    PathParam,
    to_class,
)
from MDANSE.Framework.Parsers.DCDFile import DCDFile
from MDANSE.IO.MinimalPDBReader import MinimalPDBReader
from MDANSE.MolecularDynamics.Configuration import PeriodicRealConfiguration
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter


class DCD(Converter):
    """Converts a DCD trajectory to an MDT trajectory."""

    label = "DCD"

    pdb_file = PathParam(
        mode="r",
        extensions={"PDB files": "*.pdb"},
        label="Input PDB file",
        callback=to_class(MinimalPDBReader)
    )
    dcd_file = PathParam(
        mode="r",
        extensions={"DCD files": "*.dcd"},
        label="Input DCD file",
        callback=to_class(DCDFile)
        )
    time_step = Float(default=1.0, minimum=1e-9, label="Time step (ps)")
    fold = Boolean(label="Fold coordinates into box")
    output_files = OutputTrajectory()

    def initialize(self):
        """
        Initialize the input parameters and analysis self variables
        """
        super().initialize()

        # The number of steps of the analysis.
        self.numberOfSteps = ilen(self.dcd_file.frames)

        # Create all chemical entities from the PDB file.
        self._chemical_system = self.pdb_file._chemical_system

        # A trajectory is opened for writing.
        self._trajectory = TrajectoryWriter(
            self.output_files.path,
            self._chemical_system,
            self.numberOfSteps,
            positions_dtype=self.output_files.dtype,
            chunking_limit=self.output_files.chunk_size,
            compression=self.output_files.compression,
        )

    def run_step(self, index):
        """Runs a single step of the job.

        Parameters
        ----------
        index : int
            Index of the loop.

        Returns
        -------
        tuple[int, None]
        """
        # The x, y and z values of the current frame.
        unit_cell, config = self.dcd_file.read_step()

        conf = PeriodicRealConfiguration(
            self._trajectory._chemical_system, config, unit_cell
        )

        if self.fold:
            conf.fold_coordinates()

        # The current time.
        time = index * self.time_step

        # Store a snapshot of the current configuration in the output trajectory.
        self._trajectory.dump_configuration(
            conf, time, units={"time": "ps", "unit_cell": "nm", "coordinates": "nm"}
        )

        return index, None

    def combine(self, index, x):
        """Dummy combine step.

        Parameters
        ----------
        _index : int
            Unused.
        _x : None
            Unused.
        """

        pass

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        # Close the output trajectory.
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()

        super().finalize()
