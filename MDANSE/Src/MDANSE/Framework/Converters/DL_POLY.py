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
from typing import Any, Union

import numpy as np
from more_itertools import consume as drop
from more_itertools import take

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.Core.Error import Error
from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Units import measure
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


class HistoryFile(dict):
    _dist_conversion = measure(1.0, "ang").toval("nm")
    _vel_conversion = measure(1.0, "ang/ps").toval("nm/ps")
    _grad_conversion = measure(1.0, "Da ang / ps2").toval("Da nm / ps2")

    def __init__(self, filename):
        super().__init__()

        self["filename"] = filename
        self["instance"] = open(filename, encoding="utf-8")

        drop(self["instance"], 1)
        tagline = self["instance"].readline()
        toks = tagline.split()
        self["keytrj"], self["imcon"], self["natms"] = map(int, toks[:3])

        timeline = self["instance"].readline()
        toks = timeline.split()
        self._time_step = float(toks[5])
        self._first_step = int(toks[1])

        self["instance"].seek(0)

        self["n_frames"] = sum(1 for line in self["instance"] if "timestep" in line)

        self["instance"].seek(0)

        drop(self["instance"], 2)

    def read_step(self, step: int):
        headerline = self["instance"].readline()
        current_step = int(headerline.split()[1])

        time_step = (current_step - self._first_step) * self._time_step
        lines_per_atom = 2 + self["keytrj"]

        if self["imcon"]:
            cell = " ".join(take(3, self["instance"]))
            cell = np.array(cell.split(), dtype=np.float64).reshape(3, 3).T
            cell *= self._dist_conversion
        else:
            cell = None

        charges = np.empty(self["natms"])

        data = {"positions": np.empty((self["natms"], 3))}
        if self["keytrj"] > 0:
            data["velocities"] = np.empty((self["natms"], 3))
        if self["keytrj"] > 1:
            data["gradients"] = np.empty((self["natms"], 3))

        for _ in range(self["natms"]):
            atom = iter(take(lines_per_atom, self["instance"]))
            atom_info = next(atom)
            spec, ind, mass, charge, *_rsd = atom_info.split()
            ind = int(ind) - 1
            mass = float(mass)
            charge = float(charge)
            charges[ind] = float(charge)
            for arr, val in zip(data.values(), atom):
                arr[ind] = np.array(val.split(), dtype=np.float64)

        data["positions"] *= self._dist_conversion
        if "velocities" in data:
            data["velocities"] *= self._vel_conversion
        if "gradients" in data:
            data["gradients"] *= self._grad_conversion

        return time_step, cell, data, charges

    def close(self):
        self["instance"].close()


class DL_POLY(Converter):
    """Converts a DL_POLY trajectory to an MDT trajectory."""

    label = "DL-POLY"

    settings = collections.OrderedDict()
    settings["field_file"] = (
        "FieldFileConfigurator",
        {
            "wildcard": "FIELD files (FIELD*);;All files (*)",
            "default": "INPUT_FILENAME",
            "label": "Input FIELD file",
        },
    )
    settings["history_file"] = (
        "InputFileConfigurator",
        {
            "wildcard": "HISTORY files (HISTORY*);;All files (*)",
            "default": "INPUT_FILENAME",
            "label": "Input HISTORY file",
        },
    )
    settings["atom_aliases"] = (
        "AtomMappingConfigurator",
        {
            "default": "{}",
            "label": "Atom mapping",
            "dependencies": {"input_file": "field_file"},
        },
    )
    settings["fold"] = (
        "BooleanConfigurator",
        {"default": False, "label": "Fold coordinates into box"},
    )
    # settings['output_files'] = ('output_files', {'formats':["HDFFormat"]})
    settings["output_files"] = (
        "OutputTrajectoryConfigurator",
        {
            "formats": ["MDTFormat"],
            "root": "history_file",
            "label": "MDANSE trajectory (filename, format)",
        },
    )

    def initialize(self):
        """
        Initialize the job.
        """
        super().initialize()

        self._atomic_aliases = self.configuration["atom_aliases"]["value"]
        self._field_file = self.configuration["field_file"]
        self._history_file = HistoryFile(self.configuration["history_file"]["filename"])

        # The number of steps of the analysis.
        self.n_steps = int(self._history_file["n_frames"])
        self._chemical_system = ChemicalSystem()

        self._field_file.build_chemical_system(
            self._chemical_system, self._atomic_aliases
        )

        self._trajectory = TrajectoryWriter(
            self.configuration["output_files"]["file"],
            self._chemical_system,
            self.n_steps,
            positions_dtype=self.configuration["output_files"]["dtype"],
            chunking_limit=self.configuration["output_files"]["chunk_size"],
            compression=self.configuration["output_files"]["compression"],
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
        # The x, y and z values of the current frame.
        time, unit_cell, config, charge = self._history_file.read_step(index)

        unit_cell = UnitCell(unit_cell)

        if self._history_file["imcon"]:
            conf = PeriodicRealConfiguration(
                self._trajectory.chemical_system, config["positions"], unit_cell
            )
        else:
            conf = RealConfiguration(
                self._trajectory.chemical_system, config["positions"]
            )

        if self.configuration["fold"]["value"]:
            conf.fold_coordinates()

        if "velocities" in config:
            conf["velocities"] = config["velocities"]
        if "gradients" in config:
            conf["gradients"] = config["gradients"]

        self._trajectory.dump_configuration(
            conf,
            time,
            units={
                "time": "ps",
                "unit_cell": "nm",
                "coordinates": "nm",
                "velocities": "nm/ps",
                "gradients": "Da nm/ps2",
            },
        )

        self._trajectory.write_charges(charge, index)

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

        self._history_file.close()

        # Close the output trajectory.
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()

        super().finalize()
