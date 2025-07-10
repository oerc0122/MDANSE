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

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.Framework.ConfigDescriptors import (
    AtomMapping,
    BooleanConfigDesc,
    OutputTrajectoryConfigDesc,
    PathConfigDesc,
)
from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Parsers.trj import TrjFile
from MDANSE.Framework.Parsers.xtd import XTDFile
from MDANSE.Framework.Units import measure
from MDANSE.MolecularDynamics.Configuration import (
    PeriodicBoxConfiguration,
    RealConfiguration,
)
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter

FORCE_FACTOR = measure(1.0, "kcal_per_mole/ang", equivalent=True).toval("Da nm/ps2 mol")


class Forcite(Converter):
    """Converts a Forcite trajectory to an MDT trajectory."""

    label = "Forcite"

    xtd_file = PathConfigDesc(
        mode="r",
        label="Input XTD file",
    )
    trj_file = PathConfigDesc(
        mode="r",
        label="Input TRJ file",
    )
    atom_aliases = AtomMapping(
        depends={"trajectory": "trajectory_file"}, label="Atom mapping", default={}
    )
    fold = BooleanConfigDesc(label="Fold coordinates into box")
    output_files = OutputTrajectoryConfigDesc()

    def initialize(self):
        """
        Initialize the job.
        """
        super().initialize()

        self._atomicAliases = self.atom_aliases

        self._xtdfile = XTDFile(self.xtd_file)

        self._chemical_system = ChemicalSystem()
        coordinates = np.vstack([atom.xyz for atom in self._xtdfile._atoms.values()])
        element_list = [atom.element for atom in self._xtdfile._atoms.values()]
        charges = [atom.charge for atom in self._xtdfile._atoms.values()]
        name_list = [atom.name for atom in self._xtdfile._atoms.values()]
        unique_labels = set(name_list)
        label_dict = {label: [] for label in unique_labels}
        for temp_index, atom in enumerate(self._xtdfile._atoms.values()):
            label_dict[atom.name].append(temp_index)

        self._chemical_system.initialise_atoms(element_list, name_list)
        self._chemical_system.add_bonds(self._xtdfile._bonds)
        self._chemical_system.add_labels(label_dict)
        self._chemical_system.find_clusters_from_bonds()

        if self._xtdfile.pbc:
            boxConf = PeriodicBoxConfiguration(
                self._chemical_system, coordinates, self._xtdfile.cell
            )
            real_conf = boxConf.to_real_configuration()
        else:
            coordinates *= measure(1.0, "ang").toval("nm")
            real_conf = RealConfiguration(
                self._chemical_system, coordinates, unit_cell=self._xtdfile._cell
            )

        real_conf.fold_coordinates()
        self._configuration = real_conf

        self._trjfile = TrjFile(self.trj_file)

        # The number of steps of the analysis.
        self.numberOfSteps = self._trjfile.n_frames

        self.frames = self._trjfile.frames

        if self._trjfile.velocities_written:
            self._velocities = np.zeros(
                (self._chemical_system.number_of_atoms, 3), dtype=np.float64
            )
        else:
            self._velocities = None

        if self._trjfile.gradients_written:
            self._gradients = np.zeros(
                (self._chemical_system.number_of_atoms, 3), dtype=np.float64
            )
        else:
            self._gradients = None

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

    def run_step(self, index: int) -> tuple[int, None]:
        """Runs a single step of the job.

        Parameters
        ----------
        index : int
            Index of the loop.

        Returns
        -------
        tuple[int, None]
        """

        frame = next(self.frames)

        # # If the universe is periodic set its shape with the current dimensions of the unit cell.
        conf = self._configuration
        conf["coordinates"][self._trjfile.movable_atoms, :] = frame["pos"]
        if conf.is_periodic:
            if self._trjfile.options.defcel:
                conf.unit_cell = frame["cell"]
            if self.fold:
                conf.fold_coordinates()

        if self._velocities is not None:
            self._velocities[self._trjfile.movable_atoms, :] = frame["vel"]
            conf["velocities"] = self._velocities

        if self._gradients is not None:
            self._gradients[self._trjfile.movable_atoms, :] = frame["force"]
            conf["gradients"] = self._gradients

        self._trajectory.dump_configuration(
            conf,
            frame["step_info"].time,
            units={
                "time": "ps",
                "unit_cell": "nm",
                "coordinates": "nm",
                "velocities": "nm/ps",
                "gradients": "Da nm/ps2",
            },
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

    def finalize(self):
        """
        Finalize the job.
        """
        # Close the output trajectory.
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()

        super().finalize()
