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

from collections import defaultdict
from math import isclose

import mdtraj as md

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.Framework.AtomMapping import get_element_from_mapping
from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Parameters import (
    AtomMapping,
    Boolean,
    Float,
    OutputTrajectory,
    PathParam,
    SingleChoice,
)
from MDANSE.MolecularDynamics.Configuration import (
    PeriodicRealConfiguration,
    RealConfiguration,
)
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter
from MDANSE.MolecularDynamics.UnitCell import UnitCell


def get_coords(_desc, value, deps):
    if deps["topology"] is not None:
        return md.load(
            value,
            top=deps["topology"],
            discard_overlapping_frames=deps["discard"],
        )

    return md.load(value, discard_overlapping_frames=deps["discard"])


class MDTraj(Converter):
    """Converts a trajectory to the MDT format using MDTraj.

    MDTraj reads MD trajectories by specifying trajectory files and optionally
    a topology file. Multiple files can be used for the trajectory files so that
    trajectories will be stitched together.
    """

    category = ("Converters", "General")
    label = "MDTraj"

    topology_file = PathParam(
        mode="r",
        label="Topology file",
        optional=True,
        default=None,
    )
    discard_overlapping_frames = Boolean(label="Discard overlapping frames")
    coordinate_files = PathParam(
        mode="r",
        label="Coordinate file",
        on_get_depends={
            "topology": "topology_file",
            "discard": "discard_overlapping_frames",
        },
        on_get=get_coords,
    )
    time_step = Float(
        label="Time step",
        default=0.0,
        minimum=0.0,
    )
    atom_aliases = AtomMapping(
        depends={"trajectory": "topology_file"},
        label="Atom mapping",
        default={},
    )
    fold = Boolean(label="Fold coordinates into box")
    output_files = OutputTrajectory()

    def initialize(self):
        """Load the trajectory using MDTraj and create the trajectory writer."""
        super().initialize()

        self.numberOfSteps = self.coordinate_files.n_frames
        mdtraj_to_mdanse = {}

        self._chemical_system = ChemicalSystem()
        elements, atom_names, atom_labels = [], [], defaultdict(list)
        for atnumber, at in enumerate(self.coordinate_files.topology.atoms):
            element = get_element_from_mapping(
                self.atom_aliases,
                at.name,
                symbol=at.element.symbol,
                residue=at.residue.name,
                number=at.element.number,
                mass=at.element.mass,
            )
            elements.append(element)
            atom_names.append(at.name)
            mdtraj_to_mdanse[at.index] = atnumber
            if at.residue.name:
                atom_labels[at.residue.name].append(atnumber)

        self._chemical_system.initialise_atoms(elements, atom_names)
        self._chemical_system.add_labels(atom_labels)
        bonds = [
            (mdtraj_to_mdanse[at1.index], mdtraj_to_mdanse[at2.index])
            for at1, at2 in self.coordinate_files.topology.bonds
        ]
        self._chemical_system.add_bonds(bonds)
        self._chemical_system.find_clusters_from_bonds()

        self._trajectory = TrajectoryWriter(
            self.output_files.path,
            self._chemical_system,
            self.numberOfSteps,
            positions_dtype=self.output_files.dtype,
            chunking_limit=self.output_files.chunk_size,
            compression=self.output_files.compression,
        )

    def run_step(self, index: int):
        """For the given frame, read the MDTraj trajectory data,
        convert and write it out to the MDT file.

        Parameters
        ----------
        index : int
            The frame index.

        Returns
        -------
        tuple[int, None]
            A tuple of the job index and None.
        """
        if self.coordinate_files.unitcell_vectors is None:
            conf = RealConfiguration(
                self._trajectory._chemical_system,
                self.coordinate_files.xyz[index],
            )
        else:
            conf = PeriodicRealConfiguration(
                self._trajectory._chemical_system,
                self.coordinate_files.xyz[index],
                UnitCell(
                    self.coordinate_files.unitcell_vectors[index],
                ),
            )
            if self.fold:
                conf.fold_coordinates()

        # TODO as of 11/12/2024 MDTraj does not read velocity data
        #  there is a discussion about this on GitHub
        #  (https://github.com/mdtraj/mdtraj/issues/1824).
        #  It doesn't look like they have any plans to add this but if
        #  they change their minds then we should update our code to
        #  support this.

        if self.numberOfSteps == 1:
            time = 0
        elif isclose(self.time_step, 0.0):
            time = index * self.coordinate_files.timestep
        else:
            time = index * self.time_step

        self._trajectory.dump_configuration(
            conf,
            time,
            units={
                "time": "ps",
                "unit_cell": "nm",
                "coordinates": "nm",
            },
        )

        return index, None

    def combine(self, index, x):
        pass

    def finalize(self):
        self._trajectory.close()
        super().finalize()
