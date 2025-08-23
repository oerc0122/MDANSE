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

import MDAnalysis as mda
from more_itertools import first, first_true

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
from MDANSE.Framework.Units import measure
from MDANSE.MolecularDynamics.Configuration import (
    PeriodicRealConfiguration,
    RealConfiguration,
)
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter
from MDANSE.MolecularDynamics.UnitCell import UnitCell


def get_coords(_desc, coord_files, deps):
    if len(coord_files) <= 1 or deps["coordinate_format"] is None:
        return mda.Universe(
            deps["topology"],
            *coord_files,
            continuous=deps["continuous"],
            format=deps["coordinate_format"],
            topology_format=deps["topology_format"],
        )

    coord_files = [(i, deps["coordinate_format"]) for i in coord_files]
    return mda.Universe(
        deps["topology"],
        coord_files,
        continuous=deps["continuous"],
        topology_format=deps["topology_format"],
    )


class MDAnalysis(Converter):
    """Converts a trajectory to the MDT format using MDAnalysis.

    MDAnalysis reads MD trajectories by specifying
    topology and coordinate files. Multiple files can be used for the
    coordinate files so that trajectories will be stitched together.
    For supported file formats, the continuous option ensures that
    duplicated time-frames will not be added, see
    <a href="https://userguide.mdanalysis.org/stable/reading_and_writing.html">reading and writing</a>.
    For topology and coordinate files supported by MDAnalysis see
    <a href="https://userguide.mdanalysis.org/stable/formats/index.html#formats">formats</a>.
    """

    category = ("Converters", "General")
    label = "MDAnalysis"

    topology_format = SingleChoice(choices=mda._PARSERS.keys())
    topology_file = PathParam(
        mode="r",
        label="Topology file",
    )
    coordinate_format = SingleChoice(choices=mda._PARSERS.keys())
    continuous = Boolean(label="Continuous frame stitching")
    trajectory = PathParam(
        mode="r",
        label="Coordinate file",
        get_depends={
            "topology_format": "topology_format",
            "coordinate_format": "coordinate_format",
            "topology": "topology_file",
            "continuous": "continuous",
        },
        on_get=get_coords,
    )
    time_step = Float(
        label="Time step",
        default=1.0,
        minimum=1e-9,
    )
    atom_aliases = AtomMapping(
        depends={"trajectory": "trajectory"},
        label="Atom mapping",
        default={},
    )
    fold = Boolean(label="Fold coordinates into box")
    output_files = OutputTrajectory()

    def initialize(self):
        """Load the trajectory using MDAnalysis and create the
        trajectory writer.
        """

        self.numberOfSteps = len(self.trajectory.trajectory)

        self._chemical_system = ChemicalSystem()
        element_list = []
        name_list = []
        label_dict = defaultdict(list)

        for at_number, at in enumerate(self.trajectory.atoms):
            kwargs = {
                getattr(at, arg)
                for arg in ("element", "name", "type", "resname", "mass")
                if hasattr(at, arg)
            }

            # the first out of the list above will be the main label
            (k, main_label) = first(kwargs.items())
            del kwargs[k]
            # label_list will be populated too
            tag = kwargs.get(
                first_true(("resname", "type", "name"), pred=kwargs.__contains__)
            )

            if tag:
                label_dict[tag].append(at_number)
            kwargs.pop(k)

            element = get_element_from_mapping(
                self.configuration["atom_aliases"]["value"], main_label, **kwargs
            )

            for arg in ("name", "type", "element"):
                if hasattr(at, arg):
                    name = getattr(at, arg)
                    break
            else:
                name = None

            element_list.append(element)
            name_list.append(name)
        if None in name_list:
            name_list = None
        self._chemical_system.initialise_atoms(element_list, name_list)
        self._chemical_system.add_labels(label_dict)

        self._trajectory = TrajectoryWriter(
            self.output_files.path,
            self._chemical_system,
            self.numberOfSteps,
            positions_dtype=self.output_files.dtype,
            chunking_limit=self.output_files.chunk_size,
            compression=self.output_files.compression,
            initial_charges=getattr(self.trajectory.atoms, "charges", None),
        )
        super().initialize()

    def run_step(self, index: int):
        """For the given frame, read the MDAnalysis trajectory data,
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
        self.trajectory.trajectory[index]

        # convert from MDAnalysis units to MDANSE units
        # see https://userguide.mdanalysis.org/stable/units.html for
        # default units in MDAnalysis
        if self.trajectory.trajectory.ts.triclinic_dimensions is None:
            conf = RealConfiguration(
                self._trajectory._chemical_system,
                self.trajectory.trajectory.ts.positions
                * measure(1.0, "ang").toval("nm"),
            )
        else:
            conf = PeriodicRealConfiguration(
                self._trajectory._chemical_system,
                self.trajectory.trajectory.ts.positions
                * measure(1.0, "ang").toval("nm"),
                UnitCell(
                    self.trajectory.trajectory.ts.triclinic_dimensions
                    * measure(1.0, "ang").toval("nm")
                ),
            )

            if self.fold:
                conf.fold_coordinates()

            if hasattr(self.trajectory.trajectory.ts, "velocities"):
                conf["velocities"] = self.trajectory.trajectory.ts.velocities * measure(
                    1.0, "ang/ps"
                ).toval("nm/ps")

            if hasattr(self.trajectory.trajectory.ts, "forces"):
                conf["gradients"] = self.trajectory.trajectory.ts.forces * measure(
                    1.0, "kJ/mol ang", equivalent=True
                ).toval("Da nm/ps2")

        if self.time_step == 0.0:
            time = index * self.trajectory.trajectory.ts.dt
        else:
            time = index * self.time_step

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

        return index, None

    def combine(self, index, x):
        pass

    def finalize(self):
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()
        super().finalize()
