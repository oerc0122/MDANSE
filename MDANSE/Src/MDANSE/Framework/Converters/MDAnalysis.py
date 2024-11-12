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
import collections

import MDAnalysis as mda

from MDANSE.Framework.Units import measure
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter
from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Chemistry.ChemicalEntity import ChemicalSystem, Atom
from MDANSE.Framework.AtomMapping import get_element_from_mapping
from MDANSE.MolecularDynamics.Configuration import (
    PeriodicRealConfiguration,
    RealConfiguration,
)
from MDANSE.MolecularDynamics.UnitCell import UnitCell


class MDAnalysis(Converter):

    label = "MDAnalysis"
    settings = collections.OrderedDict()
    settings["topology_file"] = (
        "TopologyFileConfigurator",
        {
            "wildcard": "All files (*)",
            "default": "INPUT_FILENAME",
            "label": "The topology file",
        },
    )
    settings["coordinate_file"] = (
        "InputFileConfigurator",
        {
            "wildcard": "All files (*)",
            "default": "INPUT_FILENAME",
            "label": "The coordinate files",
        },
    )
    settings["atom_aliases"] = (
        "AtomMappingConfigurator",
        {
            "default": "{}",
            "label": "Atom mapping",
            "dependencies": {"input_file": "topology_file"},
        },
    )
    settings["time_step"] = (
        "FloatConfigurator",
        {
            "label": "Time step (ps) - use 0.0 for MDAnalysis default",
            "default": 0.0,
            "mini": 0.0,
        },
    )
    settings["fold"] = (
        "BooleanConfigurator",
        {"default": False, "label": "Fold coordinates into box"},
    )
    settings["output_files"] = (
        "OutputTrajectoryConfigurator",
        {
            "label": "MDANSE trajectory (filename, format)",
            "formats": ["MDTFormat"],
            "root": "config_file",
        },
    )

    def initialize(self):
        self.u = mda.Universe(
            self.configuration["topology_file"]["filename"],
            self.configuration["coordinate_file"]["filename"],
        )

        self.numberOfSteps = len(self.u.trajectory)

        self._chemical_system = ChemicalSystem()

        for at in self.u.atoms:
            kwargs = {}
            for arg in ["name", "resname", "mass"]:
                if hasattr(at, arg):
                    kwargs[arg] = getattr(at, arg)
            element = get_element_from_mapping(
                self.configuration["atom_aliases"]["value"], at.type, **kwargs
            )
            at = Atom(symbol=element, name=at.type)
            self._chemical_system.add_chemical_entity(at)

        self._trajectory = TrajectoryWriter(
            self.configuration["output_files"]["file"],
            self._chemical_system,
            self.numberOfSteps,
            positions_dtype=self.configuration["output_files"]["dtype"],
            compression=self.configuration["output_files"]["compression"],
        )
        super().initialize()

    def run_step(self, index):
        self.u.trajectory[index]

        # convert from MDAnalysis units to MDANSE units
        # see https://userguide.mdanalysis.org/stable/units.html for
        # default units in MDAnalysis
        if self.u.trajectory.ts.triclinic_dimensions is None:
            conf = RealConfiguration(
                self._trajectory._chemical_system,
                self.u.trajectory.ts.positions * measure(1.0, "ang").toval("nm"),
            )
        else:
            conf = PeriodicRealConfiguration(
                self._trajectory._chemical_system,
                self.u.trajectory.ts.positions * measure(1.0, "ang").toval("nm"),
                UnitCell(
                    self.u.trajectory.ts.triclinic_dimensions
                    * measure(1.0, "ang").toval("nm")
                ),
            )
            if self.configuration["fold"]["value"]:
                conf.fold_coordinates()

        self._trajectory._chemical_system.configuration = conf

        if float(self.configuration["time_step"]["value"]) == 0.0:
            time = index * self.u.trajectory.ts.dt
        else:
            time = index * float(self.configuration["time_step"]["value"])

        self._trajectory.dump_configuration(
            time, units={"time": "ps", "unit_cell": "nm", "coordinates": "nm"}
        )

        return index, None

    def combine(self, index, x):
        pass

    def finalize(self):
        self._trajectory.close()
        super(MDAnalysis, self).finalize()
