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
from pathlib import Path

import mdtraj as md
from mdtraj.core.trajectory import _TOPOLOGY_EXTS

from MDANSE.Framework.AtomMapping import AtomLabel
from .FileWithAtomDataConfigurator import FileWithAtomDataConfigurator


class MDTrajTopologyFileConfigurator(FileWithAtomDataConfigurator):

    def parse(self) -> None:
        extension = "".join(Path(self["value"]).suffixes)[1:]

        supported = list(i[1:] for i in _TOPOLOGY_EXTS)
        if extension not in supported:
            raise ValueError(
                f"File '{extension}' not support should be one of the following: {supported}"
            )

        if not self._configurable[self._dependencies["trajectory_files"]]._valid:
            raise RuntimeError(f"Trajectory file not valid")

        trajectory_file = self._configurable[self._dependencies["trajectory_files"]][
            "filename"
        ]
        if self["filename"]:
            self.atoms = md.load(trajectory_file, top=self["filename"]).topology.atoms
        else:
            self.atoms = md.load(trajectory_file).topology.atoms

    def get_atom_labels(self) -> list[AtomLabel]:
        """
        Returns
        -------
        list[AtomLabel]
            An ordered list of atom labels.
        """
        labels = []
        for at in self.atoms:
            label = AtomLabel(
                at.name, symbol=at.element.symbol, residue=at.residue.name
            )
            if label not in labels:
                labels.append(label)
        return labels
