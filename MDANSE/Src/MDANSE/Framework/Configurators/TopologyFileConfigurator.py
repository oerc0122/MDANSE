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
import MDAnalysis as mda

from MDANSE.Framework.AtomMapping import AtomLabel
from .FileWithAtomDataConfigurator import FileWithAtomDataConfigurator


class TopologyFileConfigurator(FileWithAtomDataConfigurator):

    def parse(self) -> None:
        # TODO currently MDAnalysis guesses the atom types and masses.
        #  There is a PR https://github.com/MDAnalysis/mdanalysis/pull/3753
        #  which will give us more control over what is guessed. We may
        #  want to change the MDAnalysis guessing options in the future
        #  so that it works better with the MDANSE atom mapping.
        self.atoms = mda.Universe(self["filename"]).atoms

    def get_atom_labels(self) -> list[AtomLabel]:
        """
        Returns
        -------
        list[AtomLabel]
            An ordered list of atom labels.
        """
        labels = []
        for at in self.atoms:
            kwargs = {}
            for arg in ["element", "name", "type", "resname", "mass"]:
                if hasattr(at, arg):
                    kwargs[arg] = getattr(at, arg)
            # the first out of the list above will be the main label
            (k, main_label) = next(iter(kwargs.items()))
            kwargs.pop(k)
            label = AtomLabel(main_label, **kwargs)
            if label not in labels:
                labels.append(label)
        return labels
