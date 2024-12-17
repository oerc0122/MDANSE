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
from typing import Iterable
from abc import abstractmethod
import traceback

from MDANSE.Framework.AtomMapping import AtomLabel
from .InputFileConfigurator import InputFileConfigurator


class FileWithAtomDataConfigurator(InputFileConfigurator):

    def configure(self, filepath: str) -> None:
        """
        Parameters
        ----------
        filepath : str
            The file path.
        """
        self._original_input = filepath
        super().configure(filepath)
        if self.error_status != "OK":
            return
        try:
            self.parse()
        except Exception as e:
            self.error_status = f"File parsing error {e}: {traceback.format_exc()}"

    @abstractmethod
    def parse(self) -> None:
        """Parse the file."""
        pass

    @abstractmethod
    def atom_labels(self) -> Iterable[AtomLabel]:
        """Yields atom labels"""

    def get_labels(self) -> list[AtomLabel]:
        """
        Returns
        -------
        list[AtomLabel]
            An ordered list of atom labels.
        """
        labels = set()
        for label in self.atom_labels():
            if label not in labels:
                labels.add(label)
        return list(labels)
