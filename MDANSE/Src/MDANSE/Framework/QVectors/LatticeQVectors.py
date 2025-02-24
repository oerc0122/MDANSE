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

from MDANSE.Framework.QVectors.IQVectors import IQVectors, QVectorsError


class LatticeQVectors(IQVectors):
    is_lattice = True

    def __init__(self, atom_configuration, status=None):
        super(LatticeQVectors, self).__init__(atom_configuration, status)

        if atom_configuration is None:
            raise QVectorsError("No configuration set for the chemical system")

        if not atom_configuration.is_periodic:
            raise QVectorsError(
                "The universe must be periodic for building lattice-based Q vectors"
            )

        self._unit_cell = atom_configuration.unit_cell
