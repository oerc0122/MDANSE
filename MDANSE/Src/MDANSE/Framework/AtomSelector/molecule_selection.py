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

from typing import Union, Dict, Any, Set
from functools import reduce

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.MolecularDynamics.Trajectory import Trajectory


def select_molecules(
    trajectory: Trajectory, **function_parameters: Dict[str, Any]
) -> Set[int]:
    """Selects all the atoms in the trajectory.

    Parameters
    ----------
    selection : Set[int]
        A set of atom indices
    trajectory : Trajectory
        A trajectory instance to which the selection is applied

    Returns
    -------
    Set[int]
        Set of all the atom indices
    """
    selection = set()
    system = trajectory.chemical_system
    molecule_names = function_parameters.get("molecule_names", None)
    for molecule in molecule_names:
        if molecule in system._clusters:
            selection = selection.union(
                reduce(list.__add__, system._clusters[molecule])
            )
    return selection
