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
from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.MolecularDynamics.Trajectory import Trajectory


def select_all(
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
    return set(range(len(trajectory.chemical_system.atom_list)))


def select_none(
    trajectory: Trajectory, **function_parameters: Dict[str, Any]
) -> Set[int]:
    """Returns an empty selection.

    Parameters
    ----------
    selection : Set[int]
        A set of atom indices
    trajectory : Trajectory
        A trajectory instance to which the selection is applied

    Returns
    -------
    Set[int]
        An empty set.
    """
    return set()


def invert_selection(trajectory: Trajectory, selection: Set[int]) -> Set[int]:
    all_indices = select_all(trajectory)
    inverted = all_indices.difference(selection)
    return inverted
