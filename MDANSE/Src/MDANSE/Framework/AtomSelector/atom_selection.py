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


def select_atoms(trajectory: Trajectory, **function_parameters: Dict[str, Any]) -> Set[int]:
    """Selects specific atoms in the trajectory. These can be selected based
    on indices, atom type or trajectory-specific atom name.
    The atom type is normally the chemical element, while
    the atom name can be more specific and depend on the
    force field used.

    Parameters
    ----------
    trajectory : Trajectory
        A trajectory instance to which the selection is applied
    function_parameters : Dict[str, Any]
        may include any combination of the following
        "index_list" : List[int]
        "index_range" : a start, stop pair of indices for range(start, stop)
        "index_slice" : s start, stop, step sequence of indices for range(start, stop, step)
        "atom_types" : List[str]
        "atom_names" : List[str]

    Returns
    -------
    Set[int]
        Set of all the atom indices
    """
    selection = set()
    system = trajectory.chemical_system
    element_list = system.atom_list
    name_list = system.name_list
    indices = set(range(len(element_list)))
    index_list = function_parameters.get("index_list")
    index_range = function_parameters.get("index_range")
    index_slice = function_parameters.get("index_slice")
    if index_list is not None:
        selection |= indices & index_list
    if index_range is not None:
        selection |= indices & set(range(*index_range))
    if index_slice is not None:
        selection |= indices & set(range(*index_slice))
    atom_types = function_parameters.get("atom_types", ())
    if atom_types:
        new_indices = {index for index in indices if element_list[index] in atom_types}
        selection |= new_indices
    atom_names = function_parameters.get("atom_names", ())
    if atom_names:
        new_indices = {index for index in indices if name_list[index] in atom_names}
        selection |= new_indices
    return selection
