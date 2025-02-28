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


def select_labels(
    trajectory: Trajectory, **function_parameters: Dict[str, Any]
) -> Set[int]:
    """Selects atoms with a specific label in the trajectory.
    A residue name can be read as a label by MDANSE.

    Parameters
    ----------
    trajectory : Trajectory
        A trajectory instance to which the selection is applied
    function_parameters : Dict[str, Any]
        should include a list of string labels under key "atom_labels"

    Returns
    -------
    Set[int]
        Set of all the atom indices
    """
    system = trajectory.chemical_system
    atom_labels = function_parameters.get("atom_labels", ())
    selection = {system._labels[label] for label in atom_labels if label in system._labels}
    return selection


def select_pattern(
    trajectory: Trajectory, **function_parameters: Dict[str, Any]
) -> Set[int]:
    """Selects atoms according to the SMARTS string given as input.
    This will only work if molecules and bonds have been detected in the system.
    If the bond information was not read from the input trajectory on conversion,
    it can still be determined in a TrajectoryEditor run.

    Parameters
    ----------
    trajectory : Trajectory
        A trajectory instance to which the selection is applied
    function_parameters : Dict[str, Any]
        should include a SMARTS string under key "rdkit_pattern"

    Returns
    -------
    Set[int]
        Set of all the atom indices
    """
    selection = set()
    system = trajectory.chemical_system
    if pattern := function_parameters.get("rdkit_pattern"):
        selection = system.get_substructure_matches(pattern)
    return selection
