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

import json
from typing import Union, Dict, Any, Set
from MDANSE.MolecularDynamics.Trajectory import Trajectory
from MDANSE.Framework.AtomSelector.general_selection import (
    select_all,
    select_none,
    invert_selection,
)
from MDANSE.Framework.AtomSelector.atom_selection import select_atoms
from MDANSE.Framework.AtomSelector.molecule_selection import select_molecules
from MDANSE.Framework.AtomSelector.group_selection import select_labels, select_pattern


function_lookup = {
    function.__name__: function
    for function in [
        select_all,
        select_none,
        invert_selection,
        select_atoms,
        select_molecules,
        select_labels,
        select_pattern,
    ]
}


class ReusableSelection:
    """A reusable sequence of operations which, when applied
    to a trajectory, returns a set of atom indices based
    on the specified criteria.
    """

    def __init__(self) -> None:
        """
        Parameters
        ----------
        trajectory: Trajectory
            The chemical system to apply the selection to.
        """
        self.reset()

    def reset(self):
        self.system = None
        self.trajectory = None
        self.all_idxs = set()
        self.operations = {}

    def set_selection(
        self, number: Union[int, None] = None, function_parameters: Dict[str, Any] = {}
    ):
        if number is None:
            number = len(self.operations)
        else:
            try:
                number = int(number)
            except TypeError:
                number = len(self.operations)
        self.operations[number] = function_parameters

    def validate_selection_string(
        self, json_string: str, trajectory: Trajectory, current_selection: Set[int]
    ) -> bool:
        """Checks if the selection operation encoded in the input JSON string
        will add any new atoms to the current selection on the given trajectory.

        Parameters
        ----------
        json_string : str
            new selection operation in a JSON string
        trajectory : Trajectory
            a trajectory instance for which current_selection is defined
        current_selection : Set[int]
            set of currently selected atom indices

        Returns
        -------
        bool
            True if the selection adds atoms, False otherwise
        """
        function_parameters = json.loads(json_string)
        if len(self.operations) == 0:
            return True
        function_name = function_parameters.get("function_name", "select_all")
        if function_name == "invert_selection":
            selection = invert_selection(trajectory, current_selection)
            return True
        else:
            operation_type = function_parameters.get("operation_type", "union")
            function = function_lookup[function_name]
            temp_selection = function(trajectory, **function_parameters)
            if operation_type == "union":
                selection = selection.union(temp_selection)
            elif operation_type == "intersection":
                selection = selection.intersection(temp_selection)
            elif operation_type == "difference":
                selection = selection.difference(temp_selection)
            else:
                selection = temp_selection
        if len(selection.difference(current_selection)) > 0 and operation_type == "union":
            return True
        elif len(current_selection.difference(selection)) > 0 and operation_type != "union":
            return True
        return False

    def select_in_trajectory(self, trajectory: Trajectory) -> Set[int]:
        selection = set()
        self.all_idxs = set(range(len(trajectory.chemical_system.atom_list)))
        sequence = sorted([int(x) for x in self.operations.keys()])
        if len(sequence) == 0:
            return self.all_idxs
        for number in sequence:
            function_parameters = self.operations[number]
            function_name = function_parameters.get("function_name", "select_all")
            if function_name == "invert_selection":
                selection = self.all_idxs.difference(selection)
            else:
                operation_type = function_parameters.get("operation_type", "union")
                function = function_lookup[function_name]
                temp_selection = function(trajectory, **function_parameters)
                if operation_type == "union":
                    selection = selection.union(temp_selection)
                elif operation_type == "intersection":
                    selection = selection.intersection(temp_selection)
                elif operation_type == "difference":
                    selection = selection.difference(temp_selection)
                else:
                    selection = temp_selection
        return selection

    def convert_to_json(self) -> str:
        """For the purpose of storing the selection independent of the
        trajectory it is acting on, this method encodes the sequence
        of selection operations as a string.

        Returns
        -------
        str
            All the operations of this selection, encoded as string
        """
        return json.dumps(self.operations)

    def load_from_json(self, json_string: str):
        """_summary_

        Parameters
        ----------
        json_string : str
            A sequence of selection operations, encoded as a JSON string
        """
        json_setting = json.loads(json_string)
        for k0, v0 in json_setting.items():
            if isinstance(v0, dict):
                self.set_selection(k0, v0)
