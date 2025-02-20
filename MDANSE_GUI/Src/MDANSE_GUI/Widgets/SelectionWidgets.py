#    This file is part of MDANSE_GUI.
#
#    MDANSE_GUI is free software: you can redistribute it and/or modify
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

from typing import Dict, Any
import json

import numpy as np
from qtpy.QtCore import Signal, Slot
from qtpy.QtWidgets import (
    QGroupBox,
    QHBoxLayout,
    QPushButton,
    QComboBox,
    QLabel,
    QLineEdit,
)

from MDANSE_GUI.InputWidgets.CheckableComboBox import CheckableComboBox
from MDANSE.MolecularDynamics.Trajectory import Trajectory


class BasicSelectionWidget(QGroupBox):
    new_selection = Signal(str)

    def __init__(self, parent=None, widget_label="Atom selection widget"):
        super().__init__(parent)
        layout = QHBoxLayout()
        self.setLayout(layout)
        self.setTitle(widget_label)
        self.add_specific_widgets()
        self.add_standard_widgets()

    def parameter_dictionary(self) -> Dict[str, Any]:
        return {}

    def add_specific_widgets(self):
        return None

    def add_standard_widgets(self):
        self.mode_box = QComboBox(self)
        self.mode_box.setEditable(False)
        self.mode_box.addItems(
            ["Add (union)", "Filter (intersection)", "Remove (difference)"]
        )
        self.commit_button = QPushButton("Apply", self)
        layout = self.layout()
        layout.addWidget(self.mode_box)
        layout.addWidget(self.commit_button)
        self.commit_button.clicked.connect(self.create_selection)

    def get_mode(self) -> str:
        if self.mode_box.currentIndex() == 2:
            return "difference"
        elif self.mode_box.currentIndex() == 1:
            return "intersection"
        else:
            return "union"

    def create_selection(self):
        funtion_parameters = self.parameter_dictionary()
        funtion_parameters["operation_type"] = self.get_mode()
        self.new_selection.emit(json.dumps(funtion_parameters))


class AllAtomSelection(BasicSelectionWidget):
    def __init__(self, parent=None, widget_label="ALL ATOMS"):
        super().__init__(parent, widget_label)

    def add_specific_widgets(self):
        layout = self.layout()
        inversion_button = QPushButton("INVERT selection", self)
        inversion_button.clicked.connect(self.invert_selection)
        layout.addWidget(inversion_button)
        layout.addWidget(QLabel("Add/remove ALL atoms"))

    def invert_selection(self):
        self.new_selection.emit(json.dumps({"function_name": "invert_selection"}))

    def parameter_dictionary(self):
        return {"function_name": "select_all"}


class AtomSelection(BasicSelectionWidget):
    def __init__(
        self, parent=None, trajectory: Trajectory = None, widget_label="Select atoms"
    ):
        self.atom_types = []
        self.atom_names = []
        if trajectory:
            self.atom_types = list(np.unique(trajectory.chemical_system.atom_list))
            if trajectory.chemical_system.name_list:
                self.atom_names = list(np.unique(trajectory.chemical_system.name_list))
        self.selection_types = []
        self.selection_keyword = ""
        if self.atom_types:
            self.selection_types += ["type"]
        if self.atom_names:
            self.selection_types += ["name"]
        super().__init__(parent, widget_label)

    def add_specific_widgets(self):
        layout = self.layout()
        layout.addWidget(QLabel("Select atoms by atom"))
        self.selection_type_combo = QComboBox(self)
        self.selection_type_combo.addItems(self.selection_types)
        self.selection_type_combo.setEditable(False)
        layout.addWidget(self.selection_type_combo)
        self.selection_field = CheckableComboBox(self)
        layout.addWidget(self.selection_field)
        self.selection_type_combo.currentTextChanged.connect(self.switch_mode)
        self.selection_type_combo.setCurrentText(self.selection_types[0])
        self.switch_mode(self.selection_types[0])

    @Slot(str)
    def switch_mode(self, new_mode: str):
        self.selection_field.clear()
        if new_mode == "type":
            self.selection_field.addItems(self.atom_types)
            self.selection_keyword = "atom_types"
        elif new_mode == "name":
            self.selection_field.addItems(self.atom_names)
            self.selection_keyword = "atom_names"

    def parameter_dictionary(self):
        function_parameters = {"function_name": "select_atoms"}
        selection = self.selection_field.checked_values()
        function_parameters[self.selection_keyword] = selection
        return function_parameters


class IndexSelection(BasicSelectionWidget):
    def __init__(self, parent=None, widget_label="Index selection"):
        super().__init__(parent, widget_label)
        self.selection_keyword = "index_list"

    def add_specific_widgets(self):
        layout = self.layout()
        layout.addWidget(QLabel("Select atoms by index"))
        self.selection_type_combo = QComboBox(self)
        self.selection_type_combo.addItems(
            [
                "list",
                "range",
                "slice",
            ]
        )
        self.selection_type_combo.setEditable(False)
        layout.addWidget(self.selection_type_combo)
        self.selection_field = QLineEdit(self)
        layout.addWidget(self.selection_field)
        self.selection_type_combo.currentTextChanged.connect(self.switch_mode)

    @Slot(str)
    def switch_mode(self, new_mode: str):
        self.selection_field.setText("")
        if new_mode == "list":
            self.selection_field.setPlaceholderText("0,1,2")
            self.selection_keyword = "index_list"
            self.selection_separator = ","
        if new_mode == "range":
            self.selection_field.setPlaceholderText("0-20")
            self.selection_keyword = "index_range"
            self.selection_separator = "-"
        if new_mode == "slice":
            self.selection_field.setPlaceholderText("first:last:step")
            self.selection_keyword = "index_slice"
            self.selection_separator = ":"

    def parameter_dictionary(self):
        function_parameters = {"function_name": "select_atoms"}
        selection = self.selection_field.text()
        function_parameters[self.selection_keyword] = [
            int(x) for x in selection.split(self.selection_separator)
        ]
        return function_parameters


class MoleculeSelection(BasicSelectionWidget):
    def __init__(
        self,
        parent=None,
        trajectory: Trajectory = None,
        widget_label="Select molecules",
    ):
        self.molecule_names = []
        if trajectory:
            self.molecule_names = list(trajectory.chemical_system._clusters.keys())
        super().__init__(parent, widget_label)

    def add_specific_widgets(self):
        layout = self.layout()
        layout.addWidget(QLabel("Select molecules named: "))
        self.selection_field = CheckableComboBox(self)
        layout.addWidget(self.selection_field)
        self.selection_field.addItems(self.molecule_names)

    def parameter_dictionary(self):
        function_parameters = {"function_name": "select_molecules"}
        selection = self.selection_field.checked_values()
        function_parameters["molecule_names"] = selection
        return function_parameters


class LabelSelection(BasicSelectionWidget):
    def __init__(
        self,
        parent=None,
        trajectory: Trajectory = None,
        widget_label="Select by label",
    ):
        self.labels = []
        if trajectory:
            self.labels = list(trajectory.chemical_system._labels.keys())
        super().__init__(parent, widget_label)

    def add_specific_widgets(self):
        layout = self.layout()
        layout.addWidget(QLabel("Select atoms with label: "))
        self.selection_field = CheckableComboBox(self)
        layout.addWidget(self.selection_field)
        self.selection_field.addItems(self.labels)

    def parameter_dictionary(self):
        function_parameters = {"function_name": "select_labels"}
        selection = self.selection_field.checked_values()
        function_parameters["atom_labels"] = selection
        return function_parameters


class PatternSelection(BasicSelectionWidget):
    def __init__(
        self,
        parent=None,
        widget_label="SMARTS pattern matching",
    ):
        self.pattern_dictionary = {
            "primary amine": "[#7X3;H2;!$([#7][#6X3][!#6]);!$([#7][#6X2][!#6])](~[H])~[H]",
            "hydroxy": "[#8;H1,H2]~[H]",
            "methyl": "[#6;H3](~[H])(~[H])~[H]",
            "phosphate": "[#15X4](~[#8])(~[#8])(~[#8])~[#8]",
            "sulphate": "[#16X4](~[#8])(~[#8])(~[#8])~[#8]",
            "thiol": "[#16X2;H1]~[H]",
        }
        super().__init__(parent, widget_label)

    def add_specific_widgets(self):
        layout = self.layout()
        layout.addWidget(QLabel("Pick a group"))
        self.selection_field = QComboBox(self)
        layout.addWidget(self.selection_field)
        self.selection_field.addItems(self.pattern_dictionary.keys())
        layout.addWidget(QLabel("pattern:"))
        self.input_field = QLineEdit("", self)
        self.input_field.setPlaceholderText("can be edited")
        layout.addWidget(self.input_field)
        self.selection_field.currentTextChanged.connect(self.update_string)

    @Slot(str)
    def update_string(self, key_string: str):
        if key_string in self.pattern_dictionary:
            self.input_field.setText(self.pattern_dictionary[key_string])

    def parameter_dictionary(self):
        function_parameters = {"function_name": "select_pattern"}
        selection = self.input_field.text()
        function_parameters["rdkit_pattern"] = selection
        return function_parameters
