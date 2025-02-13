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

from qtpy.QtCore import Signal, Slot
from qtpy.QtWidgets import QGroupBox, QHBoxLayout, QPushButton, QComboBox

from MDANSE.Framework.AtomSelector.general_selection import (
    select_all,
    select_none,
    invert_selection,
)


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
        if self.mode_box.currentIndex == 0:
            return "difference"
        elif self.mode_box.currentIndex == 1:
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

    def parameter_dictionary(self):
        return {"function_name": "select_all"}
