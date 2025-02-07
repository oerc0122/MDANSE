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

from typing import TYPE_CHECKING
from contextlib import suppress

from qtpy.QtCore import Signal, Slot
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QPushButton,
    QLineEdit,
    QDoubleSpinBox,
    QSpinBox,
    QGroupBox,
)

if TYPE_CHECKING:
    from MDANSE_GUI.MolecularViewer.MolecularViewer import MolecularViewer


class TraceWidget(QWidget):

    new_atom_trace = Signal(dict)
    remove_atom_trace = Signal(int)

    def __init__(self, parent):
        super().__init__(parent)
        self._molviewer = None

        self.setWindowTitle("Add atom trace to the view")
        layout = QVBoxLayout(self)
        self.setLayout(layout)

        self._n_atoms = 0
        self._opacity = 0.5
        self._color = (0, 0.5, 0.75)
        self._iso_percentile = 90
        self.populate_layout()

    def initialise_values(self, viewer: "MolecularViewer"):
        self._molviewer = viewer
        self._n_atoms = viewer._n_atoms
        self._opacity = 0.5
        self._color = (0, 0.5, 0.75)
        self._iso_percentile = 90

    def populate_layout(self):
        layout = self.layout()
        self.add_trace_button = QPushButton("Calculate and add atom trace", self)
        self.remove_trace_button = QPushButton("Remove atom trace", self)
        self.add_trace_button.clicked.connect(self.new_trace_details)
        self.remove_trace_button.clicked.connect(self.remove_trace)
        self._atom_spinbox = QSpinBox(self)
        self._surface_spinbox = QSpinBox(self)
        self._fraction_spinbox = QSpinBox(self)
        self._grid_spinbox = QSpinBox(self)
        self._opacity_spinbox = QDoubleSpinBox(self)
        self._colour_lineedit = QLineEdit("0,128,192", self)
        for sbox in [
            self._atom_spinbox,
            self._surface_spinbox,
            self._fraction_spinbox,
            self._opacity_spinbox,
        ]:
            sbox.setMinimum(0)
            sbox.setValue(0)
        self.update_limits()
        self._fraction_spinbox.setMaximum(100)
        self._fraction_spinbox.setValue(5)
        self._grid_spinbox.setMaximum(10)
        self._grid_spinbox.setMinimum(1)
        self._grid_spinbox.setValue(3)
        self._opacity_spinbox.setMaximum(1.0)
        self._opacity_spinbox.setValue(0.5)
        self._opacity_spinbox.setSingleStep(0.01)
        for label, widget in [
            ("Selected atom index: ", self._atom_spinbox),
            ("Sampling step (1=coarse, 10=fine)", self._grid_spinbox),
            ("Trace percentile for isovalue", self._fraction_spinbox),
            ("Isosurface opacity", self._opacity_spinbox),
            ("Isosurface colour (R,G,B)", self._colour_lineedit),
        ]:
            temp_box = QGroupBox(label, self)
            temp_layout = QVBoxLayout(temp_box)
            temp_layout.addWidget(widget)
            layout.addWidget(temp_box)
        layout.addWidget(self.add_trace_button)
        for label, widget in [
            ("Remove the surface with index: ", self._surface_spinbox)
        ]:
            temp_box = QGroupBox(label, self)
            temp_layout = QVBoxLayout(temp_box)
            temp_layout.addWidget(widget)
            layout.addWidget(temp_box)
        layout.addWidget(self.remove_trace_button)

    @Slot()
    def update_limits(self):
        if self._molviewer is None:
            return
        self._atom_spinbox.setMaximum(max(self._molviewer._n_atoms - 1, 0))
        self._surface_spinbox.setMaximum(max(len(self._molviewer._surfaces) - 1, 0))
        self.enable_buttons()

    def enable_buttons(self):
        if self._molviewer is None:
            return
        self.remove_trace_button.setEnabled(len(self._molviewer._surfaces) != 0)
        self.add_trace_button.setEnabled(self._molviewer._n_atoms > 0)

    def get_values(self):
        params = {
            "atom_number": 0,
            "fine_sampling": 3,
            "surface_colour": (0, 0.5, 0.75),
            "surface_opacity": 0.5,
            "trace_cutoff": 5,
            "surface_number": -1,
        }
        params["atom_number"] = self._atom_spinbox.value()
        params["surface_number"] = self._surface_spinbox.value()
        with suppress(ValueError, TypeError):
            params["surface_colour"] = [
                float(x) / 256 for x in self._colour_lineedit.text().split(",")
            ]
        params["trace_cutoff"] = self._fraction_spinbox.value()
        params["fine_sampling"] = self._grid_spinbox.value()
        params["surface_opacity"] = self._opacity_spinbox.value()
        return params

    def new_trace_details(self):
        params = self.get_values()
        self.new_atom_trace.emit(params)

    def remove_trace(self):
        params = self.get_values()
        self.remove_atom_trace.emit(params["surface_number"])
