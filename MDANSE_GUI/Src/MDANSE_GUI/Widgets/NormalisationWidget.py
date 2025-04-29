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

from typing import TYPE_CHECKING, List

if TYPE_CHECKING:
    from MDANSE_GUI.Tabs.Models.PlottingContext import PlottingContext

import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import (
    NavigationToolbar2QT as NavigationToolbar2QTAgg,
)
from qtpy.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QComboBox,
    QLabel,
    QCheckBox,
    QSpinBox,
)
from qtpy.QtCore import Slot, Signal, Qt

from MDANSE.MLogging import LOG

from MDANSE_GUI.Widgets.RestrictedSlider import RestrictedSlider
from MDANSE_GUI.Tabs.Plotters.Plotter import (
    Plotter,
    NORMALISATION_DEFAULTS,
    str_to_enum,
    enum_to_str,
)


class NormalisationWidget(QWidget):
    new_values = Signal(dict)

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        layout = QHBoxLayout(self)
        self.setLayout(layout)
        self.combo_sum_average = QComboBox(self)
        self.minspin = QSpinBox(self)
        self.maxspin = QSpinBox(self)
        self.apply_norm = QCheckBox("Normalise data:", self)
        for widget in [
            self.apply_norm,
            QLabel("divide curves by"),
            self.combo_sum_average,
            QLabel("of points from"),
            self.minspin,
            QLabel("to"),
            self.maxspin,
            QLabel("(indices)"),
        ]:
            layout.addWidget(widget)
        self.apply_norm.setChecked(NORMALISATION_DEFAULTS["apply"])
        self.minspin.setValue(NORMALISATION_DEFAULTS["min_index"])
        self.maxspin.setValue(NORMALISATION_DEFAULTS["max_index"])
        self.combo_sum_average.addItems(["average", "sum"])
        self.combo_sum_average.setCurrentText(
            enum_to_str(NORMALISATION_DEFAULTS["operation"])
        )
        self.minspin.valueChanged.connect(self.collect_values)
        self.maxspin.valueChanged.connect(self.collect_values)
        self.apply_norm.checkStateChanged.connect(self.collect_values)
        self.combo_sum_average.currentTextChanged.connect(self.collect_values)

    @Slot(int)
    def update_spinbox_limits(self, curve_length: int):
        newmin = -abs(curve_length)
        newmax = abs(curve_length)
        for sb in [self.minspin, self.maxspin]:
            current_value = sb.value()
            sb.setMinimum(newmin)
            sb.setMaximum(newmax)
            new_value = min(current_value, newmax)
            new_value = max(new_value, newmin)
            sb.setValue(new_value)

    @Slot()
    def collect_values(self):
        results = {
            "apply": self.apply_norm.isChecked(),
            "min_index": self.minspin.value(),
            "max_index": self.maxspin.value(),
            "operation": str_to_enum(self.combo_sum_average.currentText()),
        }
        self.new_values.emit(results)
        return results
