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
from __future__ import annotations

from typing import SupportsInt

from qtpy.QtCore import Qt, Slot
from qtpy.QtGui import QIntValidator
from qtpy.QtWidgets import QCheckBox, QLineEdit, QSpinBox

from MDANSE.Framework.Parameters import Integer
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase


class IntegerWidget(WidgetBase):
    def __init__(self, *args, parameter: Integer, **kwargs):
        super().__init__(*args, parameter=parameter, **kwargs)

        self._field = QSpinBox(self._base)

        if self.choices:
            self.reset_choices()
        elif self.ranges:
            self.reset_ranges()

        if self.parameter.optional:
            label = getattr(
                self.parameter, "optional_label", f"Default ({self.parameter.default})"
            )
            self._apply_box = QCheckBox(label, self._base)
            self._apply_box.setTristate(False)
            self._apply_box.checkStateChanged.connect(self.toggle_widgets)
            self._apply_box.setChecked(False)
            self._layout.addWidget(self._apply_box)

        self.set_value(self.default)

        self._layout.addWidget(self._field)

        self.default_labels()
        self.update_labels()
        self.updateValue()

        tooltip_text = self._tooltip or "A single integer number"
        self._field.setToolTip(tooltip_text)

    @Slot()
    def reset_choices(self) -> None:
        choices = self.choices
        self._field.setMinimum(min(choices))
        self._field.setMaximum(max(choices))

    @Slot()
    def reset_ranges(self) -> None:
        mini, maxi = self.ranges
        if mini is not None:
            self._field.setMinimum(mini)
        if maxi is not None:  # Range exclusive.
            self._field.setMaximum(maxi - 1)

    @Slot()
    def toggle_widgets(self):
        if not self.parameter.optional:
            return

        self._field.setEnabled(self._apply_box.checkState() != Qt.CheckState.Checked)

    def default_labels(self):
        """Each Widget should have a default tooltip and label,
        which will be set in this method, unless specific
        values are provided in the settings of the job that
        is being configured."""
        if not self._label_text:
            self._label_text = "IntegerWidget"
        if not self._tooltip:
            self._tooltip = "A single integer number"

    def set_value(self, value: SupportsInt | None) -> None:
        if self.parameter.optional and value is None:
            self._apply_box.setChecked(False)
            return

        self._field.setValue(int(value))
        self.updateValue()

    def get_widget_value(self):
        """Collect the results from the input widgets and return the value."""
        if (
            self.parameter.optional
            and self._apply_box.checkState() != Qt.CheckState.Checked
        ):
            return None

        return self._field.value()
