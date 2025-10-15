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

from qtpy.QtCore import Qt, Slot
from qtpy.QtGui import QDoubleValidator
from qtpy.QtWidgets import QCheckBox, QDoubleSpinBox, QLineEdit

from MDANSE.Framework.Parameters import Float
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase


class FloatWidget(WidgetBase):
    def __init__(self, *args, parameter: Float, **kwargs):
        super().__init__(*args, parameter=parameter, **kwargs)

        self._default_value = self.get_default()

        if self.parameter.choices:
            field = QDoubleSpinBox(self._base)
            field.setMinimum(min(self.parameter.choices))
            field.setMaximum(max(self.parameter.choices))
            if len(self._configurator.choices) > 1:
                field.setSingleStep(
                    self._configurator.choices[1] - self._configurator.choices[0]
                )
            field.setValue(self._default_value)
        else:
            field = QLineEdit(self._base)
            validator = QDoubleValidator(field)

            if self.parameter.minimum is not None:
                validator.setBottom(self.parameter.minimum)
            if self.parameter.maximum is not None:
                validator.setTop(self.parameter.maximum)

            field.setValidator(validator)
            field.setText(str(self._default_value))
            field.setPlaceholderText(str(self._default_value))
            field.textChanged.connect(self.updateValue)

        if self.parameter.optional:
            self._apply_box = QCheckBox("Enable", self._base)
            self._apply_box.setTristate(False)
            self._apply_box.checkStateChanged.connect(self.toggle_widgets)
            self._apply_box.setChecked(False)
            self._layout.addWidget(self._apply_box)

        self._layout.addWidget(field)
        self._field = field
        self.default_labels()
        self.update_labels()
        self.updateValue()
        tooltip_text = self._tooltip or "A single floating-point number"
        field.setToolTip(tooltip_text)

    @Slot()
    def toggle_widgets(self):
        if not self.parameter.optional:
            return

        self._field.setEnabled(self._apply_box.checkState() == Qt.CheckState.Checked)

    def default_labels(self):
        """Each Widget should have a default tooltip and label,
        which will be set in this method, unless specific
        values are provided in the settings of the job that
        is being configured."""
        if not self._label_text:
            self._label_text = "FloatWidget"
        if not self._tooltip:
            self._tooltip = "A single floating-point number"

    def get_widget_value(self):
        """Collect the results from the input widgets and return the value."""
        if (
            self.parameter.optional
            and self._apply_box.checkState() == Qt.CheckState.Checked
        ):
            return None

        strval = self._field.text().strip()
        self._empty = not strval
        if self._empty:
            return self._default_value

        return strval
