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

from qtpy.QtWidgets import QDoubleSpinBox, QLabel, QSpinBox

from MDANSE.Framework.Parameters import Range
from MDANSE.MLogging import LOG
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase


class RangeWidget(WidgetBase):
    def __init__(self, *args, parameter: Range, **kwargs):
        super().__init__(
            *args, parameter=parameter, layout_type="QGridLayout", **kwargs
        )

        box_type = QSpinBox if self.parameter.dtype is int else QDoubleSpinBox

        self._labels = []
        self._fields = []
        labels = ("Start", "Stop", "Step")

        tooltip_text = (
            self._tooltip or "Values to be used, given as (First, Last, StepSize)"
        )

        for i, label_name in enumerate(labels):
            label = QLabel(label_name, self._base)
            field = box_type(self._base)

            field.setToolTip(tooltip_text)
            field.valueChanged.connect(self.updateValue)

            self._labels.append(label)
            self._fields.append(field)

            self._layout.addWidget(label, 0, 2 * i)
            self._layout.addWidget(field, 0, 2 * i + 1)

        self.trajectory_changed()
        self.default_labels()
        self.update_labels()
        self.updateValue()

    def default_labels(self):
        """Each Widget should have a default tooltip and label,
        which will be set in this method, unless specific
        values are provided in the settings of the job that
        is being configured."""
        if not self._label_text:
            self._label_text = "RangeWidget"
        if not self._tooltip:
            self._tooltip = "Values to be used, given as (First, Last, StepSize)"

    def trajectory_changed(self) -> None:
        if not self.ranges:
            for field in self._fields:
                field.setEnabled(True)
                field.setRange(None, None)

            return


        mini, maxi = self.ranges

        if self.default != "N/A":
            defaults = self.default
        elif self.ranges != (None, None) and self.parameter.dtype is not int:
            step = round((maxi - mini) / 100, 2)
            defaults = (mini, maxi, step)
        elif self.ranges != (None, None):
            step = 1
            defaults = (mini, maxi, step)
        else:
            defaults = (0, -1, 1)

        # defaults = [
        #     default if default is not None else std_default
        #     for default, std_default in zip((mini, maxi, step), (0, -1, 1), strict=True)
        # ]

        # if self.parameter.include_last:
        #     maxi += defaults[2]

        for field, val in zip(self._fields, defaults, strict=True):
            field.setValue(val)
            field.setEnabled(True)
            field.setSingleStep(step)
            field.setRange(mini, maxi)

    def get_widget_value(self):
        values = [field.value() for field in self._fields]
        if self.parameter.include_last:
            values[1] -= values[2]  # Include last

        return tuple(values)
