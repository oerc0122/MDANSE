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

from collections.abc import Sequence
from enum import Enum, EnumMeta

from qtpy.QtWidgets import QComboBox

from MDANSE.Framework.Parameters.Choices import Choice
from MDANSE.Framework.Parameters.Parameters import CustomChoices
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase
from MDANSE_GUI.Widgets.DefaultCombobox import highlight_default_value


class ComboWidget(WidgetBase):
    def __init__(self, *args, parameter: Choice, choices: Sequence = (), **kwargs):
        super().__init__(*args, parameter=parameter, **kwargs)

        default = self.default

        option_list = choices or self.choices

        if self._tooltip:
            tooltip_text = self._tooltip
        else:
            tooltip_text = (
                "A single option can be picked out of all the options listed."
            )

        self._field = QComboBox(self._base)

        if option_list:
            self._field.addItems(sorted(map(str, option_list)))
            self._field.setCurrentText(str(default))
            highlight_default_value(self._field)
        else:
            self.get_choices()

        self._field.currentTextChanged.connect(self.updateValue)
        self._field.setToolTip(tooltip_text)
        self._layout.addWidget(self._field)

        self.default_labels()
        self.update_labels()
        self.updateValue()

    def get_choices(self):
        self._field.clear()

        if any(self.parameter._bad_deps(self._configurable)):
            self.mark_error("Invalid dependencies")
            self._field.setCurrentText("Unavailable")
            self._field.setEnabled(False)
            return None

        self._field.addItems(sorted(self.choices))

        if self.default != "N/A":
            self._field.setCurrentText(self.default)
            highlight_default_value(self._field)
        else:
            self._field.setCurrentIndex(0)

    def default_labels(self):
        """Each Widget should have a default tooltip and label,
        which will be set in this method, unless specific
        values are provided in the settings of the job that
        is being configured."""
        if not self._label_text:
            self._label_text = "ComboWidget"
        if not self._tooltip:
            self._tooltip = "You only have one option. Choose wisely."

    def get_widget_value(self):
        text = self._field.currentText()
        if text == "None":
            return None
        return text
