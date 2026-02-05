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

from typing import TYPE_CHECKING

import numpy as np
from qtpy.QtCore import Qt, Slot
from qtpy.QtGui import QStandardItem, QStandardItemModel
from qtpy.QtWidgets import QAbstractScrollArea, QSizePolicy, QTableView

from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase

if TYPE_CHECKING:
    from MDANSE.Framework.Parameters import Array, Vector


class ArrayWidget(WidgetBase):
    def __init__(self, *args, parameter: Array | Vector, **kwargs):
        super().__init__(*args, parameter=parameter, **kwargs)

        self._view = QTableView(self._base)
        self._model = QStandardItemModel(self._base)

        policy = self._view.sizePolicy()
        policy.setVerticalPolicy(QSizePolicy.Policy.Minimum)
        self._view.setModel(self._model)
        self._view.setSizePolicy(policy)
        self._view.horizontalHeader().hide()
        self._view.verticalHeader().hide()
        self._view.setSizeAdjustPolicy(
            QAbstractScrollArea.SizeAdjustPolicy.AdjustToContents
        )

        if len(parameter.shape) > 2:
            raise NotImplementedError("3D data objects not supported")

        val = self.value if self.value is not None else np.zeros(self.parameter.shape)

        for row in np.atleast_2d(val):
            items = [QStandardItem(str(elem)) for elem in row]
            self._model.appendRow(items)

        self._layout.addWidget(self._view)
        self.toggle_widgets()
        self._model.itemChanged.connect(self.updateValue)

    @Slot()
    def toggle_widgets(self) -> None:
        if not self.parameter.optional:
            return

        self._view.setEnabled(self._apply_box.checkState() != Qt.CheckState.Checked)
        self.updateValue()

    def get_widget_value(self) -> np.ndarray:
        return np.array(
            [
                [
                    self._model.item(row, col).text()
                    for col in range(self._model.columnCount())
                ]
                for row in range(self._model.rowCount())
            ],
            dtype=float,
        ).reshape(self.parameter.shape)
