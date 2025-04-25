import numpy as np
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase
from qtpy.QtCore import Signal, Slot
from qtpy.QtGui import QStandardItem, QStandardItemModel
from qtpy.QtWidgets import QPushButton, QSizePolicy, QSpinBox, QTableView


class VectorListTable(QStandardItemModel):
    """Vector list container."""
    rowAdded = Signal(int)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.n_vecs = 0

    @Slot(int)
    def set_rows(self, n: int):
        if n < self.n_vecs:  # Remove
            self.setRowCount(n)

        for _ in range(n - self.n_vecs):  # Add
            self.appendRow([QStandardItem("0") for _ in "xyz"])
        self.n_vecs = n


    def get_widget_value(self):
        return [[self.item(rownum, dim).text() for dim in (0, 1, 2)]
                for rownum in range(self.rowCount())]

    @Slot()
    def add_row(self):
        self.appendRow([QStandardItem("0") for _ in "xyz"])
        self.n_vecs += 1
        self.rowAdded.emit(self.n_vecs)

class VectorListWidget(WidgetBase):
    def __init__(self, n_rows: int, *args, **kwargs):
        kwargs["layout_type"] = "QVBoxLayout"
        super().__init__(*args, **kwargs)

        self._table = VectorListTable()

        self._view = QTableView(self._base)
        self._layout.addWidget(self._view)
        self._view.setModel(self._table)
        policy = self._view.sizePolicy()
        policy.setVerticalPolicy(QSizePolicy.Policy.Minimum)
        self._view.setSizePolicy(policy)
        self._view.horizontalHeader().hide()

        self._rows = QSpinBox(self._base)
        self._rows.setRange(0, 1000)
        self._layout.addWidget(self._rows)

        self._rows.valueChanged.connect(self._table.set_rows)
        self._rows.setValue(n_rows)
        self._table.rowAdded.connect(self._rows.setValue)

        self._add_row = QPushButton("Add vector")
        self._add_row.clicked.connect(self._table.add_row)
        self._layout.addWidget(self._add_row)

        self._table.itemChanged.connect(self.updateValue)
        self._table.rowAdded.connect(self.updateValue)

    def updateValue(self):
        current_value = self._table.get_widget_value()

        try:
            self._configurator.configure(current_value)
        except Exception:
            self.mark_error(
                "COULD NOT SET THIS VALUE - you may need to change the values in other widgets"
            )

        self.value_changed.emit()

        if self._configurator.valid:
            self.clear_error()
            self.value_updated.emit()
        else:
            self.mark_error(self._configurator.error_status)
