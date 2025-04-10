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
from qtpy.QtWidgets import (
    QDialog,
    QPushButton,
    QVBoxLayout,
    QApplication,
    QMenu,
    QLineEdit,
    QTableView,
)
from qtpy.QtCore import Signal, Slot, Qt, QSortFilterProxyModel
from qtpy.QtGui import QStandardItem, QStandardItemModel

from MDANSE.Chemistry import ATOMS_DATABASE
from MDANSE.MLogging import LOG

from MDANSE_GUI.Widgets.GeneralWidgets import InputVariable, InputDialog


class ElementView(QTableView):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setSortingEnabled(True)

    def contextMenuEvent(self, event):
        menu = QMenu(self)

        Action1 = menu.addAction("New Atom")
        Action2 = menu.addAction("Copy Atoms")
        Action3 = menu.addAction("Delete Atoms")
        Action4 = menu.addAction("New Property")

        data_model = self.parent().data_model
        row_idxs = set([idx.row() for idx in self.selectionModel().selectedIndexes()])
        atm_syms = [
            data_model.verticalHeaderItem(row_idx).text() for row_idx in row_idxs
        ]
        def_atms = ATOMS_DATABASE.default_atoms_types
        enable_delete = any([atm_sym not in def_atms for atm_sym in atm_syms])

        temp_model = self.model().sourceModel()
        if temp_model is not None:
            Action1.triggered.connect(temp_model.new_line_dialog)
            Action2.triggered.connect(temp_model.copy_rows)
            if enable_delete:
                Action3.triggered.connect(temp_model.delete_rows)
            else:
                Action3.setEnabled(False)
            Action4.triggered.connect(temp_model.new_column_dialog)
        menu.exec_(event.globalPos())


class NewElementDialog(QDialog):
    got_name = Signal(str)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        layout = QVBoxLayout(self)
        edit = QLineEdit(self)
        layout.addWidget(edit)
        self.setLayout(layout)
        self.textedit = edit
        self.button = QPushButton("Accept!", self)
        self.button.clicked.connect(self.accept)
        self.accepted.connect(self.return_value)
        layout.addWidget(self.button)

    @Slot()
    def return_value(self):
        self.got_name.emit(self.textedit.text())


class ElementModel(QStandardItemModel):
    def __init__(self, *args, element_database=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.database = element_database
        self.parseDatabase()

        self.itemChanged.connect(self.write_to_database)

    @property
    def all_row_names(self):
        return self.database.atoms

    @property
    def all_column_names(self):
        return self.database.properties

    def parseDatabase(self):
        def_columns = self.database.default_atoms_properties
        def_rows = self.database.default_atoms_types

        for entry in self.all_row_names:
            row = []
            atom_info = self.database[entry]
            for key in self.all_column_names:
                item = QStandardItem(str(atom_info[key]))
                item.setEditable(not (entry in def_rows and key in def_columns))
                try:
                    intnum = int(str(atom_info[key]))
                except ValueError:
                    try:
                        floatnum = float(atom_info[key])
                    except ValueError:
                        pass
                    except TypeError:
                        pass
                    else:
                        item.setData(floatnum, role=Qt.ItemDataRole.DisplayRole)
                except TypeError:
                    pass
                else:
                    item.setData(intnum, role=Qt.ItemDataRole.DisplayRole)
                row.append(item)
            self.appendRow(row)
        self.setHorizontalHeaderLabels(self.all_column_names)
        self.setVerticalHeaderLabels(self.all_row_names)

    @Slot("QStandardItem*")
    def write_to_database(self, item: "QStandardItem"):
        data = item.data()
        text = item.text()
        row = item.row()
        column = item.column()
        LOG.info(f"data:{data}, text:{text}, row:{row}, column:{column}")
        LOG.info(
            f"column name={self.all_column_names[column]}, row name={self.all_row_names[row]}"
        )
        self.database.set_value(
            self.all_row_names[row], self.all_column_names[column], text
        )
        self.save_changes()

    @Slot()
    def save_changes(self):
        self.database.save()

    @Slot()
    def new_line_dialog(self):
        dialog_variables = [
            InputVariable(
                input_dict={
                    "keyval": "atom_name",
                    "format": str,
                    "label": "New element name",
                    "tooltip": "Type the name of the new chemical element here.",
                    "values": [""],
                }
            )
        ]
        ne_dialog = InputDialog(fields=dialog_variables)
        ne_dialog.got_values.connect(self.add_row)
        ne_dialog.show()
        _result = ne_dialog.exec()

    @Slot()
    def new_column_dialog(self):
        dialog_variables = [
            InputVariable(
                input_dict={
                    "keyval": "property_name",
                    "format": str,
                    "label": "New property name",
                    "tooltip": "Type the name of the new property here; it will be added to the table.",
                    "values": [""],
                }
            ),
            InputVariable(
                input_dict={
                    "keyval": "property_type",
                    "format": str,
                    "label": "Type of the new property",
                    "tooltip": "One of the following: int, float, str, list",
                    "values": ["int", "float", "str", "list"],
                }
            ),
        ]
        ne_dialog = InputDialog(fields=dialog_variables)
        ne_dialog.got_values.connect(self.add_new_column)
        ne_dialog.show()
        _result = ne_dialog.exec()

    def copy_row_from_database(self, new_label: str, db_key: str):
        """Copy row data from one database entry to another and update
        the table. Saves changes to the database.

        Parameters
        ----------
        new_label : str
            The new label the results are copied to.
        db_key : str
            The key of the data the new_labels data will be copied from.
        """
        row = []
        for key in self.all_column_names:
            new_value = self.database.get_value(db_key, key)
            self.database.set_value(new_label, key, new_value)
            item = QStandardItem(str(new_value))
            row.append(item)
        self.appendRow(row)
        self.setVerticalHeaderItem(self.rowCount() - 1, QStandardItem(str(new_label)))
        LOG.info(f"self.all_row_names has length: {len(self.all_row_names)}")

    @Slot()
    def copy_rows(self):
        """Update the database and table with a copied atoms."""
        view = self.parent().viewer
        row_idx = set([idx.row() for idx in view.selectionModel().selectedIndexes()])
        for idx in row_idx:
            atm_sym = self.verticalHeaderItem(idx).text()
            atm_sym_copy = atm_sym + " (copy)"
            while True:
                if atm_sym_copy not in self.database.atoms:
                    self.database.add_atom(atm_sym_copy)
                    self.copy_row_from_database(atm_sym_copy, atm_sym)
                    break
                else:
                    atm_sym_copy += " (copy)"
        self.save_changes()

    @Slot(dict)
    def add_row(self, input_variables: dict):
        """Add a new line to the table from the database.

        Parameters
        ----------
        input_variables : dict
            Variables used to create the new entry.
        """
        try:
            new_label = input_variables["atom_name"]
        except KeyError:
            return None
        if new_label not in self.database.atoms:
            self.database.add_atom(new_label)
            self.copy_row_from_database(new_label, new_label)
        self.save_changes()

    @Slot()
    def delete_rows(self):
        """Delete a rows from the table and update the database."""
        view = self.parent().viewer
        row_idx = list(
            set([idx.row() for idx in view.selectionModel().selectedIndexes()])
        )
        row_idx.sort(reverse=True)
        def_atms = ATOMS_DATABASE.default_atoms_types
        for idx in row_idx:
            atm_sym = self.verticalHeaderItem(idx).text()
            if atm_sym not in def_atms:
                view.model().removeRow(idx)
                self.database.remove_atom(atm_sym)
        self.save_changes()

    @Slot(dict)
    def add_new_column(self, input_variables: dict):
        try:
            new_label = input_variables["property_name"]
        except KeyError:
            return None
        try:
            new_type = input_variables["property_type"]
        except KeyError:
            return None
        if new_label not in self.database.atoms:
            self.database.add_property(new_label, new_type)
            column = []
            for key in self.all_row_names:
                new_value = self.database.get_value(key, new_label)
                self.database.set_value(key, new_label, new_value)
                item = QStandardItem(str(new_value))
                column.append(item)
            self.all_column_names.append(new_label)
            self.appendColumn(column)
            self.setHorizontalHeaderItem(
                self.columnCount() - 1, QStandardItem(str(new_label))
            )
            LOG.info(f"self.all_column_names has length: {len(self.all_column_names)}")


class ElementsDatabaseEditor(QDialog):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.setWindowTitle("MDANSE Chemical Elements Database Editor")

        layout = QVBoxLayout(self)

        self.setLayout(layout)

        self.viewer = ElementView(self)
        layout.addWidget(self.viewer)

        self.data_model = ElementModel(self, element_database=ATOMS_DATABASE)

        self.proxy_model = QSortFilterProxyModel(self)

        self.proxy_model.setSourceModel(self.data_model)
        self.viewer.setModel(self.proxy_model)
        self.resize(800, 600)


if __name__ == "__main__":
    import sys

    app = QApplication(sys.argv)
    root = ElementsDatabaseEditor()
    root.show()
    app.exec()
