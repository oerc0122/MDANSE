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
from qtpy.QtCore import Qt, Slot
from qtpy.QtWidgets import (
    QLineEdit,
    QPushButton,
    QDialog,
    QCheckBox,
    QVBoxLayout,
    QHBoxLayout,
    QGroupBox,
    QLabel,
    QPlainTextEdit,
    QWidget,
)
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase
from MDANSE_GUI.Tabs.Visualisers.View3D import View3D
from MDANSE_GUI.MolecularViewer.MolecularViewer import MolecularViewerWithPicking
from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from .CheckableComboBox import CheckableComboBox


class SelectionHelper(QDialog):
    """Generates a string that specifies the atom selection.

    Attributes
    ----------
    _helper_title : str
        The title of the helper dialog window.
    _cbox_text : dict
        The dictionary that maps the selector settings to text used in
        the helper dialog.
    """

    _helper_title = "Atom selection helper"
    _cbox_text = {
        "all": "All atoms (excl. dummy atoms):",
        "dummy": "All dummy atoms:",
        "hs_on_heteroatom": "Hs on heteroatoms:",
        "primary_amine": "Primary amine groups:",
        "hydroxy": "Hydroxy groups:",
        "methyl": "Methyl groups:",
        "phosphate": "Phosphate groups:",
        "sulphate": "Sulphate groups:",
        "thiol": "Thiol groups:",
        "water": "Water molecules:",
        "hs_on_element": "Hs on elements:",
        "element": "Elements:",
        "name": "Atom name:",
        "fullname": "Atom fullname:",
        "index": "Indexes:",
    }

    def __init__(
        self,
        selector,
        traj_data: tuple[str, HDFTrajectoryInputData],
        field: QLineEdit,
        parent,
        *args,
        **kwargs,
    ):
        """
        Parameters
        ----------
        selector : Selector
            The MDANSE selector initialized with the current chemical
            system.
        traj_data : tuple[str, HDFTrajectoryInputData]
            A tuple of the trajectory data used to load the 3D viewer.
        field : QLineEdit
            The QLineEdit field that will need to be updated when
            applying the setting.
        """
        super().__init__(parent, *args, **kwargs)
        self.setWindowTitle(self._helper_title)

        self.selector = selector
        self._field = field
        self.settings = self.selector.settings
        self.atm_full_names = self.selector.system.name_list

        self.selection_textbox = QPlainTextEdit()
        self.selection_textbox.setReadOnly(True)

        mol_view = MolecularViewerWithPicking()
        mol_view.picked_atoms_changed.connect(self.update_from_3d_view)
        self.view_3d = View3D(mol_view)
        self.view_3d.update_panel(traj_data)

        layouts = self.create_layouts()

        bottom = QHBoxLayout()
        for button in self.create_buttons():
            bottom.addWidget(button)

        layouts[-1].addLayout(bottom)

        helper_layout = QHBoxLayout()
        for layout in layouts:
            helper_layout.addLayout(layout)

        self.setLayout(helper_layout)
        self.update_others()

        self.all_selection = True
        self.selected = set([])

    def closeEvent(self, a0):
        """Hide the window instead of closing. Some issues occur in the
        3D viewer when it is closed and then reopened.
        """
        a0.ignore()
        self.hide()

    def create_buttons(self) -> list[QPushButton]:
        """
        Returns
        -------
        list[QPushButton]
            List of push buttons to add to the last layout from
            create_layouts.
        """
        apply = QPushButton("Use Setting")
        reset = QPushButton("Reset")
        close = QPushButton("Close")
        apply.clicked.connect(self.apply)
        reset.clicked.connect(self.reset)
        close.clicked.connect(self.close)
        return [apply, reset, close]

    def create_layouts(self) -> list[QVBoxLayout]:
        """
        Returns
        -------
        list[QVBoxLayout]
            List of QVBoxLayout to add to the helper layout.
        """
        layout_3d = QVBoxLayout()
        layout_3d.addWidget(self.view_3d)

        left = QVBoxLayout()
        for widget in self.left_widgets():
            left.addWidget(widget)

        right = QVBoxLayout()
        for widget in self.right_widgets():
            right.addWidget(widget)

        return [layout_3d, left, right]

    def right_widgets(self) -> list[QWidget]:
        """
        Returns
        -------
        list[QWidget]
            List of QWidgets to add to the right layout from
            create_layouts.
        """
        return [self.selection_textbox]

    def left_widgets(self) -> list[QWidget]:
        """
        Returns
        -------
        list[QWidget]
            List of QWidgets to add to the left layout from
            create_layouts.
        """
        match_exists = self.selector.match_exists

        select = QGroupBox("selection")
        select_layout = QVBoxLayout()

        self.check_boxes = []
        self.combo_boxes = []

        for k, v in self.settings.items():

            if isinstance(v, bool):
                check_layout = QHBoxLayout()
                checkbox = QCheckBox()
                checkbox.setChecked(v)
                checkbox.setLayoutDirection(Qt.RightToLeft)
                label = QLabel(self._cbox_text[k])
                checkbox.setObjectName(k)
                checkbox.stateChanged.connect(self.update_others)
                if not match_exists[k]:
                    checkbox.setEnabled(False)
                    label.setStyleSheet("color: grey;")
                self.check_boxes.append(checkbox)
                check_layout.addWidget(label)
                check_layout.addWidget(checkbox)
                select_layout.addLayout(check_layout)

            elif isinstance(v, dict):
                combo_layout = QHBoxLayout()
                combo = CheckableComboBox()
                items = [str(i) for i in v.keys() if match_exists[k][i]]
                # we blocksignals here as there can be some
                # performance issues with a large number of items
                combo.model().blockSignals(True)
                combo.addItems(items)
                combo.model().blockSignals(False)
                combo.setObjectName(k)
                combo.model().dataChanged.connect(self.update_others)
                label = QLabel(self._cbox_text[k])
                if len(items) == 0:
                    combo.setEnabled(False)
                    label.setStyleSheet("color: grey;")
                self.combo_boxes.append(combo)
                combo_layout.addWidget(label)
                combo_layout.addWidget(combo)
                select_layout.addLayout(combo_layout)

        invert_layout = QHBoxLayout()
        label = QLabel("Invert selection:")
        apply = QPushButton("Apply")
        apply.clicked.connect(self.invert_selection)
        invert_layout.addWidget(label)
        invert_layout.addWidget(apply)
        select_layout.addLayout(invert_layout)

        select.setLayout(select_layout)
        return [select]

    def update_others(self) -> None:
        """Using the checkbox and combobox widgets: update the settings,
        get the selection and update the textedit box with details of
        the current selection and the 3d view to match the selection.
        """
        for check_box in self.check_boxes:
            self.settings[check_box.objectName()] = check_box.isChecked()
        for combo_box in self.combo_boxes:
            for i in range(combo_box.n_items):
                txt = combo_box.text[i]
                if combo_box.objectName() == "index":
                    key = int(txt)
                else:
                    key = txt
                self.settings[combo_box.objectName()][key] = combo_box.checked[i]

        self.selector.update_settings(self.settings)
        self.selected = self.selector.get_idxs()
        self.view_3d._viewer.change_picked(self.selected)
        self.update_selection_textbox()

    def update_from_3d_view(self, selection: set[int]) -> None:
        """A selection/deselection was made in the 3d view, update the
        check_boxes, combo_boxes and textbox.

        Parameters
        ----------
        selection : set[int]
            Selection indexes from the 3d view.
        """
        self.selector.update_with_idxs(selection)
        self.settings = self.selector.settings
        self.update_selection_widgets()
        self.selected = self.selector.get_idxs()
        self.update_selection_textbox()

    def invert_selection(self):
        """Inverts the selection."""
        self.selected = self.selector.all_idxs - self.selected
        self.selector.update_with_idxs(self.selected)
        self.settings = self.selector.settings
        self.update_selection_widgets()
        self.view_3d._viewer.change_picked(self.selected)
        self.update_selection_textbox()

    def update_selection_widgets(self) -> None:
        """Updates the selection widgets so that it matches the full
        setting.
        """
        for check_box in self.check_boxes:
            check_box.blockSignals(True)
            if self.settings[check_box.objectName()]:
                check_box.setCheckState(Qt.Checked)
            else:
                check_box.setCheckState(Qt.Unchecked)
            check_box.blockSignals(False)
        for combo_box in self.combo_boxes:
            combo_box.model().blockSignals(True)
            for i in range(combo_box.n_items):
                txt = combo_box.text[i]
                if combo_box.objectName() == "index":
                    key = int(txt)
                else:
                    key = txt
                combo_box.set_item_checked_state(
                    i, self.settings[combo_box.objectName()][key]
                )
            combo_box.update_all_selected()
            combo_box.update_line_edit()
            combo_box.model().blockSignals(False)

    def update_selection_textbox(self) -> None:
        """Update the selection textbox."""
        num_sel = len(self.selected)
        text = [f"Number of atoms selected:\n{num_sel}\n\nSelected atoms:\n"]
        for idx in self.selected:
            text.append(f"{idx}  ({self.atm_full_names[idx]})\n")
        self.selection_textbox.setPlainText("".join(text))

    def apply(self) -> None:
        """Set the field of the AtomSelectionWidget to the currently
        chosen setting in this widget.
        """
        self.selector.update_settings(self.settings)
        self._field.setText(self.selector.settings_to_json())

    def reset(self) -> None:
        """Resets the helper to the default state."""
        self.selector.reset_settings()
        self.selector.settings["all"] = self.all_selection
        self.settings = self.selector.settings
        self.update_selection_widgets()
        self.selected = self.selector.get_idxs()
        self.view_3d._viewer.change_picked(self.selected)
        self.update_selection_textbox()


class AtomSelectionWidget(WidgetBase):
    """The atoms selection widget."""

    _push_button_text = "Atom selection helper"
    _default_value = '{"all": true}'
    _tooltip_text = "Specify which atoms will be used in the analysis. The input is a JSON string, and can be created using the helper dialog."

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._value = self._default_value
        self._field = QLineEdit(self._default_value, self._base)
        self._field.setPlaceholderText(self._default_value)
        self._field.setMaxLength(2147483647)  # set to the largest possible
        self._field.textChanged.connect(self.updateValue)
        traj_config = self._configurator._configurable[
            self._configurator._dependencies["trajectory"]
        ]
        traj_filename = traj_config["filename"]
        hdf_traj = traj_config["hdf_trajectory"]
        self.helper = self.create_helper((traj_filename, hdf_traj))
        helper_button = QPushButton(self._push_button_text, self._base)
        helper_button.clicked.connect(self.helper_dialog)
        self._layout.addWidget(self._field)
        self._layout.addWidget(helper_button)
        self.update_labels()
        self.updateValue()
        self._field.setToolTip(self._tooltip_text)

    def create_helper(
        self, traj_data: tuple[str, HDFTrajectoryInputData]
    ) -> SelectionHelper:
        """
        Parameters
        ----------
        traj_data : tuple[str, HDFTrajectoryInputData]
            A tuple of the trajectory data used to load the 3D viewer.

        Returns
        -------
        SelectionHelper
            Create and return the selection helper QDialog.
        """
        selector = self._configurator.get_selector()
        return SelectionHelper(selector, traj_data, self._field, self._base)

    @Slot()
    def helper_dialog(self) -> None:
        """Opens the helper dialog."""
        if self.helper.isVisible():
            self.helper.close()
        else:
            self.helper.show()

    def get_widget_value(self) -> str:
        """
        Returns
        -------
        str
            The JSON selector setting.
        """
        selection_string = self._field.text()
        if len(selection_string) < 1:
            self._empty = True
            return self._default_value
        else:
            self._empty = False
        return selection_string
