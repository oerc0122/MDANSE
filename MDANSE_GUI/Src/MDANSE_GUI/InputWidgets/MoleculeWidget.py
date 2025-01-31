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

import numpy as np

from qtpy.QtCore import Slot
from qtpy.QtWidgets import QComboBox, QPushButton
from MDANSE_GUI.InputWidgets.WidgetBase import WidgetBase
from MDANSE_GUI.InputWidgets.MoleculePreviewWidget import MoleculePreviewWidget


class MoleculeWidget(WidgetBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        configurator = kwargs.get("configurator", None)
        trajectory_configurator = kwargs.get("trajectory_configurator", None)
        default_option = ""
        if trajectory_configurator is not None:
            option_list = trajectory_configurator[
                "instance"
            ].chemical_system.unique_molecules()
            if len(option_list) > 0:
                default_option = option_list[0]
        else:
            if configurator is None:
                option_list = kwargs.get("choices", [])
            else:
                option_list = configurator.choices
                default_option = configurator.default
        traj_config = self._configurator._configurable[
            self._configurator._dependencies["trajectory"]
        ]
        hdf_traj = traj_config["hdf_trajectory"]
        unique_molecules = hdf_traj.chemical_system.unique_molecules()
        traj_bond_list = hdf_traj.chemical_system._bonds
        self.atom_database = hdf_traj
        self.mol_dict = {}
        coords_0 = hdf_traj.trajectory.coordinates(0)
        for mol_name in unique_molecules:
            no_of_molecules = len(hdf_traj.chemical_system._clusters[mol_name])
            atom_indices = hdf_traj.chemical_system._clusters[mol_name][0]
            atom_symbols = [
                hdf_traj.chemical_system.atom_list[index] for index in atom_indices
            ]
            coordinates = coords_0[atom_indices]
            unique_atoms, atom_counts = np.unique(atom_symbols, return_counts=True)
            atom_counts = {
                unique_atoms[n]: atom_counts[n] for n in range(len(unique_atoms))
            }
            bonds = [
                (coords_0[bond[0]], coords_0[bond[1]])
                for bond in traj_bond_list
                if bond[0] in atom_indices or bond[1] in atom_indices
            ]
            self.mol_dict[mol_name] = {
                "no_of_molecules": no_of_molecules,
                "atom_coordinates": coordinates,
                "atom_number": atom_counts,
                "atom_indices": atom_indices,
                "atom_symbols": atom_symbols,
                "bond_list": bonds,
            }

        self.field = QComboBox(self._base)
        self.field.addItems(option_list)
        self.field.setCurrentText(default_option)
        self.selected_name = self.field.currentText()
        self.selected_mol = self.mol_dict[self.selected_name]
        self.field.currentTextChanged.connect(self.updateValue)
        self.field.currentTextChanged.connect(self.molecule_changed)
        button = QPushButton(self._base)
        button.setText("Molecule Preview")
        button.clicked.connect(self.button_clicked)
        if self._tooltip:
            tooltip_text = self._tooltip
        else:
            tooltip_text = (
                "A single option can be picked out of all the options listed."
            )
        self.field.setToolTip(tooltip_text)
        self._field = self.field
        self._layout.addWidget(self.field)
        self._layout.addWidget(button)
        self._configurator = configurator
        self.default_labels()
        self.update_labels()
        self.updateValue()

    @Slot()
    def molecule_changed(self):
        """
        Change molecule preview and molecule information
        """
        self.selected_name = self.field.currentText()
        self.selected_mol = self.mol_dict[self.selected_name]
        self.window = MoleculePreviewWidget(
            self._base, self.selected_mol, self.selected_name, self.atom_database
        )

    @Slot()
    def button_clicked(self):
        """
        Opens a window that shows a preview of selected molecule
        """
        self.window = MoleculePreviewWidget(
            self._base, self.selected_mol, self.selected_name, self.atom_database
        )
        if self.window.isVisible():
            self.window.close()
        else:
            self.window.show()

    def configure_using_default(self):
        """This is too complex to have a default value"""

    def default_labels(self):
        """Each Widget should have a default tooltip and label,
        which will be set in this method, unless specific
        values are provided in the settings of the job that
        is being configured."""
        if self._label_text == "":
            self._label_text = "ComboWidget"
        if self._tooltip == "":
            self._tooltip = "You only have one option. Choose wisely."

    def get_widget_value(self):
        return self._field.currentText()
