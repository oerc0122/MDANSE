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

from qtpy.QtWidgets import QVBoxLayout, QLabel, QDialog
from qtpy.QtGui import QPixmap, QFont
import rdkit.Chem as chem
import rdkit.Chem.Draw as draw
from PIL.ImageQt import ImageQt


class MoleculePreviewWidget(QDialog):
    def __init__(self, parent, molecule_information, molecule_name, atom_database):
        super().__init__(parent)
        self.setWindowTitle("Molecule Preview")
        self.resize(800, 600)
        self.axes = []
        layout = QVBoxLayout(self)
        self.setLayout(layout)
        self.mol_name = molecule_name
        self.mol_info = molecule_information
        self.atom_database = atom_database
        self.image_label = QLabel(self)
        self.text_label = QLabel(self)
        font = QFont("Arial", 12)
        self.text_label.setFont(font)
        self.text_label.setWordWrap(True)
        layout.addWidget(self.image_label)
        layout.addWidget(self.text_label)
        self.update_text_label()
        self.show_formula()

    def update_text_label(self):
        info_text = f"Molecule name: {self.mol_name}\n"
        for key, value in self.mol_info["atom_number"].items():
            info_text += f"Number of {key} atoms: {value}\n"
        info_text += f"Number of such molecules in trajectory: {self.mol_info['no_of_molecules']}\n"
        self.text_label.setText(info_text)

    def show_formula(self):
        if len(self.mol_info["atom_indices"]) > 300:
            self.image_label.clear()
            self.image_label.setText("Molecule is too large for preview")
            return
        self.image_label.setText("")
        large_molecule = self.atom_database.chemical_system.rdkit_mol
        submolecule = chem.RWMol()
        mapping = {}
        for index in self.mol_info["atom_indices"]:
            new_index = submolecule.AddAtom(large_molecule.GetAtomWithIdx(index))
            mapping[index] = new_index
        for bond in self.mol_info["bond_list"]:
            new_pair = mapping[bond[0]], mapping[bond[1]]
            submolecule.AddBond(new_pair[0], new_pair[1])
        draw_options = draw.MolDrawOptions()
        draw_options.addAtomIndices = True
        pil_image = draw.MolToImage(submolecule, size=(600, 600), options=draw_options)
        qt_image = ImageQt(pil_image)
        pixmap = QPixmap.fromImage(qt_image)
        self.image_label.setPixmap(pixmap)
