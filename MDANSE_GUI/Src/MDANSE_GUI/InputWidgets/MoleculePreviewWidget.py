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
from qtpy.QtGui import QPixmap, QImage, QFont
import rdkit.Chem as chem
import rdkit.Chem.Draw as draw
import rdkit.Chem.AllChem as allchem
from PIL.ImageQt import ImageQt


class MoleculePreviewWidget(QDialog):
    def __init__(self, parent, molecule_information, molecule_name, atom_database):
        super().__init__(parent)
        self.setWindowTitle("Molecule Preview")
        self.resize(800, 600)
        self.axes = []
        self.bond_mapping = {}
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

    def prepare_rdkit_molecule(self):
        if len(self.mol_info["atom_indices"]) > 300:
            self.image_label.clear()
            self.image_label.setText("Molecule is too large for preview")
            return
        self.image_label.setText("")
        large_molecule = self.atom_database.chemical_system.rdkit_mol
        submolecule = chem.RWMol()
        mapping = {}
        self.bond_mapping = {}
        for index in self.mol_info["atom_indices"]:
            new_index = submolecule.AddAtom(large_molecule.GetAtomWithIdx(index))
            mapping[index] = new_index
        for bond in self.mol_info["bond_list"]:
            new_pair = mapping[bond[0]], mapping[bond[1]]
            bond_index = submolecule.AddBond(new_pair[0], new_pair[1])
            self.bond_mapping[bond_index] = new_pair
        allchem.Compute2DCoords(submolecule)
        return submolecule

    def show_formula(self):
        submolecule = self.prepare_rdkit_molecule()
        draw_options = draw.MolDrawOptions()
        draw_options.addAtomIndices = True
        pil_image = draw.MolToImage(submolecule, size=(600, 600), options=draw_options)
        qt_image = ImageQt(pil_image)
        pixmap = QPixmap.fromImage(qt_image)
        self.image_label.setPixmap(pixmap)

    def filter_bonds(self, selected_atoms: list[int]) -> list[int]:
        return [
            key
            for key, pair in self.bond_mapping.items()
            if pair[0] in selected_atoms and pair[1] in selected_atoms
        ]

    def show_formula_with_selection(self, selected_atoms: list[int]):
        submolecule = self.prepare_rdkit_molecule()
        optional_bonds = self.filter_bonds(selected_atoms)
        driver = draw.rdMolDraw2D.MolDraw2DCairo(600, 600)
        draw_options = draw.MolDrawOptions()
        draw_options.addAtomIndices = True
        driver.SetDrawOptions(draw_options)
        draw.rdMolDraw2D.PrepareAndDrawMolecule(
            driver,
            submolecule,
            highlightAtoms=selected_atoms,
            highlightBonds=optional_bonds,
        )
        png_data = driver.GetDrawingText()
        image = QImage.fromData(png_data)
        pixmap = QPixmap.fromImage(image)
        self.image_label.setPixmap(pixmap)
