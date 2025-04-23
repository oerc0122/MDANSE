import numpy as np

from qtpy.QtWidgets import QWidget, QSizePolicy, QVBoxLayout, QLabel, QDialog
from qtpy.Qt3DRender import QDirectionalLight
from qtpy.QtGui import QPixmap, QFont
from qtpy.QtCore import Qt as _Qt
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
        large_molecule = self.atom_database.chemical_system.rdkit_mol
        draw_options = draw.MolDrawOptions()
        draw_options.addAtomIndices = True
        pil_image = draw.MolToImage(large_molecule, options=draw_options)
        qt_image = ImageQt(pil_image)
        pixmap = QPixmap.fromImage(qt_image)
        self.image_label.setPixmap(pixmap)

    def show_coordinates(self):
        mass = []
        coords = []
        coordinates = self.mol_info["atom_coordinates"]
        _indices = self.mol_info["atom_indices"]
        atom_symbols = self.mol_info["atom_symbols"]
        bonds = self.mol_info["bond_list"]
        for at_number in range(len(coordinates)):
            x, y, z = coordinates[at_number]
            x, y, z = (20 * x - 10, 20 * y - 10, 20 * z - 10)
            symbol = atom_symbols[at_number]
            colour = self.atom_database.get_atom_property(symbol, "color")
            radius = self.atom_database.get_atom_property(symbol, "covalent_radius")
            mass.append(self.atom_database.get_atom_property(symbol, "atomic_weight"))
            coords.append(coordinates[at_number])
        mass = np.array(mass)
        coords = np.array(coords)
        com = np.einsum("i,ik->k", mass, coords) / np.sum(mass)
        x_com, y_com, z_com = 20 * com - 10
