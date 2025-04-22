import numpy as np

from qtpy.QtWidgets import QWidget, QSizePolicy, QVBoxLayout, QLabel, QDialog
from qtpy.Qt3DRender import QDirectionalLight
from qtpy.QtGui import QColor, QVector3D, QQuaternion, QFont
from qtpy.Qt3DExtras import (
    QPhongMaterial,
    QCylinderMesh,
    QSphereMesh,
    Qt3DWindow,
    QOrbitCameraController,
)
from qtpy.QtCore import Qt as _Qt
from qtpy.Qt3DCore import QEntity, QTransform


class MoleculePreviewWidget(QDialog):
    def __init__(self, parent, molecule_information, molecule_name, atom_database):
        super().__init__(parent)
        self.setWindowTitle("Molecule Preview")
        self.resize(800, 600)
        self.view = Qt3DWindow()
        self.view.defaultFrameGraph().setClearColor(
            QColor(0x4D4D4F)
        )  # molecular viewer mdansechemistry atoms.json
        container = QWidget.createWindowContainer(
            self.view
        )  # from mdanse chemistry atoms database
        screenSize = self.view.screen().size()
        container.setMinimumSize(200, 100)
        container.setMaximumSize(screenSize)
        container.setSizePolicy(
            QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.MinimumExpanding
        )
        layout = QVBoxLayout()
        layout.addWidget(container)
        self.rootEntity = QEntity()
        self.cuboidTransform = QTransform()
        self.axes = []
        mass = []
        coords = []
        info_text = f"Molecule name: {molecule_name}\n"
        for key, value in molecule_information["atom_number"].items():
            info_text += f"Number of {key} atoms: {value}\n"

        info_text += f"Number of such molecules in trajectory: {molecule_information['no_of_molecules']}\n"

        coordinates = molecule_information["atom_coordinates"]
        _indices = molecule_information["atom_indices"]
        atom_symbols = molecule_information["atom_symbols"]
        bonds = molecule_information["bond_list"]
        for at_number in range(len(coordinates)):
            x, y, z = coordinates[at_number]
            x, y, z = (20 * x - 10, 20 * y - 10, 20 * z - 10)
            symbol = atom_symbols[at_number]
            colour = atom_database.get_atom_property(symbol, "color")
            radius = atom_database.get_atom_property(symbol, "covalent_radius")
            mass.append(atom_database.get_atom_property(symbol, "atomic_weight"))
            coords.append(coordinates[at_number])
            r, g, b = [int(x) for x in colour.split(";")]
            colour = QColor(r, g, b)
            m_sphereEntity = QEntity(self.rootEntity)
            sphereMesh = QSphereMesh()
            sphereMesh.setRings(20)
            sphereMesh.setSlices(20)
            sphereMesh.setRadius(radius * 10)
            sphereTransform = QTransform()
            sphereTransform.setScale(1.0)
            sphereTransform.setTranslation(QVector3D(x, y, z))
            sphereMaterial = QPhongMaterial()
            sphereMaterial.setDiffuse(colour)
            sphereMaterial.setAmbient(colour)  # Set ambient to the same as diffuse.
            # sphereMaterial.setSpecular(QColor(0, 0, 0))  # Eliminate specular reflection.
            # sphereMaterial.setShininess(0.0)  # Eliminate shininess.
            m_sphereEntity.addComponent(sphereMesh)
            m_sphereEntity.addComponent(sphereMaterial)
            m_sphereEntity.addComponent(sphereTransform)

        _atom_information = molecule_information["atom_information"]
        for bond in bonds:
            coord1, coord2 = bond[0], bond[1]
            coord1 = (20 * coord1[0] - 10, 20 * coord1[1] - 10, 20 * coord1[2] - 10)
            coord2 = (20 * coord2[0] - 10, 20 * coord2[1] - 10, 20 * coord2[2] - 10)
            direction = QVector3D(
                coord2[0] - coord1[0], coord2[1] - coord1[1], coord2[2] - coord1[2]
            )
            length = direction.length()
            direction.normalize()
            # Compute rotation
            up_vector = QVector3D(0, 1, 0)
            axis = QVector3D.crossProduct(up_vector, direction)
            angle = float(
                np.degrees(np.arccos(QVector3D.dotProduct(up_vector, direction)))
            )

            # Create cylinder mesh
            cylinder_mesh = QCylinderMesh()
            cylinder_mesh.setRadius(radius)
            cylinder_mesh.setLength(length)

            # Create material
            material = QPhongMaterial()
            # material.setDiffuse(QColor(color))

            # Set transformation
            transform = QTransform()
            transform.setTranslation(QVector3D(*coord1) + direction * length / 2)
            transform.setRotation(QQuaternion.fromAxisAndAngle(axis, angle))

            # Create entity
            entity = QEntity(self.rootEntity)
            entity.addComponent(cylinder_mesh)
            entity.addComponent(material)
            entity.addComponent(transform)

        info_label = QLabel(info_text)
        font = QFont("Arial", 12)
        info_label.setFont(font)
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        self.setLayout(layout)
        mass = np.array(mass)
        coords = np.array(coords)
        com = np.einsum("i,ik->k", mass, coords) / np.sum(mass)
        x_com, y_com, z_com = 20 * com - 10
        # Camera
        self.camera = self.view.camera()
        self.camera.lens().setPerspectiveProjection(45.0, 16.0 / 9.0, 0.1, 1000.0)
        self.camera.setPosition(QVector3D(x_com + 5, y_com + 5, z_com + 10))
        self.camera.setViewCenter(QVector3D(x_com, y_com, z_com))
        self.camera.setUpVector(QVector3D(0.0, 0.0, 1.0))
        # add light
        lightEntity = QEntity(self.rootEntity)
        light = QDirectionalLight(lightEntity)
        light.setColor(_Qt.white)
        light.setIntensity(1)
        lightEntity.addComponent(light)
        lightTransform = QTransform(lightEntity)
        lightTransform.setTranslation(self.camera.position())
        lightEntity.addComponent(lightTransform)
        # For camera controls
        camController = QOrbitCameraController(self.rootEntity)
        camController.setLinearSpeed(-20)
        camController.setLookSpeed(-90)
        camController.setCamera(self.camera)
        self.view.setRootEntity(self.rootEntity)
        # # info_box.setStandardButtons(QMessageBox.Close)
