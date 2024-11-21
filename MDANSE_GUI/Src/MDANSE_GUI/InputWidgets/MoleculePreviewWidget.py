from qtpy.QtWidgets import QWidget, QSizePolicy, QVBoxLayout, QLabel, QDialog
from qtpy.Qt3DExtras import Qt3DWindow
from qtpy.QtGui import QColor
from qtpy.Qt3DRender import QDirectionalLight, QGeometryRenderer
from qtpy.QtGui import QColor, QVector3D, QQuaternion, QFont
from qtpy.Qt3DExtras import QPhongMaterial, QCylinderMesh, \
    QCuboidMesh, QPlaneMesh, QSphereMesh, Qt3DWindow, QOrbitCameraController
from qtpy.QtCore import Qt as _Qt
from qtpy.Qt3DCore import QEntity, QTransform
from MDANSE.Chemistry import ATOMS_DATABASE
import numpy as np

class MoleculePreviewWidget(QDialog):
    def __init__(self, parent, molecule_information, molecule_name):
        super().__init__(parent)
        self.setWindowTitle("Molecule Preview")
        self.resize(800, 600)
        self.view = Qt3DWindow()
        self.view.defaultFrameGraph().setClearColor(QColor(0x4d4d4f)) #molecular viewer mdansechemistry atoms.json
        container = QWidget.createWindowContainer(self.view) #from mdanse chemistry atoms database
        screenSize = self.view.screen().size()
        container.setMinimumSize(200, 100)
        container.setMaximumSize(screenSize)
        container.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.MinimumExpanding)
        layout = QVBoxLayout()
        layout.addWidget(container)        
        self.rootEntity = QEntity()
        self.cuboidTransform = QTransform()
        self.axes = []
        mass = []
        coords = []
        info_text = f"Molecule name: {molecule_name}\n"
        for key in molecule_information["atom_number"]:
            info_text += f"Number of {key} atoms: {molecule_information['atom_number'][key]}\n"

        info_text += f"Number of such molecules in trajectory: {molecule_information['no_of_molecules']}\n"

        for i, atom in enumerate(molecule_information["atom_information"]):
            x, y, z = atom["coords"]
            symbol = atom["symbol"]            
            colour = ATOMS_DATABASE.get_atom_property(symbol, "color")
            radius = ATOMS_DATABASE.get_atom_property(symbol, "covalent_radius")
            mass.append(ATOMS_DATABASE.get_atom_property(symbol, "atomic_weight"))
            coords.append(atom["coords"])
            r, g, b = [int(x) for x in colour.split(";")]
            colour = QColor(r, g, b)
            m_sphereEntity = QEntity(self.rootEntity)
            sphereMesh = QSphereMesh()
            sphereMesh.setRings(20)
            sphereMesh.setSlices(20)
            sphereMesh.setRadius(radius*10)
            sphereTransform = QTransform()
            sphereTransform.setScale(1.0)
            sphereTransform.setTranslation(QVector3D(20*x-10, 20*y-10, 20*z-10))
            sphereMaterial = QPhongMaterial()
            sphereMaterial.setDiffuse(colour)
            sphereMaterial.setAmbient(colour)  # Set ambient to the same as diffuse.
            # sphereMaterial.setSpecular(QColor(0, 0, 0))  # Eliminate specular reflection.
            # sphereMaterial.setShininess(0.0)  # Eliminate shininess.
            m_sphereEntity.addComponent(sphereMesh)
            m_sphereEntity.addComponent(sphereMaterial)
            m_sphereEntity.addComponent(sphereTransform)

        info_label = QLabel(info_text)
        font = QFont("Arial", 12)
        info_label.setFont(font)
        info_label.setWordWrap(True)
        layout.addWidget(info_label)
        self.setLayout(layout)
        mass = np.array(mass)
        coords = np.array(coords)
        com = np.einsum("i,ik->k", mass, coords)/np.sum(mass)
        x_com, y_com, z_com = 20*com-10
        # for n, x in enumerate(['x', 'y', 'z']):
        #     tempcyl = QCylinderMesh()
        #     tempcyl.setRadius(0.05)
        #     tempcyl.setLength(100)
        #     tempcyl.setRings(50)
        #     tempcyl.setSlices(20)
        #     tempTransform = QTransform()
        #     tempTransform.setScale(1.0)
        #     tempMaterial = QPhongMaterial()
        #     if n == 0:
        #         tempTransform.setRotation(QQuaternion.fromAxisAndAngle(QVector3D(1.0, 0.0, 0.0), 90.0))
        #         tempTransform.setTranslation(QVector3D(x_com, y_com, z_com))
        #         tempMaterial.setDiffuse(QColor(0xBB0000))
        #     elif n == 1:
        #         tempTransform.setRotation(QQuaternion.fromAxisAndAngle(QVector3D(0.0, 1.0, 0.0), 90.0))
        #         tempTransform.setTranslation(QVector3D(x_com, y_com, z_com))
        #         tempMaterial.setDiffuse(QColor(0x00BB00))
        #     else:
        #         tempTransform.setRotation(QQuaternion.fromAxisAndAngle(QVector3D(0.0, 0.0, 1.0), 90.0))
        #         tempTransform.setTranslation(QVector3D(x_com, y_com, z_com))
        #         tempMaterial.setDiffuse(QColor(0x0000BB))
        #     m_cylinderEntity = QEntity(self.rootEntity)
        #     m_cylinderEntity.addComponent(tempcyl)
        #     m_cylinderEntity.addComponent(tempMaterial)
        #     m_cylinderEntity.addComponent(tempTransform)
        #     self.axes.append(m_cylinderEntity)
        # Camera
        self.camera = self.view.camera()
        self.camera.lens().setPerspectiveProjection(45.0, 16.0/9.0, 0.1, 1000.0)
        self.camera.setPosition(QVector3D(x_com+5, y_com+5, z_com+10))
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



