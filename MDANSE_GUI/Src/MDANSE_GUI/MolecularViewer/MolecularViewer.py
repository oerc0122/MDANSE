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
from typing import List, Tuple, Dict, Any

import numpy as np
from scipy.spatial import cKDTree as KDTree

from qtpy import QtWidgets
from qtpy.QtCore import Signal, Slot
from qtpy.QtWidgets import (
    QSizePolicy,
    QDialog,
    QFormLayout,
    QPushButton,
    QLineEdit,
    QDoubleSpinBox,
    QSpinBox,
    QLabel,
)

import vtk
from vtk.util import numpy_support
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.vtkRenderingAnnotation import vtkAxesActor

from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from MDANSE.MLogging import LOG

from MDANSE_GUI.MolecularViewer.readers import hdf5wrapper
from MDANSE_GUI.MolecularViewer.Dummy import PyConnectivity
from MDANSE_GUI.MolecularViewer.Contents import TrajectoryAtomData
from MDANSE_GUI.MolecularViewer.AtomProperties import (
    AtomProperties,
    ndarray_to_vtkarray,
)


def array_to_3d_imagedata(data: np.ndarray, spacing: Tuple[float]):
    nx = data.shape[0]
    ny = data.shape[1]
    nz = data.shape[2]
    image = vtk.vtkImageData()
    image.SetDimensions(nx, ny, nz)
    dx, dy, dz = spacing
    image.SetSpacing(dx, dy, dz)
    image.SetExtent(0, nx - 1, 0, ny - 1, 0, nz - 1)
    if vtk.vtkVersion.GetVTKMajorVersion() < 6:
        image.SetScalarTypeToDouble()
        image.SetNumberOfScalarComponents(1)
    else:
        image.AllocateScalars(vtk.VTK_DOUBLE, 1)

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                image.SetScalarComponentFromDouble(i, j, k, 0, data[i, j, k])

    return image


def smear_grid(grid: np.ndarray, fine_sampling: int) -> np.ndarray:
    """Include atom radius effect in the array of atom counts on a grid.

    Parameters
    ----------
    grid : np.ndarray
        a 3D histogram of atom positions over time
    fine_sampling : int
        the fraction of the atom radius defining the binning grid

    Returns
    -------
    np.ndarray
        a 3D histogram of the volume taken by the atom
    """
    if fine_sampling < 2:
        return grid
    final_histogram = grid.copy()
    for _ in range(1, fine_sampling):
        n = 1
        new_histogram = np.zeros_like(grid)
        new_histogram[:, :, n:] += final_histogram[:, :, :-n]
        new_histogram[:, :, :-n] += final_histogram[:, :, n:]
        new_histogram[:, n:, :] += final_histogram[:, :-n, :]
        new_histogram[:, :-n, :] += final_histogram[:, n:, :]
        new_histogram[n:, :, :] += final_histogram[:-n, :, :]
        new_histogram[:-n, :, :] += final_histogram[n:, :, :]
        final_histogram += new_histogram
    return final_histogram


class MolecularViewer(QtWidgets.QWidget):
    """This class implements a molecular viewer."""

    new_max_frames = Signal(int)
    changed_trace = Signal()

    def __init__(self):
        super(MolecularViewer, self).__init__()

        self._scale_factor = 0.4

        self._datamodel = None
        self._element_database = None

        self._iren = QVTKRenderWindowInteractor(self)

        def dummy_method(self, ev=None):
            pass

        setattr(self._iren, "keyPressEvent", dummy_method)

        self._renderer = vtk.vtkRenderer()

        self._iren.GetRenderWindow().AddRenderer(self._renderer)

        self._iren.GetRenderWindow().SetPosition((0, 0))

        self._iren.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

        self._iren.Enable()

        self._iren.GetRenderWindow()

        self.axes_actor = vtkAxesActor()
        self.axes_actor.AxisLabelsOn()
        self.axes_actor.SetShaftTypeToCylinder()
        self.axes_widget = vtk.vtkOrientationMarkerWidget()
        self.axes_widget.SetOrientationMarker(self.axes_actor)
        self.axes_widget.SetInteractor(self._iren.GetRenderWindow().GetInteractor())
        self.axes_widget.SetViewport(0.0, 0.0, 0.25, 0.25)
        self.axes_widget.SetEnabled(True)
        self.axes_widget.InteractiveOff()

        layout = QtWidgets.QStackedLayout(self)
        layout.addWidget(self._iren)
        self.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Minimum)

        # create camera
        self._camera = vtk.vtkCamera()
        # associate camera to renderer
        self._renderer.SetActiveCamera(self._camera)
        self._camera.SetFocalPoint(0, 0, 0)
        self._camera.SetPosition(0, 0, 20)

        self._n_atoms = 0
        self._n_frames = 0
        self._resolution = 0

        self._atoms_visible = True
        self._bonds_visible = True
        self._axes_visible = True
        self._cell_visible = True

        self._iren.Initialize()

        self._atoms = []

        self._polydata = None
        self._polydata_bonds_exist = False
        self._uc_polydata = None

        self._surfaces = []
        self._isocontours = []

        self._reader = None

        self._current_frame = 0

        self._colour_manager = AtomProperties()
        self.dummy_size = 0.0

        self.reset_camera = False
        self.create_trace_dialog()

    def setDataModel(self, datamodel: TrajectoryAtomData):
        self._datamodel = datamodel

    def _new_trajectory_object(self, fname: str, data: HDFTrajectoryInputData):
        reader = hdf5wrapper.HDF5Wrapper(fname, data.trajectory, data.chemical_system)
        self.set_reader(reader)

    @Slot(str)
    def _new_trajectory(self, fname: str):
        data = HDFTrajectoryInputData(fname)
        reader = hdf5wrapper.HDF5Wrapper(fname, data.trajectory, data.chemical_system)
        self.set_reader(reader)

    @Slot(float)
    def _new_scaling(self, scale_factor: float):
        self._scale_factor = scale_factor
        self.update_renderer()

    def _new_visibility(self, flags: List[bool]):
        self._atoms_visible = flags[0]
        self._bonds_visible = flags[1]
        self._axes_visible = flags[2]
        self._cell_visible = flags[3]
        self.axes_widget.SetEnabled(self._axes_visible)
        result = self.set_coordinates(self._current_frame)
        if result is False:
            self.update_renderer()

    @Slot()
    def _new_atom_parameters(self):
        if self._polydata is None or self._datamodel is None:
            return

        # we need to add the new colours to LUT
        # then assign new colour NUMBERS to _atom_colours

        # we also need to assign new atom radii to _atom_scales

        scalars = ndarray_to_vtkarray(
            self._atom_colours, self._atom_scales, self._n_atoms
        )

        self._polydata.GetPointData().SetScalars(scalars)

        self.update_renderer()

    def trace_from_dialog(self, params: Dict[str, Any]):
        self._draw_isosurface(params["atom_number"], params)

    def delete_from_dialog(self, trace_number: int):
        try:
            surface = self._surfaces[trace_number]
        except IndexError:
            return
        else:
            surface.VisibilityOff()
            surface.ReleaseGraphicsResources(self._iren.GetRenderWindow())
            self._renderer.RemoveActor(surface)
            self._surfaces.pop(trace_number)
            self._iren.Render()
            self.changed_trace.emit()

    def _draw_isosurface(self, index: int, params=None):
        """Draw the isosurface of an atom with the given index"""

        if self._reader is None:
            return

        LOG.info("Computing isosurface ...")
        if params is not None:
            fine_sampling = params.get("fine_sampling", 5)
            r, g, b = params.get("surface_colour", (0, 0.5, 0.75))
            opacity = params.get("surface_opacity", 0.5)
            trace_cutoff = params.get("trace_cutoff", 90)
        else:
            fine_sampling = 5
            r, g, b = 0, 0.5, 0.75
            opacity = 0.5
            trace_cutoff = 90

        coords = self._reader.read_atom_trajectory(index)
        element = self._reader._atom_types[index]
        radius = self._reader._trajectory.get_atom_property(element, "covalent_radius")
        upper_limit = np.max(coords, axis=0) + radius
        lower_limit = np.min(coords, axis=0) - radius
        grid_step = radius / fine_sampling

        span = upper_limit - lower_limit
        grid_steps = [int(x) for x in span / grid_step]
        gdim = (grid_steps[0], grid_steps[1], grid_steps[2])
        grid = np.zeros(gdim, dtype=np.int32)

        indices = np.floor((coords - lower_limit.reshape((1, 3))) / grid_step).astype(
            int
        )
        unique_indices, counts = np.unique(indices, return_counts=True, axis=0)
        grid[tuple(unique_indices.T)] += counts
        self._atomic_trace_histogram = smear_grid(grid, fine_sampling)

        self._image = array_to_3d_imagedata(
            self._atomic_trace_histogram, (grid_step, grid_step, grid_step)
        )
        isovalue = np.percentile(self._atomic_trace_histogram, trace_cutoff)

        new_isocontour = vtk.vtkMarchingContourFilter()
        new_isocontour.UseScalarTreeOn()
        new_isocontour.ComputeNormalsOn()
        if vtk.vtkVersion.GetVTKMajorVersion() < 6:
            new_isocontour.SetInput(self.image)
        else:
            new_isocontour.SetInputData(self._image)
        new_isocontour.SetValue(0, isovalue)

        self._depthSort = vtk.vtkDepthSortPolyData()
        self._depthSort.SetInputConnection(new_isocontour.GetOutputPort())
        self._depthSort.SetDirectionToBackToFront()
        self._depthSort.SetVector(1, 1, 1)
        self._depthSort.SetCamera(self._camera)
        self._depthSort.SortScalarsOn()
        self._depthSort.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(self._depthSort.GetOutputPort())
        mapper.ScalarVisibilityOff()
        mapper.Update()

        new_surface = vtk.vtkActor()
        new_surface.SetMapper(mapper)
        new_surface.GetProperty().SetColor((r, g, b))
        new_surface.GetProperty().SetOpacity(opacity)
        new_surface.PickableOff()

        new_surface.GetProperty().SetRepresentationToSurface()

        new_surface.GetProperty().SetInterpolationToGouraud()
        new_surface.GetProperty().SetSpecular(0.4)
        new_surface.GetProperty().SetSpecularPower(10)

        self._renderer.AddActor(new_surface)

        new_surface.SetPosition(lower_limit[0], lower_limit[1], lower_limit[2])
        self._surfaces.append(new_surface)
        self._isocontours.append(new_isocontour)

        self._iren.Render()

        LOG.info("... done")
        self.changed_trace.emit()

    def create_all_actors(self):
        actors = []
        if self._polydata is None:
            return actors

        line_actor, ball_actor = self.create_traj_actors(self._polydata)
        if self._cell_visible:
            uc_actor = self.create_uc_actor()
            actors.append(uc_actor)
        if self._bonds_visible and self._polydata_bonds_exist:
            actors.append(line_actor)
        if self._atoms_visible:
            actors.append(ball_actor)
        return actors

    def create_uc_actor(self):
        uc_mapper = vtk.vtkPolyDataMapper()
        if vtk.vtkVersion.GetVTKMajorVersion() < 6:
            uc_mapper.SetInput(self._uc_polydata)
        else:
            uc_mapper.SetInputData(self._uc_polydata)
        uc_mapper.ScalarVisibilityOn()
        uc_actor = vtk.vtkLODActor()
        uc_actor.GetProperty().SetLineWidth(3 * self._scale_factor)
        uc_actor.SetMapper(uc_mapper)
        return uc_actor

    def create_traj_actors(self, polydata, line_opacity=1.0, ball_opacity=1.0):
        line_mapper = vtk.vtkPolyDataMapper()
        if vtk.vtkVersion.GetVTKMajorVersion() < 6:
            line_mapper.SetInput(polydata)
        else:
            line_mapper.SetInputData(polydata)

        line_mapper.SetLookupTable(self._colour_manager._lut)
        line_mapper.ScalarVisibilityOn()
        line_mapper.ColorByArrayComponent("scalars", 1)
        line_actor = vtk.vtkLODActor()
        line_actor.GetProperty().SetLineWidth(3 * self._scale_factor)
        line_actor.SetMapper(line_mapper)
        line_actor.GetProperty().SetAmbient(0.2)
        line_actor.GetProperty().SetDiffuse(0.5)
        line_actor.GetProperty().SetSpecular(0.3)
        line_actor.GetProperty().SetOpacity(line_opacity)

        temp_radius = float(1.0 * self._scale_factor)
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(0, 0, 0)
        sphere.SetRadius(temp_radius)
        sphere.SetThetaResolution(self._resolution)
        sphere.SetPhiResolution(self._resolution)
        glyph = vtk.vtkGlyph3D()
        if vtk.vtkVersion.GetVTKMajorVersion() < 6:
            glyph.SetInput(polydata)
        else:
            glyph.SetInputData(polydata)

        temp_scale = float(1.0 * self._scale_factor)
        glyph.SetScaleModeToScaleByScalar()
        glyph.SetColorModeToColorByScalar()
        glyph.SetScaleFactor(temp_scale)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.SetIndexModeToScalar()
        sphere_mapper = vtk.vtkPolyDataMapper()
        sphere_mapper.SetLookupTable(self._colour_manager._lut)
        sphere_mapper.SetScalarRange(polydata.GetScalarRange())
        sphere_mapper.SetInputConnection(glyph.GetOutputPort())
        sphere_mapper.ScalarVisibilityOn()
        sphere_mapper.ColorByArrayComponent("scalars", 1)
        ball_actor = vtk.vtkLODActor()
        ball_actor.SetMapper(sphere_mapper)
        ball_actor.GetProperty().SetAmbient(0.2)
        ball_actor.GetProperty().SetDiffuse(0.5)
        ball_actor.GetProperty().SetSpecular(0.3)
        ball_actor.GetProperty().SetOpacity(ball_opacity)
        ball_actor.SetNumberOfCloudPoints(30000)
        return [line_actor, ball_actor]

    def clear_trajectory(self):
        """Clear the vtk scene from atoms and bonds actors."""

        if not hasattr(self, "_actors"):
            return
        if self._actors is None:
            return

        self.on_clear_atomic_trace()
        self._actors.VisibilityOff()
        self._actors.ReleaseGraphicsResources(self.get_render_window())
        self._renderer.RemoveActor(self._actors)

        del self._actors

    def clear_panel(self) -> None:
        """Clears the Molecular Viewer panel"""
        self.clear_trajectory()

        self._reader = None

        # set everything to some empty/zero value
        self._n_atoms = 0
        self._n_frames = 0
        self.new_max_frames.emit(0)
        self._atoms = []
        self._atom_colours = []
        self._current_frame = 0
        self.reset_all_polydata()

        self.update_renderer()

        # clear the atom properties table
        self._colour_manager.removeRows(0, self._colour_manager.rowCount())

    def reset_all_polydata(self):
        self._polydata = vtk.vtkPolyData()
        self._uc_polydata = vtk.vtkPolyData()

    def update_all_polydata(self):
        self.update_polydata()
        self.update_uc_polydata()

    def update_polydata(self):
        coords = self._reader.read_frame(self._current_frame)

        if self._atoms_visible or self._bonds_visible:
            atoms = vtk.vtkPoints()
            atoms.SetData(numpy_support.numpy_to_vtk(coords))
            self._polydata.SetPoints(atoms)

        if self._bonds_visible:
            # do not bond atoms to dummy atoms
            rs = coords[self.not_du]
            bonds, bonds_exist = self.create_bond_cell_array(
                rs, self.covs[self.not_du], self.not_du
            )
            if bonds_exist:
                self._polydata.SetLines(bonds)
                self._polydata_bonds_exist = True
                return

        self._polydata_bonds_exist = False

    def create_bond_cell_array(self, rs, covs, not_du, tolerance=0.04):
        # determine and set bonds without PBC applied
        bonds = vtk.vtkCellArray()

        tree = KDTree(rs)
        contacts = tree.query_ball_point(rs, 2 * np.max(covs) + tolerance, workers=-1)
        n_dists = sum([len(i) for i in contacts])

        js = np.zeros(n_dists, dtype=int)
        ks = np.zeros(n_dists, dtype=int)
        start = 0
        for i, idxs in enumerate(contacts):
            n_idxs = len(idxs)
            if n_idxs == 0:
                continue
            js[start : start + n_idxs] = i
            ks[start : start + n_idxs] = idxs
            start += n_idxs

        mask = js < ks
        js = js[mask]
        ks = ks[mask]
        diff = rs[js] - rs[ks]
        dist = np.sum(diff * diff, axis=1)
        sum_radii = (covs[js] + covs[ks] + tolerance) ** 2
        js = js[(0 < dist) & (dist < sum_radii)]
        ks = ks[(0 < dist) & (dist < sum_radii)]
        ls = not_du[js]
        ms = not_du[ks]

        n_points = len(ls)
        idxs = np.zeros((len(ls), 3), dtype=np.int64)
        idxs[:, 0] = 2
        idxs[:, 1] = ls
        idxs[:, 2] = ms
        bonds.SetCells(n_points, numpy_support.numpy_to_vtkIdTypeArray(idxs.flatten()))

        return bonds, len(ls) > 0

    def update_uc_polydata(self):
        uc = self._reader.read_pbc(self._current_frame)
        if self._cell_visible and uc is not None:
            # update the unit cell
            a = uc.a_vector
            b = uc.b_vector
            c = uc.c_vector
            uc_points = vtk.vtkPoints()
            uc_points.SetNumberOfPoints(8)
            for i, v in enumerate([[0, 0, 0], a, b, c, a + b, a + c, b + c, a + b + c]):
                x, y, z = v
                uc_points.SetPoint(i, x, y, z)
            self._uc_polydata.SetPoints(uc_points)

            uc_lines = vtk.vtkCellArray()
            for i, j in [
                (0, 1),
                (0, 2),
                (0, 3),
                (1, 4),
                (1, 5),
                (4, 7),
                (2, 4),
                (2, 6),
                (5, 7),
                (3, 5),
                (3, 6),
                (6, 7),
            ]:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, i)
                line.GetPointIds().SetId(1, j)
                uc_lines.InsertNextCell(line)
            self._uc_polydata.SetLines(uc_lines)

    def get_atom_index(self, pid):
        """Return the atom index from the vtk data point index.

        Args:
            pid (int): the data point index
        """

        _, _, idx = (
            self.glyph.GetOutput().GetPointData().GetArray("scalars").GetTuple3(pid)
        )

        return int(idx)

    def get_render_window(self):
        """Returns the render window."""
        return self._iren.GetRenderWindow()

    @property
    def iren(self):
        return self._iren

    def on_change_atomic_trace_opacity(self, surface_number, opacity):
        """Event handler called when the opacity level is changed."""

        if len(self._surfaces) <= surface_number:
            return

        self._surfaces[surface_number].GetProperty().SetOpacity(opacity)
        self._iren.Render()

    def on_change_atomic_trace_isocontour_level(self, surface_number, level):
        """Event handler called when the user change the isocontour level."""

        if len(self._surfaces) <= surface_number:
            return

        self._isocontours[surface_number].SetValue(0, level)
        self._isocontours[surface_number].Update()
        self._iren.Render()

    def on_change_atomic_trace_rendering_type(self, surface_number, rendering_type):
        """Event handler called when the user change the rendering type for the atomic trace."""

        if len(self._surfaces) <= surface_number:
            return

        surface = self._surfaces[surface_number]

        if rendering_type == "wireframe":
            surface.GetProperty().SetRepresentationToWireframe()
        elif rendering_type == "surface":
            surface.GetProperty().SetRepresentationToSurface()
        elif rendering_type == "points":
            surface.GetProperty().SetRepresentationToPoints()
            surface.GetProperty().SetPointSize(3)
        else:
            return

        self._iren.Render()

    def on_clear_atomic_trace(self):
        """Event handler called when the user select the 'Atomic trace -> Clear' main menu item"""

        if not self._surfaces:
            return

        for surface in self._surfaces:
            surface.VisibilityOff()
            surface.ReleaseGraphicsResources(self._iren.GetRenderWindow())
            self._renderer.RemoveActor(surface)
        self._iren.Render()

        self._surfaces = []
        self.changed_trace.emit()

    def show_trace_settings_dialog(self):
        dialog = self._trace_dialog
        if dialog.isVisible():
            if dialog.isMaximized():
                dialog.showMaximized()
            else:
                dialog.showNormal()
            dialog.activateWindow()
        else:
            dialog.show()

    def create_trace_dialog(self):
        self._trace_dialog = AtomTraceDialog(self)
        self._trace_dialog.new_atom_trace.connect(self.trace_from_dialog)
        self._trace_dialog.remove_atom_trace.connect(self.delete_from_dialog)
        self.changed_trace.connect(self._trace_dialog.update_limits)

    def on_show_atomic_trace(self):
        if self._previously_picked_atom is None:
            LOG.warning("No atom selected for computing atomic trace")
            return

        self._draw_isosurface(self._previously_picked_atom[0])

    @property
    def renderer(self):
        return self._renderer

    def set_connectivity_builder(self, coords, covalent_radii):
        # Compute the bounding box of the system
        lower_bound = coords.min(axis=0)
        upper_bound = coords.max(axis=0)

        # Enlarge it a bit to not miss any atom
        lower_bound -= 1.0e-6
        upper_bound += 1.0e-6

        # Initializes the octree used to build the connectivity
        self._connectivity_builder = PyConnectivity(lower_bound, upper_bound, 0, 10, 18)

        # Add the points to the octree
        for index, xyz, radius in zip(range(self._n_atoms), coords, covalent_radii):
            self._connectivity_builder.add_point(index, xyz, radius)

    @Slot(int)
    def set_coordinates(self, frame: int):
        """Sets a new configuration.

        @param frame: the configuration number
        @type frame: integer
        """
        if self._reader is None:
            return False

        self._current_frame = frame % self._reader.n_frames

        # update the atoms
        self.update_all_polydata()

        # Update the view.
        self.update_renderer()

    def set_reader(self, reader):
        """Set the trajectory at a given frame

        Args:
            reader (IReader): the trajectory object
        """

        if (self._reader is not None) and (reader.filename == self._reader.filename):
            return

        self.reset_camera = True
        self.clear_trajectory()

        self._reader = reader

        self._element_database = self._reader._trajectory
        self._n_atoms = self._reader.n_atoms
        self._n_frames = self._reader.n_frames
        self._current_frame = min(self._current_frame, self._n_frames - 1)

        self._atoms = self._reader.atom_types

        # Hack for reducing objects resolution when the system is big
        self._resolution = int(np.sqrt(3000000.0 / self._n_atoms))
        self._resolution = 10 if self._resolution > 10 else self._resolution
        self._resolution = 4 if self._resolution < 4 else self._resolution

        self._atom_colours = self._colour_manager.reinitialise_from_database(
            self._atoms, self._element_database, self.dummy_size
        )
        # this returns a list of indices, mapping colours to atoms

        self._atom_scales = np.array(
            [
                self._element_database.get_atom_property(at, "vdw_radius")
                for at in self._atoms
            ]
        ).astype(np.float32)
        self.du_log = np.array(
            [
                self._element_database.get_atom_property(at, "dummy") == 0
                for at in self._reader.atom_types
            ]
        )
        self.not_du = np.array(
            [
                i
                for i, at in enumerate(self._reader.atom_types)
                if self._element_database.get_atom_property(at, "dummy") == 0
            ]
        )
        self.covs = np.array(
            [
                self._element_database.get_atom_property(at, "covalent_radius")
                for at in self._reader.atom_types
            ]
        )

        scalars = ndarray_to_vtkarray(
            self._atom_colours, self._atom_scales, self._n_atoms
        )

        self.reset_all_polydata()
        self._polydata.GetPointData().SetScalars(scalars)

        self._colour_manager.onNewValues()
        self.new_max_frames.emit(self._n_frames - 1)
        self._trace_dialog.update_limits()

    @Slot(object)
    def take_atom_properties(self, data):
        colours, radii, numbers = data
        scalars = ndarray_to_vtkarray(colours, radii, numbers)
        self._polydata = vtk.vtkPolyData()
        self._polydata.GetPointData().SetScalars(scalars)
        self.update_all_polydata()
        self.update_renderer()
        # self._datamodel.setReader(reader)

    def update_renderer(self):
        """
        Update the renderer
        """
        # deleting old frame
        self.clear_trajectory()

        # creating new polydata
        self._actors = vtk.vtkAssembly()
        for actor in self.create_all_actors():
            self._actors.AddPart(actor)

        # adding polydata to renderer
        self._renderer.AddActor(self._actors)

        # rendering
        if self.reset_camera:
            self._renderer.ResetCamera()
            self.reset_camera = False

        self._iren.Render()


class MolecularViewerWithPicking(MolecularViewer):
    """This class implements a molecular viewer with picking."""

    picked_atoms_changed = Signal(object)

    def __init__(self):
        super().__init__()
        # we set dummy size to something non-zero since we need to be
        # able to see it for picking purposes
        self.dummy_size = 0.1
        self._picking_domain = None
        self._picked_polydata = None
        self._picked_polydata_bonds_exist = False
        self._polydata_opacity = 0.15
        self.picked_atoms = set()
        self.build_events()

    def build_events(self):
        """Build the events."""
        self._iren.AddObserver("LeftButtonPressEvent", self.on_pick)

    def on_pick(self, obj, event=None):
        """Event handler when an atom is mouse-picked with the left mouse button"""

        if not self._reader:
            return

        if self._picking_domain is None:
            return

        picker = vtk.vtkPropPicker()

        picker.AddPickList(self._picking_domain)
        picker.PickFromListOn()

        pos = obj.GetEventPosition()
        picker.Pick(pos[0], pos[1], 0, self._renderer)

        picked_actor = picker.GetActor()
        if picked_actor is None:
            return

        picked_pos = np.array(picker.GetPickPosition())
        coords = self._reader.read_frame(self._current_frame)
        _, idx = KDTree(coords).query(picked_pos)
        self.on_pick_atom(idx)

    def on_pick_atom(self, picked_atom):
        """Change the color of a selected atom"""
        if self._reader is None:
            return

        if picked_atom < 0 or picked_atom >= self._n_atoms:
            return

        if picked_atom in self.picked_atoms:
            self.picked_atoms.remove(picked_atom)
        else:
            self.picked_atoms.add(picked_atom)

        self.update_picked_polydata()
        self.update_renderer()
        self.picked_atoms_changed.emit(self.picked_atoms)

    def change_picked(self, picked: set[int]):
        self.picked_atoms = picked
        self.update_picked_polydata()
        self.update_renderer()

    def update_picked_polydata(self):
        atoms = vtk.vtkPoints()

        if len(self.picked_atoms) == 0:
            self._picked_polydata.SetPoints(atoms)
            return

        picked = np.array(sorted(list(self.picked_atoms)))
        coords = self._reader.read_frame(self._current_frame)
        atoms.SetData(numpy_support.numpy_to_vtk(coords[picked]))
        self._picked_polydata.SetPoints(atoms)

        scalars = ndarray_to_vtkarray(
            self._colour_manager.colours[picked],
            self._colour_manager.radii[picked],
            np.arange(len(self.picked_atoms)),
        )
        self._picked_polydata.GetPointData().SetScalars(scalars)

        not_du = np.arange(len(self.picked_atoms))[self.du_log[picked]]
        if self._bonds_visible and len(not_du) >= 1:
            # do not bond atoms to dummy atoms
            rs = coords[picked][not_du]
            covs = self.covs[picked][not_du]

            bonds, bonds_exist = self.create_bond_cell_array(rs, covs, not_du)
            if bonds_exist:
                self._picked_polydata.SetLines(bonds)
                self._picked_polydata_bonds_exist = True
                return

        self._picked_polydata_bonds_exist = False

    def reset_all_polydata(self):
        super().reset_all_polydata()
        self._picked_polydata = vtk.vtkPolyData()

    def update_all_polydata(self):
        super().update_all_polydata()
        self.update_picked_polydata()

    def create_all_actors(self):
        actors = []
        if self._polydata is None or self._picked_polydata is None:
            return actors

        line_actor, ball_actor = self.create_traj_actors(
            self._polydata,
            line_opacity=self._polydata_opacity,
            ball_opacity=self._polydata_opacity,
        )
        picked_line_actor, picked_ball_actor = self.create_traj_actors(
            self._picked_polydata
        )
        if self._cell_visible:
            uc_actor = self.create_uc_actor()
            actors.append(uc_actor)
        if self._bonds_visible and self._polydata_bonds_exist:
            actors.append(line_actor)
        if self._bonds_visible and self._picked_polydata_bonds_exist:
            actors.append(picked_line_actor)
        if self._atoms_visible:
            self._picking_domain = ball_actor
            actors.append(ball_actor)
            actors.append(picked_ball_actor)
        else:
            self._picking_domain = None
        return actors

    @Slot(object)
    def take_atom_properties(self, data):
        colours, radii, numbers = data
        scalars = ndarray_to_vtkarray(colours, radii, numbers)
        self._polydata = vtk.vtkPolyData()
        self._polydata.GetPointData().SetScalars(scalars)

        picked_colours = []
        picked_radii = []
        picked_numbers = []
        for i in sorted(list(self.picked_atoms)):
            picked_colours.append(colours[i])
            picked_radii.append(radii[i])
            picked_numbers.append(numbers[i])

        scalars = ndarray_to_vtkarray(
            np.array(picked_colours),
            np.array(picked_radii),
            np.arange(len(self.picked_atoms)),
        )
        self._picked_polydata = vtk.vtkPolyData()
        self._picked_polydata.GetPointData().SetScalars(scalars)

        self.set_coordinates(self._current_frame)


class AtomTraceDialog(QDialog):

    new_atom_trace = Signal(dict)
    remove_atom_trace = Signal(int)

    def __init__(self, parent):
        super().__init__(parent)

        self.setWindowTitle("Add atom trace to the view")
        layout = QFormLayout(self)
        self.setLayout(layout)

        self._molviewer = parent
        self.initialise_values(parent)
        self.populate_layout()

    def initialise_values(self, viewer: MolecularViewer):
        self._n_atoms = viewer._n_atoms
        self._opacity = 0.5
        self._color = (0, 0.5, 0.75)
        self._iso_percentile = 90

    def populate_layout(self):
        layout = self.layout()
        self.add_trace_button = QPushButton("Calculate and add atom trace", self)
        self.remove_trace_button = QPushButton("Remove atom trace", self)
        self.add_trace_button.clicked.connect(self.new_trace_details)
        self.remove_trace_button.clicked.connect(self.remove_trace)
        self._atom_spinbox = QSpinBox(self)
        self._surface_spinbox = QSpinBox(self)
        self._fraction_spinbox = QSpinBox(self)
        self._grid_spinbox = QSpinBox(self)
        self._opacity_spinbox = QDoubleSpinBox(self)
        self._colour_lineedit = QLineEdit("0,128,192", self)
        for sbox in [
            self._atom_spinbox,
            self._surface_spinbox,
            self._fraction_spinbox,
            self._opacity_spinbox,
        ]:
            sbox.setMinimum(0)
            sbox.setValue(0)
        self.update_limits()
        self._fraction_spinbox.setMaximum(100)
        self._fraction_spinbox.setValue(5)
        self._grid_spinbox.setMaximum(10)
        self._grid_spinbox.setMinimum(1)
        self._grid_spinbox.setValue(3)
        self._opacity_spinbox.setMaximum(1.0)
        self._opacity_spinbox.setValue(0.5)
        self._opacity_spinbox.setSingleStep(0.01)
        layout.addRow("Selected atom index: ", self._atom_spinbox)
        layout.addRow("Sampling step (1=coarse, 10=fine)", self._grid_spinbox)
        layout.addRow("Trace percentile for isovalue", self._fraction_spinbox)
        layout.addRow("Isosurface opacity", self._opacity_spinbox)
        layout.addRow("Isosurface colour (R,G,B)", self._colour_lineedit)
        layout.addWidget(self.add_trace_button)
        layout.addRow("Remove the surface with index: ", self._surface_spinbox)
        layout.addWidget(self.remove_trace_button)

    @Slot()
    def update_limits(self):
        self._atom_spinbox.setMaximum(max(self._molviewer._n_atoms - 1, 0))
        self._surface_spinbox.setMaximum(max(len(self._molviewer._surfaces) - 1, 0))
        self.enable_buttons()

    def enable_buttons(self):
        if len(self._molviewer._surfaces) == 0:
            self.remove_trace_button.setEnabled(False)
        else:
            self.remove_trace_button.setEnabled(True)
        if self._molviewer._n_atoms < 1:
            self.add_trace_button.setEnabled(False)
        else:
            self.add_trace_button.setEnabled(True)

    def get_values(self):
        params = {
            "atom_number": 0,
            "fine_sampling": 3,
            "surface_colour": (0, 0.5, 0.75),
            "surface_opacity": 0.5,
            "trace_cutoff": 5,
            "surface_number": -1,
        }
        params["atom_number"] = self._atom_spinbox.value()
        params["surface_number"] = self._surface_spinbox.value()
        try:
            params["surface_colour"] = [
                float(x) / 256 for x in self._colour_lineedit.text().split(",")
            ]
        except ValueError:
            pass
        except TypeError:
            pass
        params["trace_cutoff"] = self._fraction_spinbox.value()
        params["fine_sampling"] = self._grid_spinbox.value()
        params["surface_opacity"] = self._opacity_spinbox.value()
        return params

    def new_trace_details(self):
        params = self.get_values()
        self.new_atom_trace.emit(params)

    def remove_trace(self):
        params = self.get_values()
        self.remove_atom_trace.emit(params["surface_number"])
