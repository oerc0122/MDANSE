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
from __future__ import annotations

import copy
from typing import TYPE_CHECKING

from more_itertools import nth
from qtpy.QtCore import QObject, Qt, Signal, Slot
from qtpy.QtGui import QBrush, QStandardItem, QStandardItemModel
from qtpy.QtWidgets import (
    QAbstractScrollArea,
    QComboBox,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QTableView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from MDANSE.Framework.QVectors.IQVectors import IQVectors
from MDANSE.MLogging import LOG
from MDANSE_GUI.InputWidgets.ComboWidget import ComboWidget
from MDANSE_GUI.InputWidgets.WidgetBase import Layouts, WidgetBase
from MDANSE_GUI.Tabs.Models.PlottingContext import PlottingContext
from MDANSE_GUI.Tabs.Views.PlotDataView import (
    qvector_binning_from_dict,
    vector_angular_datasets,
    vector_projection_datasets,
    vector_q_statistics_datasets,
)
from MDANSE_GUI.Tabs.Visualisers.PlotWidget import PlotWidget
from MDANSE_GUI.Utils import block_signals

if TYPE_CHECKING:
    from MDANSE.Framework.Parameters import QVectors


def numerator_suffix(index: int) -> str:
    """Return the text numerator suffix for the input number.

    For numbers ending with 1 is will return 'st' for 'first',
    except for numbers ending with 11 where it will be 'th' for 'eleventh'.
    The same reasoning follows for 2nd and 3rd.
    """
    if 10 <= index % 100 <= 20:
        return "th"
    match index % 10:
        case 1:
            return "st"
        case 2:
            return "nd"
        case 3:
            return "rd"
        case _:
            return "th"


def parameters_have_changed(qvec_configurator, plotting_object) -> bool:
    """Check if the vector generator type and input parameters have changed."""
    new_params = qvec_configurator["parameters"]
    new_type = qvec_configurator["vector_type"]
    if (
        new_type != plotting_object.last_vec_type
        or new_params != plotting_object.last_vec_params
    ):
        plotting_object.last_vec_type = new_type
        plotting_object.last_vec_params = copy.copy(new_params)
        return True
    return False


class VectorModel(QStandardItemModel):
    """Qt model for selecting vector generator types."""

    type_changed = Signal()
    unit_cell_missing = Signal()
    input_is_valid = Signal(bool)

    def __init__(self, *args, trajectory=None, **kwargs):
        super().__init__(*args, **kwargs)
        self._generator = None
        self._defaults = []
        self._trajectory = trajectory

    @Slot(str)
    def switch_qvector_type(
        self,
        vector_type: str,
        optional_settings: dict | None = None,
    ):
        """Create a vector generator of the input type and updates the input parameters.

        Parameters
        ----------
        vector_type : str
            Name of the IQVectors subclass to be created.
        optional_settings : dict | None, optional
            Dictionary of vector generator input parameters, by default None
        """
        self.clear()
        self._defaults = []
        self._generator = IQVectors.create(vector_type, self._trajectory.unit_cell(0))
        settings = self._generator.settings
        for kv in settings.items():
            name = kv[0]  # dictionary key
            if optional_settings is None:
                value = kv[1][1]["default"]  # tuple value 1: dictionary
            else:
                try:
                    value = optional_settings[name]
                except KeyError:
                    value = kv[1][1]["default"]  # tuple value 1: dictionary
            self._defaults.append(value)
            vtype = kv[1][0]  # tuple value 0: type
            items = [QStandardItem(str(x)) for x in [name, value, vtype]]
            for it in items[0::2]:
                it.setEditable(False)
            for it in items[1::2]:
                it.setData(value, role=Qt.ItemDataRole.ToolTipRole)
            self.appendRow(items)
        self.type_changed.emit()

    def params_summary(self) -> dict:
        """Validate input types and return a dictionary of input parameters."""
        params = {}
        all_inputs_are_valid = True
        for rownum in range(self.rowCount()):
            name = str(self.item(rownum, 0).text())
            value = str(self.item(rownum, 1).text())
            vtype = str(self.item(rownum, 2).text())
            try:
                params[name] = self.parse_vtype(vtype, value, name)
            except ValueError:
                params[name] = "failed"
            if params[name] == "failed":
                self.item(rownum, 1).setData(
                    QBrush(Qt.GlobalColor.red),
                    role=Qt.ItemDataRole.BackgroundRole,
                )
                all_inputs_are_valid = False
            else:
                self.item(rownum, 1).setData(0, role=Qt.ItemDataRole.BackgroundRole)
        self.input_is_valid.emit(all_inputs_are_valid)
        return params

    def parse_vtype(self, vtype: str, value: str, vname: str):
        """Validate the inputs of a type accepting multiple numbers."""
        if vtype in ("RangeConfigurator", "VectorConfigurator"):
            inner_type = self._generator.settings[vname][1]["valueType"]
            tempstring = value.strip("()[] ")
            result = [inner_type(x) for x in tempstring.split(",")]
            if len(result) == 3:
                return result
        elif vtype == "FloatConfigurator":
            return float(value)
        elif vtype == "IntegerConfigurator":
            return int(value)
        else:
            return value

        return "failed"


class ShellPanel(QWidget):
    """GUI element for plotting the vectors in a single shell.

    It stores a reference to the vector generator and uses a spinbox
    to select different shells in the generator.
    """

    def __init__(
        self,
        *args,
        qvec_configurator: QVectorsConfigurator | None = None,
        plotter_type: str = "Vectors3D",
        dialog_reference: QDialog,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.qvec_config = qvec_configurator
        self.plotter_type = plotter_type
        self.parent_dialog = dialog_reference
        self.last_vec_type = ""
        self.last_vec_params = {}
        self.current_plotted_shell = -1
        self.create_layout()
        self.set_shell(0)

    def create_layout(self):
        """Create the GUI elements of a vector shell viewer."""
        layout = QGridLayout()
        self.setLayout(layout)
        box_label = QLabel("Vector shell: ", self)
        self.shell_selector = QSpinBox(self)
        self.shell_selector.setValue(0)
        self.shell_selector.valueChanged.connect(self.set_shell)
        self.shell_information = QLabel("th shell, |q| = N/A")
        self.plot_widget = PlotWidget(self, plotter_type=self.plotter_type)
        self.plot_widget._figure.set_layout_engine("constrained")
        self.plot_widget.plot_selector.setVisible(False)
        self.plot_widget._normaliser.setVisible(False)
        self.plot_widget._sliderpack.setVisible(False)
        layout.addWidget(box_label, 0, 0)
        layout.addWidget(self.shell_selector, 0, 1)
        layout.addWidget(self.shell_information, 0, 2)
        layout.addWidget(self.plot_widget, 1, 0, 2, 3)

    def update_label(self, shell_index: int):
        """Set the GUI label text to the current shell number and q value."""
        vec_dict = self.qvec_config["q_vectors"]
        self.shell_information.setText(
            f"{numerator_suffix(shell_index)} shell, |q| = {nth(vec_dict, shell_index)}",
        )

    @Slot(int)
    def set_upper_limit(self, number_of_shells: int):
        """Set the spinbox upper limit to the last q vector shell index."""
        current_shell = min(self.shell_selector.value(), number_of_shells - 1)
        self.shell_selector.setMaximum(number_of_shells - 1)
        with block_signals(self.shell_selector):
            self.shell_selector.setValue(current_shell)
        self.set_shell(current_shell, update_limit=False)

    @Slot(int)
    def set_shell(self, shell_index: int, *, update_limit: bool = True):
        """Change the index of the current vector shell and update the plots.

        Parameters
        ----------
        shell_index : int
            Vector shell index in the vector generator.
        update_limit : bool, optional
            Update the spinbox upper limit, by default True
        """
        if (
            not parameters_have_changed(self.qvec_config, self)
            and shell_index == self.current_plotted_shell
        ):
            return
        vec_dict = self.qvec_config["q_vectors"]
        if update_limit:
            self.set_upper_limit(len(vec_dict))
        try:
            vec_key = nth(vec_dict, shell_index)
        except KeyError:
            LOG.warning(
                "Shell %i was not available in the Q vector generator",
                shell_index,
            )
            self.plot_widget.plot_blank()
            return
        else:
            if vec_key is None:
                self.set_shell(self.shell_selector.value())
                return
            self.update_label(shell_index)
        if vec_dict[vec_key] is None:
            LOG.warning("Shell %i does not contain any vectors", shell_index)
            self.plot_widget.plot_blank(
                draw_cross=False,
                override_title=f"Shell {shell_index}, q={vec_key}, contains no vectors.",
                override_label="To introduce vectors at this Q, you may have to run the simulation with a larger cell.",
            )
            return
        if not self.parent_dialog.isVisible():
            return
        new_model = PlottingContext()
        for dataset in vector_projection_datasets(self.qvec_config, vec_key):
            new_model.add_dataset(dataset)
        for dataset in vector_angular_datasets(self.qvec_config, vec_key):
            new_model.add_dataset(dataset)
        self.plot_widget.set_plotter(self.plotter_type)
        self.plot_widget.set_context(new_model)
        self.plot_widget.use_grid()
        self.plot_widget.use_legend()
        self.plot_widget.plot_data()
        self.current_plotted_shell = shell_index
        new_model.needs_an_update.connect(self.update_plot)

    @Slot()
    def update_plot(self):
        """Draw the plots for the current vector shell."""
        self.set_shell(self.shell_selector.value())


class VectorViewer(QDialog):
    """A pop-up dialog for plotting vector distribution previews.

    Attributes
    ----------
    _helper_title : str
        The title of the helper dialog window.

    """

    _helper_title = "Q Vector Viewer"

    def __init__(
        self,
        parent: QObject,
        *args,
        **kwargs,
    ):
        """Create the selection dialog.

        Parameters
        ----------
        parent : QObject
            parent object in the Qt object hierarchy
        *args : Any, ...
            catches all the arguments that may be passed to the QDialog constructor
        **kwargs : dict[str, Any]
            catches all the keyword arguments passed to the QDialog constructor

        """
        super().__init__(parent, *args, **kwargs)
        self.tab_index = {}
        self.last_vec_type = ""
        self.last_vec_params = {}
        self.setWindowTitle(self._helper_title)
        self.setWindowFlags(Qt.Window)
        layout = QVBoxLayout()
        self.main_widget = QTabWidget(self)
        layout.addWidget(self.main_widget)
        self.setLayout(layout)
        self.plot_widget = PlotWidget(self, plotter_type="Vectors")
        self.plot_widget._figure.set_layout_engine("constrained")
        self.plot_widget.plot_selector.setVisible(False)
        self.plot_widget._normaliser.setVisible(False)
        self.plot_widget._sliderpack.setVisible(False)
        self.tab_index["statistics"] = self.main_widget.addTab(
            self.plot_widget, "Vector |q| statistics"
        )
        self.shell_panel_3D = ShellPanel(
            qvec_configurator=configurator,
            plotter_type="Vectors3D",
            dialog_reference=self,
        )
        self.tab_index["shells"] = self.main_widget.addTab(
            self.shell_panel_3D, "Vector angle statistics"
        )
        self.main_widget.currentChanged.connect(self.update_plot)

    def calculate_vector_statistics(self):
        if not parameters_have_changed(self.qvec_configurator, self):
            return
        model = PlottingContext()
        modq_binning = qvector_binning_from_dict(self.qvec_configurator["parameters"])
        try:
            for qvec_dataset in vector_q_statistics_datasets(
                self.qvec_configurator,
                q_bin_limits=modq_binning,
            ):
                model.add_dataset(qvec_dataset)
        except ValueError:
            self.plot_widget.plot_blank(
                override_title="No vectors have been generated.",
                override_label="Please change the input parameters.",
            )
            return
        self.plot_widget.set_plotter("Vectors")
        self.plot_widget.set_context(model)
        self.plot_widget.use_grid()
        self.plot_widget.use_legend()
        self.plot_widget.plot_data()
        model.needs_an_update.connect(self.update_plot)

    @Slot()
    def update_plot(self):
        """Draw the vector plots using the latest data."""
        current_tab = self.main_widget.currentIndex()
        stats = self.tab_index["statistics"]
        shells = self.tab_index["shells"]
        if current_tab == stats:
            self.calculate_vector_statistics()
            self.plot_widget.plot_data(update_only=True)
        elif current_tab == shells:
            self.shell_panel_3D.update_plot()


class QVectorsParams(QObject):
    value_changed = Signal()
    _relative_size = 1

    def __init__(
        self,
        parent: QObject | None = None,
        *args,
        parameter: IQVectors,
        label: str = "Params",
        tooltip: str = "Q-Vector Params",
        layout: QLayout = QVBoxLayout,
    ):
        super().__init__(*args, parent=parent)
        from MDANSE_GUI.Tabs.Visualisers.Action import widget_lookup

        self.parameter = parameter

        self._widgets_in_layout = {}
        self._raw_widgets = {}

        self._base = QGroupBox(label, parent)
        self._base.setToolTip(tooltip)
        self._layout = layout(self._base)
        self._base.setLayout(self._layout)

        for key, desc in self.parameter.descriptors.items():
            widget_class = widget_lookup[type(desc)]
            LOG.info(f"Setting up {widget_class.__name__} for {key}")

            input_widget = widget_class(
                parent=self._base,
                base_type="QWidget",
                label=key,
                configurable=self.parameter,
                prop=key,
                parameter=desc,
            )
            self.add_widget(input_widget, key)

    def add_widget(
        self,
        input_widget: WidgetBase,
        key: str,
        *,
        add_to_layout: bool = True,
    ) -> None:
        if add_to_layout:
            self._layout.addWidget(
                input_widget._base, stretch=input_widget._relative_size
            )
        self._widgets_in_layout[key] = input_widget._base
        self._raw_widgets[key] = input_widget
        input_widget.value_changed.connect(self.value_changed.emit)


class QVectorsWidget(WidgetBase):
    def __init__(
        self,
        *args,
        parameter: QVectors,
        layout_type: Layouts = "QVBoxLayout",
        **kwargs,
    ):
        super().__init__(*args, parameter=parameter, layout_type=layout_type, **kwargs)

        self._relative_size = 3

        self.helper = None

        self._selector = ComboWidget(
            parent=self._base,
            base_type="QWidget",
            label="Generator type",
            tooltip="The q-vectors will be generated using the method chosen here.",
            choices=IQVectors.indirect_subclasses(),
            parent_widget=self,
            configurable=self.parameter,
            prop="generator",
            parameter=self.parameter.descriptors["generator"],
        )
        self._selector.value_changed.connect(self.switch_qvector_type)

        self._preview_button = QPushButton("Preview vector distribution")
        self._preview_button.clicked.connect(self.helper_dialog)

        top_bar_layout = QHBoxLayout()
        top_bar_layout.addWidget(self._selector._base, stretch=1)
        top_bar_layout.addStretch(1)
        top_bar_layout.addWidget(self._preview_button)

        self._layout.addLayout(top_bar_layout)
        self.switch_qvector_type()
        self.value_changed.connect(self.preview_vectors)

    def switch_qvector_type(self):
        if "params" in self._raw_widgets:
            self.remove_widget("params")

        params = QVectorsParams(
            parent=self._base,
            label="Generator Parameters",
            tooltip="The q-vectors will be generated using the method chosen here.",
            parameter=self.parameter.generator,
        )
        self.add_widget(params, "params")

    def get_widget_value(self):
        """Collect the results from the input widgets and return the value."""

    def create_helper(
        self,
    ) -> VectorViewer:
        """Create the vector plotting dialog.

        The PlotWidget is modified compared to the normal plotter,
        not allowing the user to change the plotting mode from 'vectors'.
        """
        dialog_instance = VectorViewer(self._base)
        self.new_shell_number.connect(dialog_instance.shell_panel_3D.set_upper_limit)
        dialog_instance.plot_widget.plot_selector.setVisible(False)
        dialog_instance.plot_widget._normaliser.setVisible(False)
        dialog_instance.plot_widget._sliderpack.setVisible(False)
        return dialog_instance

    @Slot()
    def helper_dialog(self) -> None:
        """Open the helper dialog."""
        if self.helper is None:
            self.helper = self.create_helper()

        if self.helper.isVisible():
            geometry = self.helper.saveGeometry()
            self.helper.previous_geometry = geometry
            self.helper.close()
        else:
            if hasattr(self.helper, "previous_geometry"):
                self.helper.restoreGeometry(self.helper.previous_geometry)
            self.helper.show()
            self.preview_vectors()

    @Slot()
    def preview_vectors(self):
        """Build the data sets and plot them in the dialog window.

        If the dialog is not open, this function exits immediately.
        """
        if self.helper is None:
            self.helper = self.create_helper()

        if not self.helper.isVisible():
            return

        if not self.parameter.check_status():
            self.helper.plot_widget._plotter.plot_blank()
            return
        self.helper.update_plot()

    def updateValue(self):
        temp = super().updateValue()
        self._preview_button.setEnabled(self.parameter.check_status())
        return temp
