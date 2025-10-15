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

from abc import abstractmethod
from enum import EnumMeta
from typing import TYPE_CHECKING, Any, Literal

from qtpy.QtCore import QObject, Signal, Slot
from qtpy.QtWidgets import (
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QVBoxLayout,
    QWidget,
)

from MDANSE.Framework.Parameters.Parameters import ConfigError, CustomChoices
from MDANSE.MLogging import LOG

if TYPE_CHECKING:
    from MDANSE.Framework.Parameters.Parameters import Parameter

Layouts = Literal["QHBoxLayout", "QVBoxLayout", "QGridLayout"]
Bases = Literal["QWidget", "QGroupBox"]

WARNING_STYLE = (
    "QWidget#InputWidget { background-color:rgb(220,210,30); font-weight: bold }"
)
ERROR_STYLE = (
    "QWidget#InputWidget { background-color:rgb(180,20,180); font-weight: bold }"
)


class WidgetBase(QObject):
    """Object to serve as an ABC to GUI widgets in MDANSE.

    Parameters
    ----------
    parent : Optional[QObject]
        Parent of widget.
    parameter : Parameter
        Parameter controlled by this widget.
    label : str
        Label of widget.
    tooltip : str
        Hover-over tooltip.
    base_type : QWidget | QGroupBox
        Base type to use.
    layout_type : QHBoxLayout | QVBoxLayout | QGridLayout
        Layout to add.
    """

    valid_changed = Signal()
    value_updated = Signal()
    value_changed = Signal()

    def __init__(
        self,
        parent: QObject | None = None,
        *args,
        parameter: Parameter,
        configurable: object,
        prop: str,
        label: str = "",
        tooltip: str = "",
        base_type: Bases = "QGroupBox",
        layout_type: Layouts = "QHBoxLayout",
        **kwargs,
    ):
        super().__init__(*args, parent=parent)
        self._parent = parent
        self._value = None
        self._relative_size = 1
        self._label_text = label
        self._tooltip = tooltip
        self._base_type = base_type
        self._layout_type = layout_type

        self.parameter = parameter
        self._configurable = configurable
        self._property = prop

        for dep in self.get_widget_deps().values():
            dep.value_changed.connect(self.updateValue)

        match self._layout_type:
            case "QHBoxLayout":
                layoutClass = QHBoxLayout
            case "QVBoxLayout":
                layoutClass = QVBoxLayout
            case "QGridLayout":
                layoutClass = QGridLayout
            case _:
                raise NotImplementedError(
                    f"Cannot create layout of type {self._layout_type}."
                )

        match self._base_type:
            case "QWidget":
                base = QWidget(parent)
                layout = layoutClass(base)
                base.setLayout(layout)
                self._label = QLabel(self._label_text, base)
                layout.addWidget(self._label)
            case "QGroupBox":
                base = QGroupBox(self._label_text, parent)
                base.setToolTip(self._tooltip)
                layout = layoutClass(base)
                base.setLayout(layout)
            case _:
                raise NotImplementedError(
                    f"Cannot create base of type {self._base_type}."
                )

        self._base = base
        self._base.setObjectName("InputWidget")
        self._layout = layout
        self._parent_dialog = parent
        self._empty = False
        self.has_warning = False

    def update_labels(self):
        """Update contained labels (dependent on base_type)."""

        match self._base_type:
            case "QWidget":
                self._label.setText(self._label_text)
            case "QGroupBox":
                self._base.setTitle(self._label_text)

        for widget in self._layout.children():
            if issubclass(widget, QWidget):
                widget.setToolTip(self._tooltip)

    def default_labels(self):
        """Provide default labels.

        Each Widget should have a default tooltip and label,
        which will be set in this method, unless specific
        values are provided in the settings of the job that
        is being configured.
        """
        if not self._label_text:
            self._label_text = "Base Widget"
        if not self._tooltip:
            self._tooltip = "A specific tooltip for this widget SHOULD have been set"

    @abstractmethod
    def value_from_configurator(self):
        """
        Set the widgets to the values of the underlying configurator object.

        Should also check for dependencies of the configurator.
        """

    @abstractmethod
    def get_widget_value(self):
        """Collect the results from the input widgets and return the value."""

    def mark_error(self, error_text: str, *, silent: bool = False):
        """Highlight an erroneous entry and display given error_text.

        Parameters
        ----------
        error_text : str
            Message displayed on hover-over.
        silent : bool
            If True, update the widget's error without sending signals

        """
        self._base.setStyleSheet(ERROR_STYLE)
        self._base.setToolTip(error_text)
        if not silent:
            self.valid_changed.emit()

    def mark_warning(self, warning_text: str):
        """If the input caused a warning, display warning_text and highlight the widget.

        If warning_text is an empty string, this method will clear errors instead.

        Parameters
        ----------
        warning_text : str
            Message displayed on hover-over.
        """
        if warning_text:
            self.has_warning = True
            self._base.setStyleSheet(WARNING_STYLE)
            self._base.setToolTip(warning_text)
            self.valid_changed.emit()
            return
        self.has_warning = False
        self.clear_error()

    def clear_error(self):
        """Remove error highlighting."""
        self._base.setStyleSheet("")
        self._base.setToolTip("")
        self.valid_changed.emit()

    def set_parameter(self) -> None:
        setattr(self._configurable, self._property, self.get_widget_value())

    @property
    def valid(self) -> bool:
        return getattr(self._configurable, self.parameter.configured_var, False)

    @property
    def raw(self) -> Any:
        return getattr(self._configurable, self.parameter.raw_name, None)

    @property
    def default(self):
        return self._configurable.default_parameters[self._property]

    @property
    def value(self):
        return getattr(self._configurable, self._property)

    @property
    def choices(self):
        if isinstance(self.parameter.choices, EnumMeta):
            option_list = (
                member.name.capitalize() for member in self.parameter.choices
            )
            return option_list

        if isinstance(self.parameter, CustomChoices):
            if self.parameter._bad_deps(self._configurable):
                self.mark_error("Invalid dependencies")
                return ()
            deps = self.parameter._get_deps(self._configurable)
            return self.parameter.get_choices(deps)

        if self.parameter.choices:
            return self.parameter.choices

    def get_widget_deps(self):
        return {
            key: self._parent._raw_widgets[val]
            for key, val in self.parameter.depends.items()
        }

    @Slot()
    def updateValue(self):
        try:
            self.set_parameter()
            self.clear_error()
        except ConfigError as err:
            self.mark_error(str(err))

        self.value_changed.emit()

    @property
    def default_path(self) -> str:
        """Default path of parent dialog."""
        return self._parent_dialog.default_path

    @default_path.setter
    def default_path(self, value: str) -> None:
        self._parent_dialog.default_path = value
