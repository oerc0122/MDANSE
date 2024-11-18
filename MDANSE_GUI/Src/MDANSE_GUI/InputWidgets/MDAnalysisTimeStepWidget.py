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
import MDAnalysis as mda

from MDANSE.MLogging import LOG
from .FloatWidget import FloatWidget


class MDAnalysisTimeStepWidget(FloatWidget):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # file input widgets should be loaded into the _widgets list
        # before this one
        for widget in self.parent()._widgets:
            if (
                widget._configurator
                is self._configurator._configurable[
                    self._configurator._dependencies["topology_file"]
                ]
            ):
                self._topology_file_widget = widget
                self._topology_file_widget.value_changed.connect(self.update_from_files)
            if (
                widget._configurator
                is self._configurator._configurable[
                    self._configurator._dependencies["coordinate_files"]
                ]
            ):
                self._coordinates_file_widget = widget
                self._coordinates_file_widget.value_changed.connect(
                    self.update_from_files
                )

    def update_from_files(self) -> None:
        """Updates the time step field from the topology and coordinates
        files if possible else set it back to the default value.
        """
        if (
            self._topology_file_widget._configurator.valid
            and self._coordinates_file_widget._configurator.valid
        ):
            try:
                coord_format = self._coordinates_file_widget._configurator["format"]
                coord_files = self._coordinates_file_widget._configurator["filenames"]

                if len(coord_files) <= 1 or coord_format is None:
                    value = mda.Universe(
                        self._topology_file_widget._configurator["filename"],
                        *coord_files,
                        format=coord_format,
                        topology_format=self._topology_file_widget._configurator[
                            "format"
                        ],
                    ).trajectory.ts.dt
                else:
                    coord_files = [(i, coord_format) for i in coord_files]
                    value = mda.Universe(
                        self._topology_file_widget._configurator["filename"],
                        coord_files,
                        topology_format=self._topology_file_widget._configurator[
                            "format"
                        ],
                    ).trajectory.ts.dt

                self._field.setText(str(value))
                return
            except Exception as e:
                LOG.warning(f"Failed to determine time step from MDAnalysis: {e}")

        self._field.setText(str(self._default_value))
