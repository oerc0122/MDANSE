#    This file is part of MDANSE.
#
#    MDANSE is free software: you can redistribute it and/or modify
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

from MDANSE.Framework.Configurators.FloatConfigurator import FloatConfigurator


class MDAnalysisTimeStepConfigurator(FloatConfigurator):

    _default = 0.0

    def configure(self, value):
        # if the value is not valid then we use the MDAnalysis
        # default values which maybe the time step in the inputted
        # files or 1 ps
        try:
            value = float(value)
        except (TypeError, ValueError):
            pass

        if value is None or value == "" or value == 0.0:
            file_configurator = self._configurable[self._dependencies["topology_file"]]
            files_configurator = self._configurable[
                self._dependencies["coordinate_files"]
            ]
            if file_configurator._valid and files_configurator._valid:
                try:
                    value = mda.Universe(
                        file_configurator["filename"], *files_configurator["filenames"]
                    ).trajectory.ts.dt
                except (AttributeError, ValueError) as e:
                    self.error_status = (
                        f"Unable to determine a time step from MDAnalysis: {e}"
                    )
                    return
            else:
                self.error_status = f"Unable to determine a time step from MDAnalysis"
                return

        super().configure(value)
