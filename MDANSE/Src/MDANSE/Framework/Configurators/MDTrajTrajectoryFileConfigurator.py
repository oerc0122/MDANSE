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
from pathlib import Path

from mdtraj.formats.registry import FormatRegistry


from .InputFileConfigurator import InputFileConfigurator


class MDTrajTrajectoryFileConfigurator(InputFileConfigurator):

    def configure(self, value):
        super().configure(value)

        extension = "".join(Path(self["value"]).suffixes)[1:]

        supported = list(i[1:] for i in FormatRegistry.loaders.keys())
        if extension not in supported:
            self.error_status = f"File '{extension}' not support should be one of the following: {supported}"
            return
