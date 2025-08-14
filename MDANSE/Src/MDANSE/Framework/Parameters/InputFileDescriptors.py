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
from __future__ import annotations

import h5py

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .BaseTypesDescriptor import PathParam


class MDANSETrajectory(PathParam):
    def __init__(
        self,
        mode: None = None,
        *,
        extensions: None = None,
        directory: None = None,
        **params,
    ):
        super().__init__(
            mode="r",
            extensions={"MDANSE trajectory": "mdt", "HDF5 file": "h5"},
            **params,
        )

    def validate(self, value, *_) -> Trajectory:
        value = super().validate(value)

        return Trajectory(value)
