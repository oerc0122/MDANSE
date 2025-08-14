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

from typing import Any

import h5py

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .AbsConfigDesc import ConfigError, ConfigureDescriptor
from .BaseTypesDescriptor import PathParam


class MDANSETrajectory(ConfigureDescriptor[Trajectory]):
    default_tooltip = "Input MDANSE trajectory from converter or h5md file."
    default_label = "Trajectory to use."

    def __init__(
        self,
        **params,
    ):
        super().__init__(**params)
        self.extension = {
            "MDANSE trajectory": "mdt",
            "HDF5 file": "h5",
            "All files": "*",
        }

    def validate_choices(
        self,
        value: Trajectory,
        choices: set[str] | None = None,
    ) -> bool:
        choices = self.choices if choices is None else choices
        return all(value.has_variable(choice) for choice in choices)

    def validate_exclude(
        self,
        value: Trajectory,
        exclude: set[str] | None = None,
    ) -> bool:
        excludes = self.exclude if exclude is None else exclude
        return not any(value.has_variable(exclude) for exclude in excludes)

    def validate(self, value, deps: dict[str, Any]) -> Trajectory:
        value = PathParam.validate(self, value, deps)

        try:
            out = Trajectory(value)
        except Exception:
            raise ConfigError(f"Could not load file ({value}) as Trajectory")

        out = super().validate(out)

        return out


class MDANSEResult(ConfigureDescriptor[h5py.File]):
    default_tooltip = "Input MDANSE result from prior analysis."

    def __init__(
        self,
        **params,
    ):
        super().__init__(**params)
        self.extension = {"MDANSE result": "mda", "HDF5 file": "h5", "All files": "*"}

    def validate_choices(
        self,
        value: h5py.File,
        choices: set[str] | None = None,
    ) -> bool:
        choices = self.choices if choices is None else choices
        return all(choice in value for choice in choices)

    def validate_exclude(
        self,
        value: h5py.File,
        exclude: set[str] | None = None,
    ) -> bool:
        excludes = self.exclude if exclude is None else exclude
        return not any(exclude in value for exclude in excludes)

    def validate(self, value, deps: dict[str, Any]) -> Trajectory:
        value = PathParam.validate(self, value, deps)

        try:
            out = h5py.File(value)
        except Exception:
            raise ConfigError(f"Could not load file ({value}) as hdf5.")

        out = super().validate(out, deps)

        return out
