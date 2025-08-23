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

from pathlib import Path
from typing import Any

import h5py

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .AbsConfigDesc import ConfigError, ConfigureDescriptor
from .BaseTypesDescriptor import PathParam
from .UtilTypes import Depends, DescID


def set_up_trajectory(self, trajectory: Trajectory, deps):
    """Apply operations to the trajectory instance, if present.

    Atom selection, atom transmutation and result grouping are all
    applied to the Trajectory object. If the job works on a trajectory,
    the Trajectory instance is now saved as an attribute of this IJob
    instance.

    These operations were previously handled by IConfigurator subclasses.
    """
    if (selection := deps.get("atom_selection")) is not None:
        trajectory.set_selection(selection)
    if (transmutation := deps.get("atom_transmutation")) is not None:
        trajectory.set_transmutation(transmutation.transmutation)
    if (grouping := deps.get("grouping_level")) is not None:
        trajectory.set_grouping(grouping["level"])

    return trajectory

class MDANSETrajectory(ConfigureDescriptor[str | Path, Trajectory]):
    default_tooltip = "Input MDANSE trajectory from converter or h5md file."
    default_label = "Trajectory to use."

    def __init__(
        self,
        *,
        selection: str | None = None,
        transmutation: str | None = None,
        grouping: str | None = None,
        on_get_depends: None = None,
        on_get: None = None,
        **params,
    ):
        if on_get_depends is not None:
            raise ConfigError(f"Cannot set `on_get_depends` in {type(self).__name__}.")
        if on_get is not None:
            raise ConfigError(f"Cannot set `on_get` in {type(self).__name__}.")

        ogd = {
            "atom_selection": selection,
            "atom_transmutation": transmutation,
            "grouping_level": grouping,
        }

        super().__init__(on_get_depends=ogd, on_get=set_up_trajectory, **params)
        self.extension = {
            "MDANSE trajectory": "mdt",
            "HDF5 file": "h5",
            "All files": "*",
        }

    def _validate_choices(
        self,
        value: Trajectory,
        choices: set[str] | None = None,
    ) -> bool:
        choice_in = self.choices if choices is None else choices
        return all(value.has_variable(choice) for choice in choice_in)

    def _validate_exclude(
        self,
        value: Trajectory,
        exclude: set[str] | None = None,
    ) -> bool:
        exclude_in = self.exclude if exclude is None else exclude
        return not any(value.has_variable(exclude) for exclude in exclude_in)

    def validate(self, value, deps: Depends) -> Trajectory:
        try:
            value = Path(value).expanduser()
        except TypeError as error:
            raise ConfigError(f"Value ({value}) is not a valid Path.") from error

        if not value.exists():
            raise ConfigError(f"File at ({value}) does not exist.")

        try:
            out = Trajectory(value)
        except Exception:
            raise ConfigError(f"Could not load file ({value}) as Trajectory")

        out = super().validate(out, deps)

        return out


class MDANSEResult(ConfigureDescriptor[str | Path | h5py.File, h5py.File]):
    default_tooltip = "Input MDANSE result from prior analysis."

    def __init__(
        self,
        **params,
    ):
        super().__init__(**params)
        self.extension = {"MDANSE result": "mda", "HDF5 file": "h5", "All files": "*"}

    def _validate_choices(
        self,
        value: h5py.File,
        choices: set[str] | None = None,
    ) -> bool:
        choice_in = self.choices if choices is None else choices
        return all(choice in value for choice in choices)

    def _validate_exclude(
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
