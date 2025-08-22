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
from typing import TYPE_CHECKING, Any, ParamSpec, cast

from MDANSE.IO.IOUtils import json_handler
from MDANSE.Mathematics.Signal import (
    DEFAULT_FILTER,
    FILTER_MAP,
    FILTERS,
    Filter,
    filter_default_attributes,
)

from .AbsConfigDesc import ConfigError, CustomConfig
from .BaseTypesDescriptor import Dict
from .ChoiceConfigDesc import SingleChoice
from .UtilTypes import Depends, DescID

FilterType = type[Filter]


class FilterParamDict(Dict[str, Any]):
    def __init__(self, *args, fixed: None = None, **kwargs):
        super().__init__(*args, fixed=True, **kwargs)
        self.last_choices = set()

    @property
    def choices(self) -> set[str]:
        return self.last_choices

    @choices.setter
    def choices(self, value) -> None: ...

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("filter")}

    def get_choices(self, deps: Depends, /):
        return deps["filter"].default_settings.keys()


class TrajectoryFilter(CustomConfig):
    default_tooltip = "Choose trajectory filter."
    default_label = "Trajectory filter"

    filter_type = SingleChoice[str | FilterType, FilterType](
        choices=FILTERS,
        aliases=FILTER_MAP,
    )
    params = FilterParamDict(depends={"filter": "filter_type"})

    def __init__(
        self,
        filter: FilterType | str = DEFAULT_FILTER,
        params: dict[str, Any] | None = None,
    ):
        filter = FILTER_MAP.get(filter, filter)
        if filter not in FILTERS:
            raise ConfigError("Config not in keys.")

        if isinstance(filter, str):
            filter = FILTER_MAP[filter]

        self.filter_type = filter

        self.params = params or filter_default_attributes(self.filter_type)

    @property
    def filter(self) -> Filter:
        filt = cast(FilterType, self.filter_type)
        return filt(**self.params)

    @filter.setter
    def filter(self, value: dict | str | Path | tuple[str, dict[str, Any]] | Filter):
        if isinstance(value, (dict, str, Path)):
            try:
                value = json_handler(value)
            except Exception as err:
                raise ConfigError("Couldn't set filter from json") from err
            self.filter_type = value.pop("filter")
            self.params = value
        elif isinstance(value, tuple):
            self.filter_type, self.params = value
        elif isinstance(value, Filter):
            self.filter_type = type(value).__name__
            self.params = {
                key: val for key, val in value.__dict__ if key in value.default_settings
            }
