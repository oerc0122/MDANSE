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
from typing import TYPE_CHECKING, Any, cast

from MDANSE.IO.IOUtils import json_handler
from MDANSE.Mathematics.Signal import (
    DEFAULT_FILTER,
    FILTER_MAP,
    FILTERS,
    filter_default_attributes,
)
from MDANSE.Mathematics.Signal import Filter as SFilter

from .AbsConfigDesc import ConfigError, CustomConfig
from .BaseTypesConfigDesc import Dict
from .ChoiceConfigDesc import SingleChoice

if TYPE_CHECKING:
    FilterType = type[SFilter]


class FilterParamDict(Dict[str, Any]):
    def __init__(self, *args, fixed: None = None, **kwargs):
        super().__init__(*args, fixed=True, **kwargs)

    def required_deps(self) -> set[str]:
        return super().required_deps() | {"filter"}

    def get_choices(self, value, deps: dict[str, Any]):
        return deps["filter"].default_settings.keys()


class Filter(CustomConfig):
    default_tooltip = "Choose trajectory filter."
    default_label = "Trajectory filter"

    filter_type = SingleChoice[FilterType](
        choices=FILTERS,
        aliases=FILTER_MAP,
    )
    params = FilterParamDict(depends={"filter": "filter_type"})

    def __init__(
        self,
        filter: FilterType | str = DEFAULT_FILTER,
        params: dict[str, Any] | None = None,
    ):
        if filter not in FILTER_MAP:
            raise ConfigError("Config not in keys.")

        if isinstance(filter, str):
            filter = FILTER_MAP[filter]

        self.filter_type = filter

        self.params = params or filter_default_attributes(self.filter_type)

    @property
    def filter(self) -> FilterType:
        return self.filter_type(**self.params)

    @filter.setter
    def filter(self, value: dict | str | Path | tuple[str, dict[str, Any]] | SFilter):
        if isinstance(value, (dict, str, Path)):
            try:
                value = json_handler(value)
            except Exception as err:
                raise ConfigError("Couldn't set filter from json") from err
            self.filter_type = value.pop("filter")
            self.params = value
        elif isinstance(value, tuple):
            self.filter_type, self.params = value
        elif isinstance(value, SFilter):
            self.filter_type = type(value).__name__
            self.params = {
                key: val for key, val in value.__dict__ if key in value.default_settings
            }
