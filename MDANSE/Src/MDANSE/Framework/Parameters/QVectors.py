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

from collections.abc import Callable
from typing import Any

from MDANSE.Framework.QVectors.IQVectors import IQVectors

from .BaseTypes import Dict
from .Choices import SingleChoice
from .Parameters import ConfigError
from .UtilTypes import Depends, DescID


class QVectorsSelect(SingleChoice[str | IQVectors, IQVectors]):
    def __init__(self, *args, choices: None = None, aliases: None = None, **kwargs):
        if choices is not None:
            raise ConfigError(
                f"Not allowed to provide choices to {type(self).__name__}."
            )
        if aliases is not None:
            raise ConfigError(
                f"Not allowed to provide aliases to {type(self).__name__}."
            )

        super().__init__(
            *args,
            choices=IQVectors.indirect_subclass_dictionary().values(),
            aliases=IQVectors.indirect_subclass_dictionary(),
            **kwargs,
        )

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def validate(self, value: str | type[IQVectors], deps: Depends, /) -> IQVectors:
        value: type[IQVectors] = super().validate(value, deps)

        return value(deps["trajectory"])


class QVectors(Dict[str, Any]):
    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("generator")}

    def get_choices(self, deps: Depends):
        return deps["generator"].parameters

    def validate(self, value: dict[str, Any], deps: Depends, /) -> IQVectors:
        value = super().validate(value, deps)

        gen = deps["generator"]
        gen.configuration = value
        return gen
