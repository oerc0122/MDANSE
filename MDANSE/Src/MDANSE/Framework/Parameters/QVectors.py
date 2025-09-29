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

from MDANSE.Framework.QVectors.IQVectors import IQVectors

from .BaseTypes import Dict
from .Choices import SingleChoice
from .Parameters import ConfigError, CustomConfig, HasDependencies
from .UtilTypes import Depends, DescID


class QVectorsSelect(SingleChoice[str | IQVectors, IQVectors]):
    def __init__(
        self,
        *args,
        choices: None = None,
        aliases: None = None,
        on_get: None = None,
        **kwargs,
    ):
        if choices is not None:
            raise ConfigError(
                f"Not allowed to provide choices to {type(self).__name__}."
            )
        if aliases is not None:
            raise ConfigError(
                f"Not allowed to provide aliases to {type(self).__name__}."
            )
        if on_get is not None:
            raise ConfigError(
                f"Not allowed to provide on_get to {type(self).__name__}."
            )

        super().__init__(
            *args,
            choices=IQVectors.indirect_subclass_dictionary().values(),
            aliases=IQVectors.indirect_subclass_dictionary(),
            on_get=self._generate,
            **kwargs,
        )

    @staticmethod
    def _generate(_desc, out, _deps) -> IQVectors:
        try:
            out.generate()
        except Exception as err:
            raise ConfigError("Invalid or incomplete QVector configuration.") from err
        return out

    def __set_name__(self, owner: type, name: str):
        self.name = name
        self.private_name = "_" + name
        self.configured_var = self.private_name + "_configured"

        setattr(owner, self.configured_var, self.optional)

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def validate(self, value: str | type[IQVectors], deps: Depends, /) -> IQVectors:
        value: type[IQVectors] = super().validate(value, deps)
        return value(deps["trajectory"].configuration(0))


class QVectorsParams(HasDependencies):
    def __set__(self, owner: object, value):
        generator = getattr(owner, self.depends["generator"])
        generator.configuration = value

    def __get__(self, owner: object, objtype: type | None):
        return getattr(owner, self.depends["generator"]).configuration

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("generator")}


class QVectors(CustomConfig):
    generator = QVectorsSelect(depends={"trajectory": "parent"})
    params = QVectorsParams(depends={"generator": "generator"})

    def required_deps(self) -> set[str]:
        return super().required_deps() | {"trajectory"}

    def __set__(self, owner: object, value):
        self.last_deps = self._get_deps(owner)

        if isinstance(value, tuple):
            self.generator, self.params = value
            return

        if isinstance(value, str):
            self.generator = value
            return

        if isinstance(value, dict):
            self.params = value
            return

        self.configuration = value.configuration


    def preview_output_axis(self):
        """Output the values of |Q| from current parameters.

        Returns
        -------
        list[float]
            Values of |Q|.
        str
            Physical unit of Q.

        """
        return self.generator.preview_output_axis()
