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

from collections.abc import Container, Sequence
from typing import Any, TypeVar

from MDANSE.Core.get_deep_attr import get_deep_attr

from .AbsConfigDesc import ConfigError, ConfigureDescriptor, MultipleValues

T = TypeVar("T")


class Choice(ConfigureDescriptor):
    """
    Select an option from a multiple choice set.
    """

    def __init__(
        self,
        *,
        n_choices: int | None,
        aliases: dict[Any, T] | None = None,
        **params,
    ):
        super().__init__(**params)

        self.n_choices = n_choices
        self.aliases = aliases if aliases is not None else {}

    @property
    def n_choices(self) -> int:
        """Returns the minimum value allowed for an input float.

        Returns
        -------
        int
            the minimum value allowed for an input value float.
        """
        return self._n_choices if self._n_choices is not None else len(self.choices)

    @n_choices.setter
    def n_choices(self, value: int | None):
        if value is None:
            self._n_choices = value
            return

        if value < 0:
            raise ConfigError(f"Invalid n_choices ({value}) must be >0")
        self._n_choices = int(value)


class MultipleChoice(MultipleValues, Choice):
    def __init__(self, choices: Container[T], n_choices: int | None, **params):
        super().__init__(choices=choices, n_choices=n_choices, **params)

    def validate(self, values: Sequence[T], *_) -> Sequence[T]:
        values = [self.aliases.get(value, value) for value in values]
        values = super().validate(values)

        if len(values) > self.n_choices:
            raise ConfigError(f"Too many options selected ({values}).")

        return values


class SingleChoice(Choice):
    def __init__(self, choices: Container[T], **params):
        super().__init__(choices=choices, n_choices=1, **params)

    def validate(self, value: T, *_) -> T:
        value = self.aliases.get(value, value)
        return super().validate(value)


class DynamicSingleChoice(SingleChoice):
    def __init__(self, choices: str, depends: dict[str, str], **params):
        super().__init__(choices=choices, depends=depends, **params)
        self.last_choices = set()

    def required_deps(self) -> set[str]:
        return super().required_deps() | {"choices"}

    @property
    def choices(self) -> set[T]:
        return self.last_choices

    @choices.setter
    def choices(self, value: str) -> None:
        self._choices = value

    def get_choices(self, deps) -> set[T]:
        choices = set(get_deep_attr(deps["choices"], self._choices))
        if not choices:
            raise ConfigError(f"No valid choices at deps['choices'].{self._choices}")

        return choices


class DynamicMultiChoice(MultipleValues, DynamicSingleChoice):
    pass
