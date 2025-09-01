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
from enum import EnumMeta
from typing import Any, TypeVar, cast

from MDANSE.Core.get_deep_attr import get_deep_attr

from .Parameters import ConfigError, ConfigureDescriptor
from .UtilTypes import Depends, DescID

P = TypeVar("P")
T = TypeVar("T")


class Choice(ConfigureDescriptor[P, T]):
    """
    Select an option from a multiple choice set.
    """

    def __init__(
        self,
        *,
        n_choices: int | None,
        aliases: dict[Any, P] | None = None,
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


class MultipleChoice(Choice[Sequence[P], Sequence[T]]):
    def __init__(self, *args, choices: Sequence[T], n_choices: int | None, **params):
        super().__init__(*args, choices=choices, n_choices=n_choices, **params)

    def _validate_choices(
        self, value: Sequence[T], choices: set[T] | None = None
    ) -> bool:
        choice = self.choices if choices is None else choices
        assert not isinstance(choice, EnumMeta)
        return set(value) <= choice

    # def _validate_exclude(
    #     self, value: Sequence[T], exclude: set[T] | None = None
    # ) -> bool:
    #     excludes = self.exclude if exclude is None else exclude
    #     return set(value) > excludes

    def validate(self, values: Sequence[P], deps: Depends, /) -> Sequence[T]:
        dealiased = cast(list[P], [self.aliases.get(value, value) for value in values])
        out = super().validate(dealiased, deps)

        if len(values) > self.n_choices:
            raise ConfigError(f"Too many options selected ({values}).")

        return out


class SingleChoice(Choice[P, T]):
    def __init__(self, choices: Container[T], n_choices: None = None, **params):
        if n_choices is not None:
            raise ConfigError(f"Cannot define n_choices in {type(self).__name__}")
        super().__init__(choices=choices, n_choices=1, **params)

    def validate(self, value: P, deps: Depends, /) -> T:
        value = self.aliases.get(value, value)
        return super().validate(value, deps)


class DynamicSingleChoice(SingleChoice[P, T]):
    def __init__(self, choices: str, **params):
        super().__init__(choices=(), **params)
        self.getter = choices
        self.last_choices = set()

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("choices")}

    @property
    def choices(self) -> set[T]:
        return self.last_choices

    @choices.setter
    def choices(self, value: Any) -> None:
        self._choices = set()

    def get_choices(self, deps) -> set[T]:
        choices = set(get_deep_attr(deps["choices"], self.getter))
        if not choices:
            raise ConfigError(f"No valid choices at deps['choices'].{self.getter}")

        return choices


class DynamicMultiChoice(MultipleChoice[Sequence[P], Sequence[T]]):
    def __init__(self, choices: str, n_choices: int | None, **params):
        super().__init__(choices=(), n_choices=n_choices, **params)
        self.getter = choices
        self.last_choices = set()

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("choices")}

    @property
    def choices(self) -> set[T]:
        return self.last_choices

    @choices.setter
    def choices(self, value: Any) -> None:
        self._choices = set()

    def get_choices(self, deps) -> set[T]:
        choices = set(get_deep_attr(deps["choices"], self.getter))
        if not choices:
            raise ConfigError(f"No valid choices at deps['choices'].{self.getter}")

        return choices
