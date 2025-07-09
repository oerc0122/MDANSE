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


class ChoiceConfigDesc(ConfigureDescriptor):
    """
    Select an option from a multiple choice set.
    """

    def __init__(
        self,
        *,
        n_choices: int | None,
        **params,
    ):
        super().__init__(**params)

        self.n_choices = n_choices

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


class MultipleChoiceConfigDesc(MultipleValues, ChoiceConfigDesc):
    def __init__(self, choices: Container[T], n_choices: int | None, **params):
        super().__init__(choices=choices, n_choices=n_choices, **params)

    def validate(self, values: Sequence[T], *_) -> Sequence[T]:
        values = super().validate(values)

        if len(values) > self.n_choices:
            raise ConfigError(f"Too many options selected ({values}).")

        return values


class SingleChoiceConfigDesc(ChoiceConfigDesc):
    def __init__(self, choices: Container[T], **params):
        super().__init__(choices=choices, n_choices=1, **params)


class DynamicSingleChoiceConfigDesc(SingleChoiceConfigDesc):
    def __init__(self, choices: str, depends: dict[str, str], **params):
        if "choices" not in depends:
            raise ConfigError(f"{type(self).__name__} requires 'choices' in `depends`")
        super().__init__(choices=choices, depends=depends, **params)
        self.last_choices = set()

    @property
    def choices(self) -> set[T]:
        return self.last_choices

    @choices.setter
    def choices(self, value: str) -> None:
        self._choices = value

    def validate(self, value: T, deps: dict[str, Any]) -> T:
        choices = set(get_deep_attr(deps["choices"], self._choices))

        if not self.validate_exclude(value):
            raise ConfigError(
                f"Value ({value!r}) in excluded values ({', '.join(self.exclude)})."
            )

        if not self.validate_choices(value, choices):
            raise ConfigError(
                f"Value ({value!r}) not in choices ({', '.join(self.choices)})."
            )

        self.last_choices = choices
        return value

class DynamicMultiChoiceConfigDesc(MultipleValues, DynamicSingleChoiceConfigDesc):
    pass
