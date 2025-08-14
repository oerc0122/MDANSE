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

import json
from abc import ABC, abstractmethod
from collections.abc import Callable, Collection, Container, Iterator
from enum import Enum
from typing import Any, Generic, TypeVar
from warnings import Warning

import numpy as np
from more_itertools import value_chain

from MDANSE.Core.Error import Error
from MDANSE.IO.IOUtils import MDANSEEncoder

SENTINEL = object()
T = TypeVar("T")
CD = TypeVar("CD", bound="ConfigureDescriptor")


class ConfigError(Error):
    pass


class ConfigWarning(Warning):
    pass


class Parameter(ABC):
    """Abstract mixin for classes which can have a GUI Component."""

    default_label = ""
    default_tooltip = ""

    def __init__(
        self,
        *,
        label: str | None = None,
        tooltip: str | None = None,
    ):
        self.label = label or self.default_label
        self.tooltip = tooltip or self.default_tooltip

        self.__doc__ = self.__doc__ or ""
        self.__doc__ += f"""
GUI Description
---------------
{self.label}

{self.tooltip}
        """


class Configurable:
    """Allows any object that derives from it to be configurable within the MDANSE framework.

    Within that framework, to be configurable, a class must:
        - Derive from this class
    """

    @property
    def configuration(self) -> dict[str, Any]:
        return {name: getattr(self, name) for name in self.descriptors}

    @property
    def descriptors(self) -> dict[str, Parameter]:
        """Get all descriptors owned by self."""
        return {
            name: param
            for cls in value_chain(type(self).__bases__, type(self))
            for name, param in cls.__dict__.items()
            if isinstance(param, Parameter)
        }

    @property
    def parameters(self) -> list[str]:
        return list(self.descriptors)

    def to_json(self) -> str:
        to_json = {
            name: obj.to_json() if isinstance(obj, Parameter) else obj
            for name, obj in self.configuration.items()
        }
        return json.dumps(to_json, cls=MDANSEEncoder)

    output_configuration = to_json

    def check_status(self) -> bool:
        try:
            for param in self.parameters:
                getattr(self, param)
        except ConfigError:
            return False

        return True

    def __str__(self) -> str:
        out = f"{type(self).__name__}(\n"
        for param in self.parameters:
            try:
                val = getattr(self, param)
            except ConfigError:
                val = "Undefined or invalid"

            out += f"  {param} = {val},\n"
        out += ")"
        return out


class CustomConfig(Parameter, Configurable):
    """Abstract mixin for classes which contain Descriptors, but are themselves to be configured."""


class ConfigureDescriptor(Parameter, Generic[T]):
    """Abstract configure descriptor.

    Parameters
    ----------
    default : T
        Default value if not configured.
    optional : bool
        Whether this value can be disabled (None).
    choices : Container[T], optional
        Valid options for this value.
    exclude : Container[T]
        Invalid options for this value.
    mutex : Container[str]
        Other variables with which this is mutually exclusive.
    depends : Container[str]
        Other variables which must be set for this to be set.
    callback : Callable
        Custom function to be called post-validation.
    label : str
        GUI label to present to user.
    tooltip : str
        GUI tooltip to present to user.
    """

    def __init__(
        self,
        *,
        default: T = SENTINEL,
        optional: bool | None = None,
        choices: Container[T] | None = None,
        exclude: Container[T] = (),
        mutex: Container[str] = (),
        depends: dict[str, str] | None = None,
        callback: Callable[[CD, T, dict[str, str]], T] | None = None,
        label: str = "",
        tooltip: str = "",
        **params,
    ):
        self.default: T = default
        self.optional = optional if optional is not None else default is not SENTINEL

        self.choices = choices
        self.exclude = exclude

        self.mutex = mutex

        self.depends: dict[str, str] = depends if depends is not None else {}
        self.callback = callback

        if missing := (self.required_deps() - self.depends.keys()):
            raise ConfigError(f"Required deps ({', '.join(missing)}) missing.")

        super().__init__(label=label, tooltip=tooltip)

    def __set_name__(self, owner: object, name: str):
        self.name = name
        self.private_name = "_" + name
        self.configured_var = self.private_name + "_configured"
        setattr(owner, self.configured_var, False)

    def __get__(self, owner: object, objtype: type | None = None) -> T:
        if owner is None:
            return self

        if not self.optional and not getattr(owner, self.configured_var):
            raise ConfigError(f"Non-optional value ({self.name}) has not been set")

        value = getattr(owner, self.private_name, SENTINEL)
        return value if value is not SENTINEL else self.default

    def _bad_mutex(self, owner: object) -> Iterator[bool]:
        return (getattr(owner, f"_{ex}_configured") for ex in self.mutex)

    def _bad_deps(self, owner: object) -> Iterator[bool]:
        return (
            not getattr(owner, f"_{dep}_configured") for dep in self.depends.values()
        )

    def __set__(self, owner: object, value):
        setattr(owner, self.configured_var, False)
        setattr(owner, self.private_name, SENTINEL)

        if any(self._bad_deps(owner)):
            raise ConfigError(
                f"Dependencies ({', '.join(self._bad_deps)}) are not correctly defined."
            )
        if any(self._bad_mutex(owner)):
            raise ConfigError(
                f"Mutually exclusive value ({', '.join(self._bad_mutex)}) is also configured."
            )

        deps = {dep: getattr(owner, key) for dep, key in self.depends.items()}
        value = self.validate(value, deps)

        # Custom
        if self.callback is not None:
            value = self.callback(self, value, deps)

        # For grouped config descriptors
        if hasattr(owner, "validate"):
            value = owner.validate(self, value)

        setattr(owner, self.private_name, value)
        setattr(owner, self.configured_var, True)

    def required_deps(self) -> set[str]:
        return set()

    @property
    def choices(self) -> set[T]:
        """
        Returns the set of values allowed for an input.
        """
        return self._choices

    @choices.setter
    def choices(self, value: Container[T] | Enum) -> None:
        if value is None:
            self._choices = set()
            return
        if isinstance(value, Enum):
            self._choices = value
            return
        self._choices = set(value)
        self.last_choices = self._choices

    @property
    def exclude(self) -> set[int]:
        """
        Returns the set of values which are not forbidden.

        Returns
        -------
        set[int]
            Forbidden values.
        """
        return self._exclude

    @exclude.setter
    def exclude(self, value: Container[int]) -> None:
        if value is None:
            self._exclude = set()
            return
        self._exclude = set(value)

    def validate_choices(self, value: T, choices: set[T] | Enum | None = None) -> bool:
        choices = self.choices if choices is None else choices
        return not choices or value in choices

    def validate_exclude(self, value: T, exclude: set[T] | None = None) -> bool:
        exclude = self.exclude if exclude is None else exclude
        return not exclude or value not in exclude

    @abstractmethod
    def validate(self, value: T, deps: dict[str, Any]) -> T:
        """
        Ensure that the passed variable is of the right type.
        """
        # Choices
        if self.choices and isinstance(self.choices, Enum):  # Enum
            try:
                value = self.choices(value)
            except ValueError:
                raise ConfigError(
                    f"Value ({value!r}) not in choices ({', '.join(self.choices)})."
                )

        elif hasattr(self, "get_choices"):  # Function
            self.last_choices = self.get_choices(deps)

            if not self.validate_choices(value, self.last_choices):
                raise ConfigError(
                    f"Value ({value!r}) not in choices ({', '.join(self.last_choices)})."
                )

        elif not self.validate_choices(value):  # Set
            raise ConfigError(
                f"Value ({value!r}) not in choices ({', '.join(self.choices)})."
            )

        # Exclude
        if not self.validate_exclude(value):
            raise ConfigError(
                f"Value ({value!r}) in excluded values ({', '.join(self.exclude)})."
            )

        return value


class MinMax(ABC, Generic[T]):
    def __init__(
        self, minimum: T | None = None, maximum: T | None = None, *args, **kwargs
    ):
        self.minimum = minimum
        self.maximum = maximum
        super().__init__(*args, **kwargs)

    @property
    def minimum(self) -> T | None:
        """
        Returns the minimum value allowed for an input int.

        Returns
        -------
        int or None
            The minimum value allowed for an input value int.
        """
        return self._minimum

    @minimum.setter
    def minimum(self, value: T | None) -> None:
        self._minimum = value if value is not None else -np.inf

    @property
    def maximum(self) -> T | None:
        """
        Returns the maximum value allowed for an input int.

        Returns
        -------
        int or None
            The maximum value allowed for an input value int.
        """
        return self._maximum

    @maximum.setter
    def maximum(self, value: T | None) -> None:
        self._maximum = value if value is not None else np.inf

    def get_ranges(self, value: T, deps: dict[str, Any]):
        return self.minimum, self.maximum

    def validate_range(self, value: T, ranges: tuple[T, T] | None = None) -> None:
        if ranges is not None:
            mini, maxi = ranges
        else:
            mini, maxi = self.minimum, self.maximum

        if mini > value > maxi:
            raise ConfigError(
                f"Value ({value}) outside of valid range ({mini}, {maxi})"
            )

    def validate_minimum(self, value: T, mini: T | None = None) -> None:
        mini = mini if mini is not None else self.minimum

        if mini > value:
            raise ConfigError(f"Value ({value}) less than minimum ({mini})")

    def validate_maximum(self, value: T, maxi: T | None = None) -> None:
        maxi = maxi if maxi is not None else self.maximum

        if maxi > value:
            raise ConfigError(f"Value ({value}) greater than maximum ({maxi})")


class MultipleValues(ABC):
    def validate_choices(
        self, values: Collection[T], choices: set[T] | None = None
    ) -> bool:
        choices = self.choices if choices is None else choices
        return set(values) <= choices

    def validate_exclude(
        self, values: Collection[T], exclude: set[T] | None = None
    ) -> bool:
        exclude = self.exclude if exclude is None else exclude
        return set(values) > exclude
