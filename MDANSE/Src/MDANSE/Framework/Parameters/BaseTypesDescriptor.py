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

import re
from collections.abc import Sequence
from enum import Enum, auto
from math import isclose
from numbers import Number
from pathlib import Path
from typing import Any, SupportsFloat, SupportsInt

import numpy as np
from more_itertools import numeric_range

from MDANSE.IO.IOUtils import json_handler

from .AbsConfigDesc import ConfigError, ConfigureDescriptor, MinMax

cjoin = ", ".join

__all__ = [
    "Boolean",
    "Float",
    "Integer",
    "String",
    "PathParam",
    "Range",
    "NumericRange",
    "Array",
    "Vector",
]


def _nop(x):
    return x


class StrCases(Enum):
    PRESERVE = _nop
    UPPER = str.upper
    LOWER = str.lower
    TITLE = str.title

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, str):
            value = value.upper()

        for member in cls:
            if member.name == value or member.value is value:
                return member
        return NotImplemented


class Boolean(ConfigureDescriptor[bool]):
    """
    This Configurator allows to input a Boolean Value (True or False).

    The input value can be directly provided as a Python boolean or by the using the following (standard)
     representation of a boolean: 'true'/'false', 'yes'/'no', 'y'/'n', '1'/'0', 1/0
    """

    _alias = {
        True: True,
        "true": True,
        "yes": True,
        "y": True,
        "1": True,
        False: False,
        "false": False,
        "no": False,
        "n": False,
        "0": False,
    }

    def __init__(self, default: bool = False, **params):
        super().__init__(default=default, **params)

    def validate(self, value: bool | str, deps: dict[str, Any]) -> bool:
        if isinstance(value, str):
            value = value.lower()

        if value not in self._alias:
            raise ConfigError(f"Unable to convert value ({value!r}) to bool.")

        return self._alias[value]


class Float(MinMax[float], ConfigureDescriptor[float]):
    """
    This configurator allows to input a float value.
    """

    def __init__(
        self,
        non_zero: float | None = None,
        **params,
    ):
        super().__init__(**params)
        self.non_zero = non_zero

    def validate(self, value: SupportsFloat, deps: dict[str, Any]) -> float:
        try:
            value = float(value)
        except ValueError as error:
            raise ConfigError(f"Value ({value}) is not a valid float.") from error

        super().validate(value)

        ranges = self.get_ranges(self, value, deps)

        self.validate_range(value, ranges)

        if self.non_zero is not None and isclose(value, 0.0, abs_tol=self.non_zero):
            raise ConfigError(f"Non-null val ({value}) is too close to zero.")

        return value


class Integer(MinMax[int], ConfigureDescriptor[int]):
    """Configurator takes an integer input."""

    def __init__(self, **params):
        super().__init__(**params)

    def validate(self, value: SupportsInt, deps: dict[str, Any]) -> int:
        try:
            value = int(value)
        except ValueError as error:
            raise ConfigError(f"Value ({value}) is not a valid integer.") from error

        ranges = self.get_ranges(self, value, deps)

        self.validate_range(value, ranges)

        self.validate_range(value)

        return value


class Dict(ConfigureDescriptor[dict]):
    """
    This configurator allows a user to input a dict value.

    `choices` are keys which *must* be present in the dict.

    `exclude` are keys which *must not* be present in the dict.
    """

    def __init__(self, *, fixed: bool = False, **params):
        super().__init__(**params)

        self.fixed = fixed

    def validate(self, value, deps: dict[str, Any]):
        if isinstance(value, Sequence):
            value = dict(value)

        try:
            json_handler(value)
        except Exception as err:
            raise ConfigError("Failed to build dict.") from err

        if missing := self.choices - value.keys():
            raise ConfigError(f"Required keys ({cjoin(missing)}) are missing.")
        if erroneous := self.exclude & value.keys():
            raise ConfigError(f"Forbidden keys ({cjoin(erroneous)}) are present.")

        if self.fixed and (extraneous := value.keys() - self.choices):
            raise ConfigError(
                f"Extra keys are present ({cjoin(extraneous)}), "
                f"expected: {cjoin(self.choices)}."
            )

        return value


class String(ConfigureDescriptor[str]):
    """
    This configurator allows to input a string value.
    """

    def __init__(
        self,
        *,
        case: StrCases = StrCases.PRESERVE,
        regex: str | re.Pattern | None = None,
        **params,
    ):
        super().__init__(**params)
        self.case = StrCases(case)
        self.regex = regex

    @property
    def regex(self) -> re.Pattern | None:
        return self._regex

    @regex.setter
    def regex(self, value: re.Pattern | None) -> None:
        if value is None:
            self._regex = None
            return
        self.regex = re.compile(value)

    def validate(self, value: str, deps: dict[str, Any]) -> float:
        try:
            value = str(value)
        except ValueError as error:
            raise ConfigError(f"Value ({value}) is not a valid str.") from error

        value = self.case.value(value)

        super().validate(value)

        if self.regex is not None and not self.regex.match(value):
            raise ConfigError(
                f"Value ({value}) does not match requirement {self.regex.pattern}."
            )

        return value


class PathParam(ConfigureDescriptor[Path]):
    """
    This configurator allows to input a path value.
    """

    class FileModes(Enum):
        MUST_EXIST = auto()
        MAY_EXIST = auto()
        MUST_NOT_EXIST = auto()

        R = MUST_EXIST
        W = MAY_EXIST
        X = MUST_NOT_EXIST
        A = MUST_EXIST

    def __init__(
        self,
        mode: FileModes | str,
        *,
        extensions: dict[str, str] | None = None,
        directory: bool = False,
        **params,
    ):
        super().__init__(**params)

        self.mode = (
            self.FileModes[mode.upper()]
            if isinstance(mode, str)
            else self.FileModes(mode)
        )

        self.extension = extensions or {}
        self.extension.setdefault("All files", "*")
        self.directory = directory

    def validate(self, value, deps: dict[str, Any]) -> Path:
        try:
            value = Path(value).expanduser()
        except TypeError as error:
            raise ConfigError(f"Value ({value}) is not a valid Path.") from error

        super().validate(value)

        if self.mode is self.FileModes.MAY_EXIST:
            pass
        elif self.mode is self.FileModes.MUST_EXIST and not value.exists():
            raise ConfigError(f"File at ({value}) does not exist.")
        elif self.mode is self.FileModes.MUST_NOT_EXIST and value.exists():
            raise ConfigError(f"File at ({value}) must not exist.")

        return value


class Range(MinMax[int], ConfigureDescriptor[range]):
    """
    This configurator allows a user to input a range value.
    """

    def __init__(
        self,
        *args,
        include_last: bool = False,
        **params,
    ):
        super().__init__(*args, **params)
        self.include_last = include_last

    def validate(
        self, value: int | Sequence | dict | range, deps: dict[str, Any]
    ) -> range:
        if isinstance(value, int):
            value = range(value)
        elif isinstance(value, Sequence):
            value = range(*value)
        elif isinstance(value, dict):
            value = range(**value)
        elif not isinstance(value, range):
            raise ConfigError(
                f"Do not know how to convert {type(value).__name__} to range"
            )

        if self.include_last:
            value = range(value.start, value.stop + value.step, value.step)

        ranges = self.get_ranges(value, deps)

        self.validate_range(value.start, ranges)
        self.validate_range(value.stop, ranges)

        return value


class NumericRange(MinMax[Number], ConfigureDescriptor[numeric_range]):
    """
    This configurator allows a user to input a generalised range.
    """

    def __init__(
        self,
        *args,
        include_last: bool = False,
        **params,
    ):
        super().__init__(*args, **params)
        self.include_last = include_last

    def validate(
        self,
        value: float | int | Sequence | dict | range | numeric_range,
        deps: dict[str, Any],
    ) -> numeric_range:
        if isinstance(value, (float, int)):
            value = numeric_range(value)
        elif isinstance(value, Sequence):
            value = numeric_range(*value)
        elif isinstance(value, dict):
            value = numeric_range(**value)
        elif isinstance(value, range):
            value = numeric_range(value.start, value.stop, value.step)
        elif not isinstance(value, numeric_range):
            raise ConfigError(
                f"Do not know how to convert {type(value).__name__} to range"
            )

        assert value is numeric_range

        if self.include_last:
            value = numeric_range(value._start, value._stop + value._step, value._step)

        ranges = self.get_ranges(value, deps)

        self.validate_range(value._start, ranges)
        self.validate_range(value._stop, ranges)

        return value


class Array(ConfigureDescriptor[np.ndarray]):
    def __init__(
        self,
        *,
        shape: tuple[int, ...],
        dtype: np.dtype = np.float64,
        non_zero: bool = False,
        **params,
    ):
        super().__init__(**params)

        self.shape = shape
        self.non_zero = non_zero
        self.dtype = dtype

    def validate(self, value, deps: dict[str, Any]) -> np.ndarray:
        if isinstance(value, Sequence):
            value = np.array(value, dtype=self.dtype)
        elif isinstance(value, str):
            value = np.fromstring(value, dtype=self.dtype)

        value = super().validate(value, deps)

        if self.non_zero and (
            (value.ndims == 2 and not np.isclose(np.linalg.det(value), 0.0))
            or (value.ndims == 1 and not np.isclose(np.linalg.norm(value), 0.0))
        ):
            raise ConfigError("Non-zero Array has zero critical value.")

        return value


class Vector(Array):
    def __init__(
        self,
        *,
        length: int = 3,
        shape: None = None,
        normalise: bool = False,
        **params,
    ):
        if shape is not None:
            raise ConfigError("Cannot pass explicit `shape`, use `length`.")

        super().__init__(shape=(length,), **params)
        self.normalise = normalise

    def validate(self, value, deps: dict[str, Any]) -> np.ndarray:
        value = super().validate(value, deps)

        if self.normalise:
            value /= np.linalg.norm(value)

        return value
