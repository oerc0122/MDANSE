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
from collections.abc import Iterator, Sequence
from enum import Enum, auto
from functools import singledispatchmethod
from math import isclose
from numbers import Number
from pathlib import Path
from typing import Literal, SupportsFloat, SupportsInt, TypeVar, overload, Generic

import numpy as np
import numpy.typing as npt
from more_itertools import numeric_range

from MDANSE.IO.IOUtils import json_handler

from .AbsConfigDesc import ConfigError, ConfigureDescriptor, MinMax
from .UtilTypes import Depends

K = TypeVar("K")
V = TypeVar("V")
T = TypeVar("T")

cjoin = ", ".join

__all__ = [
    "Boolean",
    "Float",
    "Integer",
    "String",
    "PathParam",
    "Range",
    "Array",
    "Vector",
]


def _nop(x: str) -> str:
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


class NumericRange(Generic[T]):
    """Wrapper for :class:`more_itertools.numeric_range` for non-private start stop."""

    def __init__(self, *params: T):
        self.iterator = numeric_range(*params)

    def __iter__(self) -> Iterator[T]:
        return iter(self.iterator)

    def __getitem__(self, index: int) -> T:
        return self.iterator[index]

    def __len__(self) -> int:
        return len(self.iterator)

    @property
    def start(self) -> T:
        return self.iterator._start

    @property
    def stop(self) -> T:
        return self.iterator._stop

    @property
    def step(self) -> T:
        return self.iterator._step

    def __str__(self) -> str:
        return str(self.iterator)

    def __repr__(self) -> str:
        return repr(self.iterator)

    def __int__(self) -> NumericRange[int]:
        return NumericRange(int(self.start), int(self.stop), int(self.step))


Bools = Literal["true", "false", "yes", "no", "y", "n", "1", "0", 1, 0, True, False]


class Boolean(ConfigureDescriptor[Bools, bool]):
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

    def __init__(
        self,
        default: Bools = False,
        **params,
    ):
        super().__init__(default=default, **params)

    def validate(self, value: Bools, deps: Depends, /) -> bool:
        test = value.lower() if isinstance(value, str) else value

        if test not in self._alias:
            raise ConfigError(f"Unable to convert value ({value!r}) to bool.")

        return self._alias[test]


class Float(MinMax[float], ConfigureDescriptor[SupportsFloat, float]):
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

    def validate(self, value: SupportsFloat, deps: Depends, /) -> float:
        try:
            value = float(value)
        except ValueError as error:
            raise ConfigError(f"Value ({value}) is not a valid float.") from error

        super().validate(value, deps)

        ranges = self.get_ranges(deps)

        self.validate_range(value, ranges)

        if self.non_zero is not None and isclose(value, 0.0, abs_tol=self.non_zero):
            raise ConfigError(f"Non-null val ({value}) is too close to zero.")

        return value


class Integer(MinMax[int], ConfigureDescriptor[SupportsInt, int]):
    """Configurator takes an integer input."""

    def __init__(self, **params):
        super().__init__(**params)

    def validate(self, value: SupportsInt, deps: Depends, /) -> int:
        try:
            value = int(value)
        except ValueError as error:
            raise ConfigError(f"Value ({value}) is not a valid integer.") from error

        ranges = self.get_ranges(deps)

        self.validate_range(value, ranges)

        return value


class Dict(
    ConfigureDescriptor[Sequence[tuple[K, V]] | dict[K, V] | str | Path, dict[K, V]]
):
    """
    This configurator allows a user to input a dict value.

    `choices` are keys which *must* be present in the dict.

    `exclude` are keys which *must not* be present in the dict.
    """

    def __init__(self, *, fixed: bool = False, choices: Sequence[K] = (), **params):
        super().__init__(choices=choices, **params)

        self.fixed = fixed

    def _validate_choices(
        self, value: dict[K, V], choices: set[K] | None = None
    ) -> bool:
        choice = self.choices if choices is None else choices
        if missing := choice - value.keys():
            raise ConfigError(
                f"Required keys ({cjoin(map(str, missing))}) are missing."
            )
        return True

    def _validate_exclude(
        self, value: dict[K, V], exclude: set[K] | None = None
    ) -> bool:
        excl = self.exclude if exclude is None else exclude
        if erroneous := excl & value.keys():
            raise ConfigError(
                f"Forbidden keys ({cjoin(map(str, erroneous))}) are present."
            )
        return True

    def validate(
        self,
        value: Sequence[tuple[K, V]] | dict[K, V] | str | Path,
        deps: Depends, /,
    ) -> dict[K, V]:
        if not isinstance(value, (dict, str, Path)):
            value = dict(value)
        else:
            try:
                value = json_handler(value)
            except Exception as err:
                raise ConfigError("Failed to build dict.") from err

        value = super().validate(value, deps)

        if self.fixed and (extraneous := value.keys() - self.choices):
            raise ConfigError(
                f"Extra keys are present ({cjoin(map(str, extraneous))}), "
                f"expected: {cjoin(map(str, self.choices))}."
            )

        return value


class String(ConfigureDescriptor[str, str]):
    """
    This configurator allows to input a string value.
    """

    def __init__(
        self,
        *,
        case: Callable[[str], str] = StrCases.PRESERVE,
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
    def regex(self, value: str | re.Pattern | None) -> None:
        if value is None:
            self._regex = None
            return
        self.regex = re.compile(value)

    def validate(self, value: str, deps: Depends, /) -> str:
        try:
            value = str(value)
        except ValueError as error:
            raise ConfigError(f"Value ({value}) is not a valid str.") from error

        value = self.case.value(value)

        super().validate(value, deps)

        if self.regex is not None and not self.regex.match(value):
            raise ConfigError(
                f"Value ({value}) does not match requirement {self.regex.pattern}."
            )

        return value


class PathParam(ConfigureDescriptor[str | Path, Path]):
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

    def validate(self, value, deps: Depends, /) -> Path:
        try:
            value = Path(value).expanduser()
        except TypeError as error:
            raise ConfigError(f"Value ({value}) is not a valid Path.") from error

        super().validate(value, deps)

        if self.mode is self.FileModes.MAY_EXIST:
            pass
        elif self.mode is self.FileModes.MUST_EXIST and not value.exists():
            raise ConfigError(f"File at ({value}) does not exist.")
        elif self.mode is self.FileModes.MUST_NOT_EXIST and value.exists():
            raise ConfigError(f"File at ({value}) must not exist.")

        return value


T = TypeVar("T", int, float)


class Range(
    MinMax[T],
    ConfigureDescriptor[
        int
        | float
        | list[T]
        | tuple[T, ...]
        | range
        | numeric_range
        | NumericRange[T],
        NumericRange[T],
    ],
):
    """
    This configurator allows a user to input a range value.
    """

    def __init__(
        self,
        *args,
        include_last: bool = False,
        dtype: type = float,
        **params,
    ):
        super().__init__(*args, **params)
        self.dtype = dtype
        self.include_last = include_last

    @singledispatchmethod
    def validate(self, value, deps) -> NumericRange:
        raise ConfigError(
            f"Do not know how to convert {type(value).__name__} to NumericRange"
        )

    @validate.register(NumericRange)
    def _(self, value: NumericRange[T], deps: Depends, /) -> NumericRange[T]:
        if self.include_last:
            value = NumericRange(value.start, value.stop + value.step, value.step)

        ranges = self.get_ranges(deps)

        self.validate_range(value.start, ranges)
        self.validate_range(value.stop, ranges)

        return value

    @validate.register
    def _(self, value: range, deps: Depends, /) -> NumericRange[int]:
        cls = NumericRange(value.start, value.stop, value.step)
        return self.validate(cls, deps)

    @validate.register
    def _(self, value: numeric_range, deps: Depends, /) -> NumericRange:
        cls = NumericRange(value._start, value._stop, value._step)
        return self.validate(cls, deps)

    @validate.register(int)
    @validate.register(float)
    def _(self, value: T, deps: Depends, /) -> NumericRange[T]:
        cls = NumericRange(self.dtype(value))
        return self.validate(cls, deps)

    @validate.register(list)
    @validate.register(tuple)
    def _(self, value: list[T] | tuple[T, ...], deps: Depends, /) -> NumericRange[T]:
        cls = NumericRange(*map(self.dtype, value))
        return self.validate(cls, deps)


class Array(ConfigureDescriptor[str | npt.ArrayLike, np.ndarray]):
    def __init__(
        self,
        *,
        shape: tuple[int, ...],
        dtype: type = np.float64,
        non_zero: bool = False,
        **params,
    ):
        super().__init__(**params)

        self.shape = shape
        self.non_zero = non_zero
        self.dtype = dtype

    def validate(self, value: str | npt.ArrayLike, deps: Depends, /) -> np.ndarray:
        if isinstance(value, str):
            value = np.fromstring(value, dtype=self.dtype)
        elif isinstance(value, Sequence):
            value = np.array(value, dtype=self.dtype)

        value = super().validate(value, deps)

        assert isinstance(value, np.ndarray)

        if self.non_zero and (
            (value.ndim == 2 and not np.isclose(np.linalg.det(value), 0.0))
            or (value.ndim == 1 and not np.isclose(np.linalg.norm(value), 0.0))
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

    def validate(self, value, deps: Depends, /) -> np.ndarray:
        value = super().validate(value, deps)

        if self.normalise:
            value /= np.linalg.norm(value)

        return value
