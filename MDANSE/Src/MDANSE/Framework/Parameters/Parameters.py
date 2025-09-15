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
from collections.abc import Callable, Iterable, Mapping, Sequence
from contextlib import suppress
from enum import Enum
from functools import singledispatchmethod
from itertools import compress
from numbers import Real
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    Protocol,
    TypeGuard,
    TypeVar,
    cast,
    overload,
    runtime_checkable,
)

import numpy as np
from more_itertools import first_true, value_chain

from MDANSE.Core.Error import Error
from MDANSE.IO.IOUtils import MDANSEEncoder
from MDANSE.MLogging import LOG

from .UtilTypes import Depends, DescID

if TYPE_CHECKING:
    from typing import Self

SENTINEL = object()
P = TypeVar("P")
T = TypeVar("T")
Num = TypeVar("Num", bound=Real)
K = TypeVar("K")
V = TypeVar("V")
CD = TypeVar("CD", bound="ConfigureDescriptor")

cjoin = ", ".join


@runtime_checkable
class Validatable(Protocol):
    def validate(self, desc: ConfigureDescriptor, value: T) -> T: ...


@runtime_checkable
class CustomChoices(Protocol):
    def get_choices(self, deps: Depends) -> set: ...


@runtime_checkable
class HasOption(Protocol):
    @property
    def choices(self) -> set: ...

    @property
    def exclude(self) -> set: ...


def is_enum(x: Any) -> TypeGuard[type[Enum]]:
    return isinstance(x, type) and issubclass(x, Enum)


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

    @classmethod
    def _get_default_parameters(cls):
        return {
            key: desc.default
            if hasattr(desc, "default") and desc.default is not SENTINEL
            else "N/A"
            for key, desc in cls._get_descriptors().items()
        }

    default_parameters = property(lambda self: type(self)._get_default_parameters())

    @property
    def configuration(self) -> dict[DescID, Any]:
        return {name: getattr(self, name) for name in self.descriptors}

    @configuration.setter
    def configuration(self, value: dict[DescID, Any]):
        if extra := value.keys() - self.parameters:
            raise ConfigError(
                f"Unrecognised parameters in dict for {type(self).__name__}: {cjoin(extra)}."
            )

        # Preserve definition order
        for name in self.parameters:
            if name in value:
                setattr(self, name, value[name])

    @classmethod
    def _get_descriptors(cls) -> dict[DescID, Parameter]:
        return {
            DescID(name): param
            for par in value_chain(cls.__bases__, cls)
            for name, param in par.__dict__.items()
            if isinstance(param, Parameter)
        }

    descriptors = property(lambda self: type(self)._get_descriptors())

    @property
    def parameters(self) -> list[DescID]:
        return list(self.descriptors)

    def to_json(self) -> str:
        to_json = {
            name: obj.to_json() if isinstance(obj, Parameter) else obj
            for name, obj in self.configuration.items()
        } | {"class_type": type(self).__name__}
        return json.dumps(to_json, cls=MDANSEEncoder)

    output_configuration = to_json

    def check_status(self) -> bool:
        try:
            for param in self.parameters:
                getattr(self, param)
        except ConfigError:
            return False

        return True

    @classmethod
    def build_doc_example(cls) -> str:
        docstring = ":Example:\n\n"
        docstring += ">>> \n"
        docstring += ">>> \n"
        docstring += ">>> parameters = {}\n"
        for k, v in cls._get_default_parameters().items():
            docstring += f">>> parameters[{k!r}]={v!r}\n"
        docstring += ">>> \n"
        docstring += f">>> job = IJob.create({cls.__name__!r})\n"
        docstring += ">>> job.setup(parameters)\n"
        docstring += ">>> job.run()\n"
        return docstring

    @classmethod
    def build_doc_texttable(cls, doclist: list[dict[str, str]]) -> str:
        docstring = "\n**Job input configurators:** \n\n"

        columns = ["Parameter", "Default value", "Description"]

        sizes = [len(v) for v in columns]

        for v in doclist:
            sizes[0] = max(sizes[0], len(v["Parameter"]))
            sizes[1] = max(sizes[1], len(v["Default value"]))
            # Case of Description field: has to be splitted and parsed for inserting sphinx "|" keyword for multiline
            v["Description"] = v["Description"].strip()
            v["Description"] = v["Description"].splitlines()
            v["Description"] = ["| " + vv.strip() for vv in v["Description"]]
            sizes[2] = max(value_chain(sizes[2], map(len, v["Description"])))

        data_line = "| " + "| ".join(f"{{}}:<{size}" for size in sizes) + "|\n"
        sep_line = "+" + "+".join("-" * (size + 1) for size in sizes) + "+\n"

        docstring += sep_line
        docstring += data_line.format(*columns)
        docstring += sep_line.replace("-", "=")

        for v in doclist:
            docstring += data_line.format(
                v["Parameter"], v["Default value"], v["Description"][0]
            )
            if len(v["Description"]) > 1:
                for descr in v["Description"][1:]:
                    data_line.format("", "", descr)
            docstring += sep_line

        docstring += "\n"
        return docstring

    @classmethod
    def build_doc_htmltable(cls, doclist: list[dict[str, str]]) -> str:
        docstring = "\n**Job input configurators:**"

        columns = ["Parameter", "Default value", "Description"]

        for v in doclist:
            # Case of Description field: has to be splitted and parsed for inserting sphinx "|" keyword for multiline
            v["Description"] = v["Description"].strip()
            v["Description"] = v["Description"].split("\n")
            v["Description"] = ["" + vv.strip() for vv in v["Description"]]

        docstring += "<table>\n"
        docstring += "<tr>"
        for col in columns:
            docstring += f"<th>{col}</th>"
        docstring += "</tr>\n"

        for v in doclist:
            docstring += "<tr>"
            for item in [
                v["Parameter"],
                v["Default value"],
                v["Description"][0],
            ]:
                docstring += f"<td>{item}</td>"
            docstring += "</tr>\n"

        docstring += "</table>\n"
        return docstring

    @classmethod
    def build_doc(cls, use_html_table: bool = False) -> str:
        """Return the documentation about a configurable class based on its configurators contents.

        Parameters
        ----------
        cls : MDANSE.Framework.Parameters.Configurable
            The configurable class for which documentation should be built.
        use_html_table : bool
             Use an HTML table.

        Returns
        -------
        str
            The documentation about the configurable class.

        """

        doclist = []

        defaults = cls._get_default_parameters()

        for name, desc in cls._get_descriptors().items():
            descr = getattr(desc, "description", "")
            descr += "\n" + str(desc.__doc__)
            doclist.append(
                {
                    "Parameter": name,
                    "Default value": repr(defaults[name]),
                    "Description": descr,
                }
            )

        docstring = cls.build_doc_example()

        if use_html_table:
            docstring += cls.build_doc_htmltable(doclist)
        else:
            docstring += cls.build_doc_texttable(doclist)

        return docstring

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
    """Mixin for classes which contain Descriptors, but are to be configured."""

    def __set_name__(self, owner: object, name: str):
        self.name = name
        self._owner = owner

        # Set descriptors to true owner and remove config from self.
        for desc_name, desc in self.descriptors.items():
            with suppress(AttributeError):
                delattr(type(self), desc.configured_var)
            desc.__set_name__(owner, f"{self.name}_{desc_name}")

    def __get__(self, owner: object, objtype: type | None = None):
        # Make descriptor for unique per-instance use.
        return self

    def __getattribute__(self, name: str):
        # If dealing with a descriptor call the get using config's parent
        # rather than descriptor's.
        desc = super().__getattribute__("descriptors")
        if name in desc:
            owner = self._owner
            return desc[name].__get__(owner, type(owner))
        # Fallback for config properties
        return super().__getattribute__(name)

    def __setattr__(self, name: str, value: Any):
        # If dealing with a descriptor call the get using config's parent
        # rather than descriptor's.
        if name in self.descriptors:
            self.descriptors[name].__set__(self._owner, value)
            return
        super().__setattr__(name, value)


def to_class(cls):
    """Simplified wrapper for callbacks."""

    def tc(*args):
        return cls(args[1])

    return tc


class ConfigureDescriptor(Parameter, Generic[P, T]):
    """Abstract configure descriptor.

    Parameters
    ----------
    default : P
        Default value if not configured.
    optional : bool
        Whether this value can be disabled (None).
    choices : Sequence[T], optional
        Valid options for this value.
    exclude : Sequence[T]
        Invalid options for this value.
    mutex : Sequence[str]
        Other variables with which this is mutually exclusive.
    depends : Sequence[str]
        Other variables which must be set for this to be set.
    on_set : Callable
        Custom function to be called post-validation.
    on_get : Callable
        Custom function to be called before returning value.
    label : str
        GUI label to present to user.
    tooltip : str
        GUI tooltip to present to user.
    """

    def __init__(
        self,
        *,
        default: P | object = SENTINEL,
        optional: bool | None = None,
        choices: Sequence[T] | None = None,
        exclude: Sequence[T] = (),
        mutex: Sequence[DescID] = (),
        depends: dict[str, DescID] | None = None,
        on_set: Callable[[Self, T, Depends], T] | None = None,
        on_get_depends: dict[str, DescID] | None = None,
        on_get: Callable[[Self, T, Depends], T] | None = None,
        label: str = "",
        tooltip: str = "",
        **params,
    ):
        self.default: P = default
        self.optional = optional if optional is not None else default is not SENTINEL

        self.choices = choices
        self.exclude = exclude

        self.mutex = mutex

        self.depends: dict[str, DescID] = depends if depends is not None else {}
        self.get_depends: dict[str, DescID] = (
            on_get_depends if on_get_depends is not None else {}
        )

        self.on_set = on_set
        self.on_get = on_get

        self.dependents: set[ConfigureDescriptor] = set()

        if missing := (self.required_deps() - self.depends.keys()):
            raise ConfigError(
                f"Required deps ({cjoin(missing)}) missing for {type(self).__name__}."
            )

        super().__init__(label=label, tooltip=tooltip)

    def _find_dep_class(self, typ: type) -> ConfigureDescriptor:
        return first_true(
            self.dependents,
            pred=lambda x: isinstance(x, typ),
        )

    def _bad_mutex(self, owner: object) -> Iterable[bool]:
        return (getattr(owner, f"_{ex}_configured") for ex in self.mutex)

    def _bad_deps(
        self, owner: object, depends: Depends | None = None
    ) -> Iterable[bool]:
        depends = depends if depends is not None else self.depends
        return (not getattr(owner, f"_{dep}_configured") for dep in depends.values())

    def _get_deps(self, owner: object, depends: Depends | None = None) -> Depends:
        depends = depends if depends is not None else self.depends
        return cast(Depends, {dep: getattr(owner, key) for dep, key in depends.items()})

    def __set_name__(self, owner: object, name: str):
        self.name = name
        self.private_name = "_" + name
        self.configured_var = self.private_name + "_configured"

        # Update the depend checks
        for dep in self.depends.values():
            owner.__dict__[dep].dependents.add(self)

        setattr(owner, self.configured_var, self.optional)

    def __get__(self, owner: object, objtype: type | None = None) -> Self | T | None:
        if owner is None:
            return self

        if any(self._bad_deps(owner)):
            raise ConfigError(
                f"Dependencies ({cjoin(compress(self.depends, self._bad_deps(owner)))}) "
                "are not correctly defined."
            )

        if not self.optional and not getattr(owner, self.configured_var):
            raise ConfigError(f"Non-optional value ({self.name}) has not been set")

        value = getattr(owner, self.private_name, SENTINEL)

        if self.optional and value is SENTINEL:
            deps = self._get_deps(owner)
            value = self.validate(self.default, deps)

        out = cast(T, value)

        # Custom
        if self.on_get is not None:
            cb_deps = self._get_deps(owner, self.get_depends)
            out = self.on_get(self, out, cb_deps)

        return out

    def __set__(self, owner: object, value: P) -> None:
        setattr(owner, self.configured_var, False)
        setattr(owner, self.private_name, SENTINEL)

        if self.optional and value is None:
            return

        if any(self._bad_deps(owner)):
            raise ConfigError(
                f"Dependencies for {self.name} ({cjoin(compress(self.depends, self._bad_deps(owner)))}) "
                "are not correctly defined."
            )
        if any(self._bad_mutex(owner)):
            raise ConfigError(
                f"Mutually exclusive value ({cjoin(compress(self.mutex, self._bad_mutex(owner)))}) "
                "is also configured."
            )

        deps = self._get_deps(owner)
        out = self.validate(value, deps)

        # For grouped config descriptors
        if isinstance(owner, Validatable):
            out = owner.validate(self, out)

        # Custom
        if self.on_set is not None:
            out = self.on_set(self, out, deps)

        setattr(owner, self.private_name, out)
        setattr(owner, self.configured_var, True)

        for dependent in self.dependents:
            name = dependent.name
            try:
                val = getattr(owner, name)
                setattr(owner, name, val)
            except ConfigError as err:
                if err._msg and "Non-optional" not in err._msg:
                    LOG.info(f"Invalidated {name!r}. Reason: {err}")
                    setattr(owner, dependent.configured_var, False)

    def required_deps(self) -> set[DescID]:
        return set()

    @property
    def choices(self) -> set[T] | type[Enum]:
        """
        Returns the set of values allowed for an input.
        """
        return self._choices

    @choices.setter
    def choices(self, value: Sequence[T] | type[Enum] | None) -> None:
        self._choices: set[T] | type[Enum]
        if value is None:
            self._choices = set()
        elif is_enum(value):
            self._choices = value
        else:
            value = cast(Sequence[T], value)
            self._choices = set(value)

        self.last_choices = self._choices

    @property
    def exclude(self) -> set[T]:
        """
        Returns the set of values which are not forbidden.

        Returns
        -------
        set[int]
            Forbidden values.
        """
        return self._exclude

    @exclude.setter
    def exclude(self, value: Sequence[T]) -> None:
        if value is None:
            self._exclude = set()
            return
        self._exclude = set(value)

    @overload  # Standard choice
    def _validate_choices(self, value: T, choices: set[T] | None = None) -> bool: ...
    @overload  # Require keys in mapping
    def _validate_choices(
        self, value: Mapping[K, V], choices: set[K] | None = None
    ) -> bool: ...
    @overload  # Enum choices
    def _validate_choices(
        self, value: Enum, choices: type[Enum] | None = None
    ) -> bool: ...
    @overload  # Multi-choice
    def _validate_choices(
        self, value: Sequence[T], choices: set[T] | None = None
    ) -> bool: ...

    def _validate_choices(self, value, choices=None) -> bool:
        test_choices = self.choices if choices is None else choices
        return not test_choices or value in test_choices

    @singledispatchmethod
    def _validate_exclude(self, value: T, exclude: set[T] | None = None) -> bool:
        exclude = self.exclude if exclude is None else exclude
        return not exclude or value not in exclude

    @_validate_exclude.register(Sequence)
    def _(self, value: Sequence[T], exclude: set[T] | None = None) -> bool:
        excludes = self.exclude if exclude is None else exclude
        return set(value) > excludes

    @_validate_exclude.register(Mapping)
    def _(self, value: Mapping[K, V], exclude: set[K] | None = None) -> bool:
        excl = self.exclude if exclude is None else exclude
        if erroneous := excl & value.keys():
            raise ConfigError(
                f"Forbidden keys ({cjoin(map(str, erroneous))}) are present."
            )
        return True

    @abstractmethod
    def validate(self, value: P, deps: Depends, /) -> T:
        """
        Ensure that the passed variable is of the right type.
        """
        # Assume that the variable is the right type at this point
        out = cast(T, value)

        # Choices
        if isinstance(self, CustomChoices):  # Function
            self.last_choices = self.get_choices(deps)

            if not self._validate_choices(out, self.last_choices):
                raise ConfigError(
                    f"Value ({out!r}) not in choices ({cjoin(self.last_choices)})."
                )

        elif self.choices and is_enum(self.choices):  # Enum
            try:
                out = self.choices(out)
            except ValueError:
                raise ConfigError(
                    f"Value ({out!r}) not in choices ({cjoin(choice.name for choice in self.choices)})."
                )

        elif not self._validate_choices(out):  # Set
            raise ConfigError(
                f"Value ({out!r}) not in choices ({cjoin(map(str, self.choices))})."
            )

        # Exclude
        if not self._validate_exclude(out):
            raise ConfigError(
                f"Value ({out!r}) in excluded values ({cjoin(map(str, self.exclude))})."
            )

        return out


class MinMax(ABC, Generic[Num]):
    def __init__(
        self, *args, minimum: Num | None = None, maximum: Num | None = None, **kwargs
    ):
        self.minimum = minimum
        self.maximum = maximum
        super().__init__(*args, **kwargs)

    def get_ranges(self, deps: Depends) -> tuple[Num | None, Num | None]:
        return self.minimum, self.maximum

    def validate_range(
        self, value: Num, ranges: tuple[Num | None, Num | None] | None = None
    ) -> None:
        if ranges is not None:
            mini, maxi = ranges
        else:
            mini, maxi = self.minimum, self.maximum

        try:
            self.validate_minimum(value, mini)
            self.validate_maximum(value, maxi)
        except ConfigError:
            raise ConfigError(
                f"Value ({value}) outside of valid range ({mini}, {maxi})"
            )

    def validate_minimum(self, value: Num, mini: Num | None = None) -> None:
        minimum = mini if mini is not None else self.minimum

        if minimum is None:
            return

        if value < minimum:
            raise ConfigError(f"Value ({value}) less than minimum ({minimum})")

    def validate_maximum(self, value: Num, maxi: Num | None = None) -> None:
        maximum = maxi if maxi is not None else self.maximum

        if maximum is None:
            return

        if value > maximum:
            raise ConfigError(f"Value ({value}) greater than maximum ({maximum})")
