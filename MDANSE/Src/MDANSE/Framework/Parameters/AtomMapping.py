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
from collections.abc import Collection, Sequence
from enum import Enum, auto
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from MDANSE.Chemistry import ATOMS_DATABASE
from MDANSE.Framework.AtomMapping import check_mapping_valid, fill_remaining_labels
from MDANSE.Framework.AtomSelector.selector import ReusableSelection
from MDANSE.IO.IOUtils import UCEnum, json_handler
from MDANSE.MLogging import LOG

from .Choices import SingleChoice
from .Parameters import ConfigError, ConfigureDescriptor
from .UtilTypes import Depends, DescID

if TYPE_CHECKING:
    from MDANSE.MolecularDynamics.Trajectory import Trajectory


class Mapper(ABC):
    def __init__(self):
        self._original_map = {}
        self._new_map = {}

    def get_full_setting(self) -> dict[int, float]:
        """
        Returns
        -------
        dict[int, float]
            The full partial charge setting.
        """
        return self._original_map | self._new_map

    def get_setting(self) -> dict[int, float]:
        """
        Returns
        -------
        dict[int, float]
            The minimal partial charge setting.
        """
        return {
            key: self._new_map[key] for key in self._original_map.keys() & self._new_map
        }

    def get_json_setting(self) -> str:
        """
        Returns
        -------
        str
            A json string of the minimal transmutation setting.
        """
        return json.dumps(self.get_setting())

    def reset_setting(self) -> None:
        """Resets the transmutation setting."""
        self._new_map = {}
        self.selector.reset()

    @staticmethod
    def get_selector(selection_string: str | Path | dict) -> ReusableSelection:
        selector = ReusableSelection()
        try:
            selection = json_handler(selection_string)
        except Exception as err:
            raise ConfigError("Failed to get selection dict.") from err
        selector.load_from_dict(selection)

        return selector

    def reset(self) -> None:
        """Resets the setting."""
        self._new_map = {}

    @abstractmethod
    def apply(self, selection_string: Path | str | dict, *args, **kwargs): ...


class AtomMapping(ConfigureDescriptor[dict | str, dict[str, dict[str, str]]]):
    default_tooltip = "Mapping of index to new species"

    def __init__(self, default: dict[str, dict[str, str]] | str = {}, **kwargs):
        super().__init__(default=default, **kwargs)

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def validate(
        self, value: dict | str, deps: Depends, /
    ) -> dict[str, dict[str, str]]:
        if not value or value == "{}":
            return {}

        file_info = deps["trajectory"]
        try:
            value = json_handler(value)
        except Exception as err:
            raise ConfigError("Failed to get replacement values.") from err

        labels = file_info.labels

        try:
            fill_remaining_labels(value, labels)
        except AttributeError as err:
            raise ConfigError("Unable to map all atoms.") from err

        if not check_mapping_valid(value, labels):
            raise ConfigError("Atom mapping not valid.")

        return value


class AtomSelection(ConfigureDescriptor[Collection[int], list[int]]):
    """Selects atoms in trajectory based on the input string.

    This configurator allows the selection of a specific set of
    atoms on which the analysis will be performed. The defaults setting
    selects all atoms.

    Attributes
    ----------
    _default : str
        The defaults selection setting.

    """

    def __init__(self, default: Collection[int] | dict | None = None, **kwargs):
        super().__init__(default=default, **kwargs)

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def validate(
        self,
        value: Sequence[int] | dict | str | Path,
        deps: Depends,
        /,
    ) -> list[int]:
        trajectory = deps["trajectory"]

        if value is None or value == "all":
            return list(range(len(trajectory.atom_types)))
        if isinstance(value, list | tuple | set):
            return sorted(value)

        selector = ReusableSelection()
        try:
            selection = json_handler(value)
        except Exception as err:
            raise ConfigError("Failed to get selection dict.") from err

        try:
            selector.load_from_dict(selection)
        except TypeError as err:
            raise ConfigError("Unable to load from dictionary.") from err

        indices = selector.select_in_trajectory(trajectory)

        return sorted(indices)


class PartialCharge(ConfigureDescriptor[str | Path | dict, dict[int, float]]):
    """The partial charge mapper. Updates an atom partial charge map
    with applications of the update_charges method with a selection
    setting and partial charge."""

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def validate(
        self,
        value: str | Path | dict,
        deps: Depends,
        /,
    ) -> dict[int, float]:
        try:
            value = json_handler(value)
        except Exception as err:
            raise ConfigError("Failed to get partial charge dict.") from err

        processed_values = {}
        for key, val in value.items():
            if isinstance(val, Collection):
                try:
                    charge = float(key)
                except (TypeError, ValueError) as err:
                    raise ConfigError(
                        f"Wrong charge {key} in the charge dictionary"
                    ) from err

                processed_values.update(dict.fromkeys(map(int, val), charge))
                continue

            try:
                index, charge = int(key), float(val)
            except (TypeError, ValueError) as err:
                raise ConfigError(
                    f"Index/charge pair {key}/{val} is not a valid int/float pair."
                ) from err
            processed_values[index] = charge

        system = deps["trajectory"].chemical_system

        if not processed_values.keys() <= set(system._atom_indices):
            LOG.warning("At least one atom index not found in the current system.")

        return {
            ind: processed_values[ind]
            for ind in processed_values.keys() & system._atom_indices
        }


class PartialChargeMapper(Mapper):
    """The partial charge mapper. Updates an atom partial charge map
    with applications of the update_charges method with a selection
    setting and partial charge."""

    def __init__(self, trajectory: Trajectory):
        """
        Parameters
        ----------
        system : ChemicalSystem
            The chemical system object.
        """
        super().__init__()
        system = trajectory.chemical_system
        self._traj_charges = trajectory.charges(0)[:]
        self._current_trajectory = trajectory
        self._original_map = {
            at_num: self._traj_charges.get(at_num, 0)
            for at_num, at in enumerate(system.atom_list)
        }

    def apply(self, selection_string: Path | str | dict, charge: float) -> None:
        """With the selection dictionary update the selector and then
        update the partial charge map.

        Parameters
        ----------
        selection_dict : dict[str, Union[bool, dict]]
            The selection setting to get the indices to map the input
            partial charge.
        charge : float
            The partial charge to map the selected atoms to.
        """
        selector = self.get_selector(selection_string)
        indices = selector.select_in_trajectory(self._current_trajectory)
        self._new_map.update(dict.fromkeys(indices, charge))

    def get_grouped_setting(self) -> dict[tuple[int], float]:
        """Return a dict of groups of indices with the same charge.

        Returns
        -------
        dict[tuple[int], float]
            The minimal partial charge setting.
        """
        new_charges = self._traj_charges.copy() | self._new_map
        valid_indices = np.where(~np.isclose(new_charges, self._traj_charges))
        unique_charges = np.unique(new_charges[valid_indices])

        groups = {}
        for charge in unique_charges:
            charge_indices = set(np.where(np.isclose(new_charges, charge))[0])
            key = charge_indices.intersection(valid_indices[0])
            groups[charge] = [int(x) for x in key]
        return groups

    def get_json_setting(self) -> str:
        """
        Returns
        -------
        str
            A json string of the minimal partial charge setting.
        """
        return json.dumps(self.get_grouped_setting())


class GroupingLevels(UCEnum):
    ATOM = auto()
    MOLECULE = auto()

    def __repr__(self) -> str:
        return self.name.title()


class GroupingLevel(SingleChoice[GroupingLevels | str, GroupingLevels]):
    base_aliases = {
        "each atom": GroupingLevels.ATOM,
        "each molecule": GroupingLevels.MOLECULE,
    }

    def __init__(
        self,
        *args,
        choices: None = None,
        default: GroupingLevels | str = GroupingLevels.ATOM,
        **kwargs,
    ):
        if choices is not None:
            raise ConfigError(f"Cannot provide choices to {type(self).__name__}")

        aliases = kwargs.pop("aliases", {}) | self.base_aliases

        super().__init__(
            *args,
            choices=GroupingLevels,
            default=default,
            aliases=aliases,
            **kwargs,
        )

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def validate(
        self,
        value: GroupingLevels | str,
        deps: Depends,
        /,
    ) -> str:
        value = super().validate(value, deps)

        if (
            value is GroupingLevels.MOLECULE
            and not deps["trajectory"].chemical_system.unique_molecules()
        ):
            raise ConfigError("Requested molecule grouping, but no molecules exist")

        return value


class AtomTransmuter(Mapper):
    """The atom transmuter setting generator. Updates an atom
    transmutation setting with applications of the apply_transmutation
    method with a selection setting and symbol."""

    def __init__(self, trajectory: Trajectory) -> None:
        """
        Parameters
        ----------
        system : ChemicalSystem
            The chemical system object.
        """
        super().__init__()
        self._current_trajectory = trajectory
        self._original_map = {
            number: element
            for number, element in enumerate(trajectory.chemical_system.atom_list)
        }

    def apply(self, selection_string: str, symbol: str) -> None:
        """With the selection dictionary update selector and then
        update the transmutation map.

        Parameters
        ----------
        selection_string : str
            the JSON string of the selection operation to use.
        symbol : str
            The element to map the selected atoms to.
        """
        if symbol not in ATOMS_DATABASE:
            raise ValueError(f"{symbol} not found in the atom database.")

        selector = self.get_selector(selection_string)
        indices = selector.select_in_trajectory(self._current_trajectory)
        self._new_map.update(dict.fromkeys(indices, symbol))


class AtomTransmutation(ConfigureDescriptor[dict | str | Path, dict]):
    def __init__(self, default: dict[int, str] | None = None, **kwargs):
        if default is None:
            default = {}
        super().__init__(default=default, **kwargs)

    def validate(
        self,
        value: dict | str | Path,
        deps: Depends,
        /,
    ):
        try:
            value = json_handler(value)
        except Exception as err:
            raise ConfigError("Failed to get transmutation dict.") from err

        if not value:
            return value

        try:
            value = {int(idx): element for idx, element in value.items()}
        except ValueError:
            raise ConfigError("Keys of transmutation map should be castable to `int`")

        system = deps["trajectory"].chemical_system
        idxs = range(system.total_number_of_atoms)

        if missing := value.keys() - idxs:
            raise ConfigError(
                f"Input setting not valid - atom indices ({', '.join(missing)}) not found in the current system."
            )

        elements = set(value.values())
        known_atoms = set(deps["trajectory"].atoms) | set(ATOMS_DATABASE.atoms)

        if missing := elements - known_atoms:
            raise ConfigError(
                f"Elements ({', '.join(missing)}) not registered in the database"
            )

        return value
