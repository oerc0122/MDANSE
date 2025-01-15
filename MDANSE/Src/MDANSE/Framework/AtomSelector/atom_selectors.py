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
from typing import Union
from MDANSE.MolecularDynamics.Trajectory import Trajectory


__all__ = [
    "select_element",
    "select_dummy",
    "select_atom_name",
    "select_atom_fullname",
    "select_hs_on_element",
    "select_hs_on_heteroatom",
    "select_index",
]


def select_element(
    trajectory: Trajectory, symbol: str, check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects all atoms for the input element.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    symbol : str
        Symbol of the element.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        The atom indices of the matched atoms.
    """
    system = trajectory.chemical_system
    pattern = f"[#{trajectory.get_atom_property(symbol, 'atomic_number')}]"
    if check_exists:
        return system.has_substructure_match(pattern)
    else:
        return system.get_substructure_matches(pattern)


def select_dummy(
    trajectory: Trajectory, check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects all dummy atoms in the chemical system.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        All dummy atom indices or a bool if checking match.
    """
    system = trajectory.chemical_system
    dummy_list = ["Du", "dummy"]
    if check_exists:
        for atm in system.atom_list:
            if atm in dummy_list:
                return True
            elif trajectory.get_atom_property(atm, "dummy"):
                return True
        return False
    else:
        for atm in system._unique_elements:
            if trajectory.get_atom_property(atm, "dummy"):
                dummy_list.append(atm)
        return set(
            [
                index
                for index, element in enumerate(system.atom_list)
                if element in dummy_list
            ]
        )


def select_atom_name(
    trajectory: Trajectory, name: str, check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects all atoms with the input name in the chemical system.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    name : str
        The name of the atom to match.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        All atom indices or a bool if checking match.
    """
    system = trajectory.chemical_system
    if check_exists:
        if name in system.atom_list:
            return True
        return False
    else:
        return set(
            [index for index, element in enumerate(system.atom_list) if element == name]
        )


def select_atom_fullname(
    trajectory: Trajectory, fullname: str, check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects all atoms with the input fullname in the chemical system.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    fullname : str
        The fullname of the atom to match.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        All atom indices or a bool if checking match.
    """
    system = trajectory.chemical_system
    if check_exists:
        if fullname in system.name_list:
            return True
        return False
    else:
        return set(
            [index for index, name in enumerate(system.name_list) if name == fullname]
        )


def select_hs_on_element(
    trajectory: Trajectory, symbol: str, check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects all H atoms bonded to the input element.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    symbol : str
        Symbol of the element that the H atoms are bonded to.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        The atom indices of the matched atoms.
    """
    system = trajectory.chemical_system
    num = trajectory.get_atom_property(symbol, "atomic_number")
    if check_exists:
        return system.has_substructure_match(f"[#{num}]~[H]")
    else:
        xh_matches = system.get_substructure_matches(f"[#{num}]~[H]")
        x_matches = system.get_substructure_matches(f"[#{num}]")
        return xh_matches - x_matches


def select_hs_on_heteroatom(
    trajectory: Trajectory, check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects all H atoms bonded to any atom except carbon and
    hydrogen.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        The atom indices of the matched atoms.
    """
    system = trajectory.chemical_system
    if check_exists:
        return system.has_substructure_match("[!#6&!#1]~[H]")
    else:
        xh_matches = system.get_substructure_matches("[!#6&!#1]~[H]")
        x_matches = system.get_substructure_matches("[!#6&!#1]")
        return xh_matches - x_matches


def select_index(
    trajectory: Trajectory, index: Union[int, str], check_exists: bool = False
) -> Union[set[int], bool]:
    """Selects atom with index - just returns the set with the
    index in it.

    Parameters
    ----------
    system : ChemicalSystem
        The MDANSE chemical system.
    index : int or str
        The index to select.
    check_exists : bool, optional
        Check if a match exists.

    Returns
    -------
    Union[set[int], bool]
        The index in a set or a bool if checking match.
    """
    if check_exists:
        return True
    else:
        return {int(index)}
