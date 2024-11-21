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
    if check_exists:
        for atm in system.atom_list:
            if atm.element == "dummy":
                return True
        return False
    else:
        return set([at.index for at in system.atom_list if at.element == "dummy"])


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
        for atm in system.atom_list:
            if atm.name == name:
                return True
        return False
    else:
        return set([at.index for at in system.atom_list if at.name == name])


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
        for atm in system.atom_list:
            if atm.full_name == fullname:
                return True
        return False
    else:
        return set([at.index for at in system.atom_list if at.full_name == fullname])


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
