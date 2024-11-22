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

import abc
from ast import literal_eval
from typing import Union, TYPE_CHECKING, List, Tuple, Dict
import h5py
import numpy as np
from rdkit import Chem
from numpy.typing import NDArray
from MDANSE.Chemistry import (
    ATOMS_DATABASE,
)
from MDANSE.Mathematics.Geometry import superposition_fit, center_of_mass
from MDANSE.Mathematics.LinearAlgebra import delta, Quaternion, Tensor, Vector
from MDANSE.Mathematics.Transformation import Rotation, RotationTranslation, Translation
from MDANSE.MLogging import LOG

if TYPE_CHECKING:
    from MDANSE.MolecularDynamics.Configuration import _Configuration


class UnknownAtomError(Exception):
    pass


class InvalidVariantError(Exception):
    pass


class InvalidChemicalEntityError(Exception):
    pass


class InconsistentChemicalSystemError(Exception):
    pass


class ChemicalEntityError(Exception):
    pass


class CorruptedFileError(Exception):
    pass


class _ChemicalEntity(metaclass=abc.ABCMeta):
    """Abstract base class for other chemical entities."""

    def __init__(self):
        self._parent = None

        self._name = ""

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state

    @property
    @abc.abstractmethod
    def atom_list(self):
        pass

    @abc.abstractmethod
    def copy(self):
        pass

    @property
    def full_name(self) -> str:
        """The full name of this chemical entity, which includes the names of all parent chemical entities."""
        full_name = self.name
        parent = self._parent
        while parent is not None:
            full_name = "{}.{}".format(parent.name, full_name)
            parent = parent.parent

        return full_name

    @property
    @abc.abstractmethod
    def number_of_atoms(self):
        pass

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, parent):
        self._parent = parent

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @abc.abstractmethod
    def serialize(self, h5_file):
        pass

    def group(self, name: str) -> list["Atom"]:
        """
        Finds all atoms in this chemical entity that are a part of the provided group.

        :param name: The name of the group whose atoms are being searched for.
        :type name: str

        :return: List of atoms that are a part of the provided group.
        :rtype: list
        """
        selected_atoms = []
        for at in self.atom_list:
            if not hasattr(at, "groups"):
                continue

            if name in at.groups:
                selected_atoms.append(at)

        return selected_atoms

    def center_of_mass(self, configuration: _Configuration) -> NDArray[np.float64]:
        """
        Determines the coordinates of the centre of mass of this chemical entity.

        :param configuration: The configuration corresponding to the chemical system whose part this chemical enitity is
        :type configuration: any subclass of MDANSE.MolecularDynamics.Configuration._Configuration

        :return: The coordinates of the centre of mass of this chemical entity.
        :rtype: numpy.ndarray
        """
        coords = configuration["coordinates"]
        masses = [
            ATOMS_DATABASE.get_atom_property(at.symbol, "atomic_weight")
            for at in self.atom_list
        ]
        indices = [at.index for at in self.atom_list]

        return center_of_mass(coords[indices], masses)

    def centre_of_mass(self, configuration: _Configuration) -> NDArray[np.float64]:
        """Wrapper around the :py:meth: `center_of_mass()` method."""
        return self.center_of_mass(configuration)

    @property
    def mass(self) -> float:
        """The mass of this chemical entity. This is the sum of all non-ghost atoms in this chemical entity."""
        return sum(self.masses)

    @property
    def masses(self) -> list[float]:
        """A list of masses of all non-ghost atoms within this chemical entity."""
        return [
            ATOMS_DATABASE.get_atom_property(at.symbol, "atomic_weight")
            for at in self.atom_list
        ]

    def find_transformation_as_quaternion(
        self, conf1: _Configuration, conf2: Union[_Configuration, None] = None
    ) -> tuple[Quaternion, Vector, Vector, float]:
        """
        Finds a linear transformation that, when applied to this chemical entity with its coordinates defined in
        configuration conf1, minimizes the RMS distance to the conformation in conf2. Alternatively, if conf2 is None,
        a linear transformation from the current configuration of the ChemicalSystem to conf1 is returned.

        Unlike :py:method: `find_transformation()`, this method returns a Quaternion corresponding to the rotation.
        Under the hood, this method calls :func: `MDANSE.Mathematics.Geometry.superposition_fit()`.

        :param conf1: the configuration which is considered the initial configuration for the transformation if conf2 is
                      set, and the final configuration if it is not
        :type conf1: :class: `~MDANSE.MolecularDynamics.Configuration.Configuration`

        :param conf2: the configuration which is considered the target configuration of the transformation, or None
        :type conf2: :class: `~MDANSE.MolecularDynamics.Configuration.Configuration` or NoneType

        :return: The quaternion corresponding to the rotation, vectors of the centres of mass of both configurations and
            the minimum root-mean-square distance between the configurations.
        :rtype: tuple of (:class: `MDANSE.Mathematics.LinearAlgebra.Quaternion`,
                          :class: `MDANSE.Mathematics.LinearAlgebra.Vector`,
                          :class: `MDANSE.Mathematics.LinearAlgebra.Vector`,
                          float)
        """
        chemical_system = self.root_chemical_system
        if chemical_system is None:
            raise ChemicalEntityError(
                "Only chemical entities which are a part of a ChemicalSystem can be transformed."
                f' The provided chemical entity, "{str(self)}", has its top level parent '
                f"{str(self.top_level_chemical_entity)} which is not a part of any ChemicalSystem"
                f".\nFull chemical entity information: {repr(self)}"
            )

        if chemical_system.configuration.is_periodic:
            raise ValueError(
                "superposition in periodic configurations is not defined, therefore the configuration of "
                "the root chemical system of this chemical entity must not be periodic."
            )

        if conf1.chemical_system != chemical_system:
            raise ValueError(
                "conformations come from different chemical systems: the root chemical system of this "
                f'chemical entity is "{chemical_system.name}" but the chemical system registered with the '
                f'provided configuration (conf1) is "{conf1.chemical_system.name}".\nRoot chemical system:'
                f" {repr(chemical_system)}\nConfiguration chemical system: {repr(conf1.chemical_system)}"
            )

        if conf2 is None:
            conf2 = conf1
            conf1 = chemical_system.configuration
        else:
            if conf2.chemical_system != chemical_system:
                raise ValueError(
                    "conformations come from different chemical systems: the root chemical system of this"
                    f' chemical entity is "{chemical_system.name}" but the chemical system registered with'
                    f' the provided configuration (conf2) is "{conf2.chemical_system.name}".\nRoot '
                    f"chemical system: {repr(chemical_system)}\nConfiguration chemical system: "
                    f"{repr(conf2.chemical_system)}"
                )

        weights = chemical_system.masses

        return superposition_fit(
            [
                (
                    weights[a.index],
                    Vector(*conf1["coordinates"][a.index, :]),
                    Vector(*conf2["coordinates"][a.index, :]),
                )
                for a in self.atom_list
            ]
        )

    def find_transformation(
        self, conf1: _Configuration, conf2: _Configuration = None
    ) -> tuple[RotationTranslation, float]:
        """
        Finds a linear transformation that, when applied to this chemical entity with its coordinates defined in
        configuration conf1, minimizes the RMS distance to the conformation in conf2. Alternatively, if conf2 is None,
        a linear transformation from the current configuration of the ChemicalSystem to conf1 is returned.

        It uses the :method: `find_transformation_as_quaternion()` to do this, and then transforms its results to return
        an :class: `MDANSE.Mathematics.Transformation.RotationTranslation` object.

        :param conf1: the configuration in which this chemical entity is
        :type conf1: :class:`~MDANSE.MolecularDynamics.Configuration.Configuration`
        `
        :param conf1: the configuration which is considered the initial configuration for the transformation if conf2 is
                      set, and the final configuration if it is not
        :type conf1: :class: `~MDANSE.MolecularDynamics.Configuration.Configuration`

        :param conf2: the configuration which is considered the target configuration of the transformation, or None
        :type conf2: :class: `~MDANSE.MolecularDynamics.Configuration.Configuration` or NoneType

        :returns: the linear transformation corresponding to the transformation from conf1 to conf2, or from current
                  configuration to conf1, and the minimum root-mean-square distance
        :rtype: tuple of (:class: `MDANSE.Mathematics.Transformation.RotationTranslation`, float)
        """
        q, cm1, cm2, rms = self.find_transformation_as_quaternion(conf1, conf2)
        return Translation(cm2) * q.asRotation() * Translation(-cm1), rms

    def center_and_moment_of_inertia(
        self, configuration: _Configuration
    ) -> tuple[Vector, Tensor]:
        """
        Calculates the centre of masses and the inertia tensor of this chemical entity.

        :param configuration: a configuration that contains coordinates of this chemical entity
        :type configuration: :class:`~MDANSE.MolecularDynamics.Configuration.Configuration`

        :returns: the center of mass and the moment of inertia tensor in the given configuration
        :rtype: tuple of (:class: `MDANSE.Mathematics.LinearAlgebra.Vector`,
        :class: `MDANSE.Mathematics.LinearAlgebra.Tensor`)
        """

        m = 0.0
        mr = Vector(0.0, 0.0, 0.0)
        t = Tensor(3 * [3 * [0.0]])
        for atom in self.atom_list:
            ma = ATOMS_DATABASE.get_atom_property(atom.symbol, "atomic_weight")
            r = Vector(configuration["coordinates"][atom.index, :])
            m += ma
            mr += ma * r
            t += ma * r.dyadic_product(r)
        cm = mr / m
        t -= m * cm.dyadic_product(cm)
        t = t.trace() * delta - t
        return cm, t

    def centre_and_moment_of_inertia(
        self, configuration: _Configuration
    ) -> tuple[Vector, Tensor]:
        """Wrapper around :meth: `center_and_moment_of_inertia()`."""
        return self.center_and_moment_of_inertia(configuration)

    def normalizing_transformation(
        self, configuration: _Configuration, representation: str = None
    ) -> RotationTranslation:
        """
        Calculate a linear transformation that shifts the center of mass  of the object to the coordinate origin and
        makes its principal axes of inertia parallel to the three coordinate axes.

        :param configuration: a configuration that contains coordinates of this chemical entity
        :type configuration: :class:`~MDANSE.MolecularDynamics.Configuration.Configuration`

        :param representation: the specific representation for axis alignment:
          Ir    : x y z <--> b c a
          IIr   : x y z <--> c a b
          IIIr  : x y z <--> a b c
          Il    : x y z <--> c b a
          IIl   : x y z <--> a c b
          IIIl  : x y z <--> b a c
        :type representation: str

        :returns: the normalizing transformation
        :rtype: :class: `MDANSE.Mathematics.Transformation.RotationTranslation`
        """

        cm, inertia = self.center_and_moment_of_inertia(configuration)
        ev, diag = np.linalg.eig(inertia.array)
        diag = np.transpose(diag)
        if np.linalg.det(diag) < 0:
            diag[0] = -diag[0]

        if representation is not None:
            seq = np.argsort(ev)
            if representation == "Ir":
                seq = np.array([seq[1], seq[2], seq[0]])
            elif representation == "IIr":
                seq = np.array([seq[2], seq[0], seq[1]])
            elif representation == "IIIr":
                pass
            elif representation == "Il":
                seq = seq[2::-1]
            elif representation == "IIl":
                seq[1:3] = np.array([seq[2], seq[1]])
            elif representation == "IIIl":
                seq[0:2] = np.array([seq[1], seq[0]])
            else:
                raise ValueError(
                    f"invalid input for parameter repr: a value of {repr(representation)} was provided, "
                    'but only the following values are accepted: "Ir", "IIr", "IIIr", "Il", "IIl", "IIIl"'
                )
            diag = np.take(diag, seq)

        return Rotation(diag) * Translation(-cm)

    @property
    def root_chemical_system(self) -> Union["ChemicalSystem", None]:
        """
        The :class: `MDANSE.Chemistry.ChemicalEntity.ChemicalSystem` of which part this chemical entity is, or None if
        it is not a part of one.
        """
        if isinstance(self, ChemicalSystem):
            return self
        else:
            try:
                return self._parent.root_chemical_system
            except AttributeError:
                return None

    @property
    def top_level_chemical_entity(self) -> "_ChemicalEntity":
        """
        The ultimate non-system parent of this chemical entity, i.e. the parent of a parent etc. until an entity that
        is directly a child of a :class: `MDANSE.Chemistry.ChemicalEntity.ChemicalSystem`, wherein the child is returned.
        If this chemical entity is not a part of a ChemicalSystem, the entity whose parent is None is returned.
        """
        if isinstance(self._parent, ChemicalSystem):
            return self
        else:
            try:
                return self._parent.top_level_chemical_entity
            except AttributeError:
                return self

    @property
    @abc.abstractmethod
    def total_number_of_atoms(self):
        pass


class Atom(_ChemicalEntity):
    """A representation of atom in a trajectory."""

    def __init__(
        self,
        symbol: str = "H",
        name: str = None,
        bonds: list["Atom"] = None,
        groups: list[str] = None,
        ghost: bool = False,
        **kwargs,
    ):
        """
        :param symbol: The chemical symbol of the Atom. It has to be registered in the ATOMS_DATABASE.
        :type symbol: str

        :param name: The name of the Atom. If this is not provided, the symbol is used as the name as well
        :type nameL str

        :param bonds: List of Atom objects that this Atom is chemically bonded to.
        :type bonds: list

        :param groups: List of groups that this Atom is a part of, e.g. sidechain.
        :type groups: list

        :param ghost:
        :type ghost: bool

        :param kwargs: Any additional parameters that should be set during instantiation.
        """

        super(Atom, self).__init__()

        self._symbol = symbol

        if self._symbol not in ATOMS_DATABASE:
            raise UnknownAtomError("The atom {} is unknown".format(self.symbol))

        self._name = name if name else symbol

        self._bonds = bonds if bonds else []

        self._groups = groups if groups else []

        self._ghost = ghost

        self._index = kwargs.pop("index", None)

        self._parent = None

        self.element = ATOMS_DATABASE.get_atom_property(self._symbol, "element")

        for k, v in kwargs.items():
            try:
                setattr(self, k, v)
            except AttributeError:
                raise AttributeError(
                    f"Could not set attribute {k} to value {v}, probably because this is a protected "
                    f"attribute of this class."
                )

    def __hash__(self) -> int:
        text = self._symbol
        number = self._index
        temp = text + "_" + str(number)
        return temp.__hash__()

    def copy(self) -> "Atom":
        """
        Creates a copy of the current instance of this class.

        :return: A copy of the current instance.
        :rtype: MDANSE.Chemistry.ChemicalEntity.Atom
        """

        a = Atom(symbol=self._symbol)

        for k, v in self.__dict__.items():
            setattr(a, k, v)

        a._bonds = [bat.name for bat in self._bonds]

        return a

    def restore_bonds(self, atom_dict: dict[str, Atom]) -> List[Tuple[int, int]]:
        """After copying, the Atom._bonds is filled with atom NAMES.
        This method uses a dictionary of name: Atom pairs to
        replace the names with Atom instances.
        Trying to copy the Atom instances directly will result
        in infinite recursion.

        Arguments:
            atom_dict -- dictionary of str: atom pairs,
                where the key is the name of the Atom instance
        """
        new_bonds = [atom_dict[atm.name] for atm in self._bonds]
        self._bonds = new_bonds
        return [(self.index, other.index) for other in self._bonds]

    def __eq__(self, other):
        if not isinstance(other, Atom):
            return False
        if not self._index == other._index:
            return False
        if not self._symbol == other._symbol:
            return False
        if not self._name == other._name:
            return False
        return True

    def __getitem__(self, item):
        return getattr(self, item)

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__ = state

    def __str__(self):
        return self.full_name

    def __repr__(self):
        contents = ""
        for key, value in self.__dict__.items():
            key = key[1:] if key[0] == "_" else key
            if key == "bonds":
                bonds = ", ".join(
                    [
                        f'Atom({atom.name if hasattr(atom, "name") else atom})'
                        for atom in self.bonds
                    ]
                )
                contents += f"bonds=[{bonds}]"
            elif isinstance(value, _ChemicalEntity):
                class_name = str(type(value)).replace("<class '", "").replace("'>", "")
                contents += f"{key}={class_name}({value.name})"
            else:
                contents += f"{key}={repr(value)}"
            contents += ", "

        return f"MDANSE.Chemistry.ChemicalEntity.Atom({contents[:-2]})"

    @property
    def atom_list(self) -> list["Atom"]:
        return [self] if not self.ghost else []

    @property
    def total_number_of_atoms(self) -> int:
        return 1

    @property
    def number_of_atoms(self) -> int:
        return int(self.ghost)

    @property
    def bonds(self) -> list["Atom"]:
        """A list of atoms to which this atom is chemically bonded."""
        return self._bonds

    @bonds.setter
    def bonds(self, bonds: list["Atom"]) -> None:
        self._bonds = bonds

    @property
    def ghost(self) -> bool:
        return self._ghost

    @ghost.setter
    def ghost(self, ghost: bool) -> None:
        self._ghost = ghost

    @property
    def groups(self) -> list[str]:
        """A list of groups to which this atom belongs, e.g. sidechain."""
        return self._groups

    @property
    def index(self) -> int:
        """The index of the atom in the trajectory. Once set, it cannot be changed anymore."""
        return self._index

    @index.setter
    def index(self, index: int) -> None:
        if self._index is not None:
            return
        self._index = index

    @property
    def name(self) -> str:
        """The full name of the atom, e.g. Hydrogen."""
        return self._name

    @name.setter
    def name(self, name: str) -> None:
        self._name = name

    @property
    def symbol(self) -> str:
        """The chemical symbol of the atom. The symbol must be registered in ATOMS_DATABASE."""
        return self._symbol

    @symbol.setter
    def symbol(self, symbol: str) -> None:
        if symbol not in ATOMS_DATABASE:
            raise UnknownAtomError("The atom {} is unknown".format(symbol))

        self._symbol = symbol

    @classmethod
    def build(
        cls,
        h5_contents: Union[None, dict[str, list[list[str]]]],
        symbol: str,
        name: str,
        index: str,
        ghost: bool,
    ) -> Atom:
        """
        Creates an instance of the Atom class. This method is meant to be used when loading a trajectory from disk and
        so may be called when the :class: `ChemicalSystem`.load() method is called.

        :param h5_contents: A parameter present solely for compatibility with other classes' build() methods. It is not
            used.
        :type h5_contents: None or dict

        :param symbol: The chemical symbol of the Atom. It has to be registered in the ATOMS_DATABASE.
        :type symbol: str

        :param name: The name of the Atom. If this is not provided, the symbol is used as the name as well
        :type name: str

        :param index: The unique atom index. Since AtomCluster saves it, Atom must save it too.
        :type index: str

        :param ghost:
        :type ghost: bool
        """
        return cls(symbol=symbol, name=name, index=int(index), ghost=ghost)

    def serialize(self, h5_contents: dict[str, list[list[str]]]) -> tuple[str, int]:
        """
        Serializes the Atom object into a string in preparation of the object being stored on disk.

        :param h5_contents: A dictionary that stores all serialized information for the whole ChemicalSystem.
        :type h5_contents: dict

        :return: A tuple containing the string 'atoms' and the index of the serialization data of this Atom in the
            provided dictionary.
        :rtype: tuple
        """
        h5_contents.setdefault("atoms", []).append(
            [repr(self.symbol), repr(self.name), str(self.index), str(self.ghost)]
        )

        return "atoms", len(h5_contents["atoms"]) - 1


class AtomGroup(_ChemicalEntity):
    """
    An arbitrary selection of atoms that belong to the same chemical system. Unlike in Molecule and AtomCluster, the
    atoms in an AtomGroup do not have to be related in any way other than that they have to exist in the same system.
    Further, AtomGroup does not have the serialize() or the build() method defined, and so it will not be saved on disk
    and cannot be loeaded from disk.
    """

    def __init__(self, atoms: list[Atom]):
        """

        :param atoms: The list of atoms that form this AtomGroup
        :type atoms: list
        """
        super(AtomGroup, self).__init__()

        s = set([at.root_chemical_system for at in atoms])
        if len(s) != 1:
            raise ChemicalEntityError("The atoms comes from different chemical systems")

        self._atoms = atoms

        self._chemical_system = list(s)[0]

    def __repr__(self):
        contents = ""
        for key, value in self.__dict__.items():
            key = key[1:] if key[0] == "_" else key
            if isinstance(value, _ChemicalEntity) and not isinstance(value, Atom):
                class_name = str(type(value)).replace("<class '", "").replace("'>", "")
                contents += f"{key}={class_name}({value.name})"
            else:
                contents += f"{key}={repr(value)}"
            contents += ", "

        return f"MDANSE.MolecularDynamics.ChemicalEntity.AtomGroup({contents[:-2]})"

    def __str__(self):
        return f"AtomGroup consisting of {self.total_number_of_atoms} atoms"

    @property
    def atom_list(self) -> list[Atom]:
        """The list of all non-ghost atoms in the AtomGroup."""
        return list([at for at in self._atoms if not at.ghost])

    def copy(self) -> None:
        """The copy method is not defined for the AtomGroup class; instances of it cannot be copied."""
        pass

    @property
    def number_of_atoms(self) -> int:
        """The number of all non-ghost atoms in the AtomGroup."""
        return len([at for at in self._atoms if not at.ghost])

    @property
    def total_number_of_atoms(self) -> int:
        """The total number of atoms in the AtomGroup, including ghosts."""
        return len(self._atoms)

    @property
    def root_chemical_system(self) -> "ChemicalSystem":
        """
        :return: The chemical system whose part all the atoms in this AtomGroup are.
        :rtype: MDANSE.Chemistry.ChemicalEntity.ChemicalSystem
        """
        return self._chemical_system

    def serialize(self, h5_contents: dict) -> None:
        """The serialize method is not defined for the AtomGroup class; it cannot be saved to disk."""
        pass


class BasicCluster:

    def __init__(self, index_list: List[int], **kwargs):
        self._name = kwargs.get("name", "unknown")
        self._atoms = index_list


class ChemicalSystem(_ChemicalEntity):
    """A collection of all chemical compounds in a trajectory."""

    def __init__(self, name: str = ""):
        """

        :param name: The name of the ChemicalSystem
        :type name: str
        """
        super(ChemicalSystem, self).__init__()

        self._chemical_entities = []

        self._configuration = None

        self._number_of_atoms = 0

        self._total_number_of_atoms = 0

        self._name = name

        self._bonds = []

        self._clusters = {}

        self._labels = {}

        self._atoms = None

        self.rdkit_mol = Chem.RWMol()

    def __repr__(self):
        contents = []
        for key, value in self.__dict__.items():
            if key == "rdkit_mol":
                continue
            contents.append(f'{key[1:] if key[0] == "_" else key}={repr(value)}')

        contents = ", ".join(contents)
        return f"MDANSE.MolecularDynamics.ChemicalEntity.ChemicalSystem({contents})"

    def __str__(self):
        return f"ChemicalSystem {self.name} consisting of {len(self._chemical_entities)} chemical entities"

    def initialise_atoms(self, element_list: List[str]):
        self._atom_indices = [self.add_atom(ATOMS_DATABASE.get_atom_property(symbol, "atomic_number")) for symbol in element_list]
        self._atom_types = element_list
        self._total_number_of_atoms = len(self._atom_indices)

    def add_atom(self, atm_num: int) -> int:
        rdkit_atm = Chem.Atom(atm_num)
        rdkit_atm.SetNumExplicitHs(0)
        rdkit_atm.SetNoImplicit(True)
        return self.rdkit_mol.AddAtom(rdkit_atm)        

    def add_bonds(self, pair_list: List[Tuple[int]]):
        self._bonds += pair_list
        for pair in pair_list:
            self.rdkit_mol.AddBond(pair[0], pair[1], Chem.rdchem.BondType.UNSPECIFIED)
    
    def add_clusters(self, group_list: List[List[int]]):
        for group in group_list:
            atom_list = [self._atom_types[index] for index in group]
            unique_atoms, counts = np.unique(atom_list, return_counts=True)
            name = " ".join([str(unique_atoms[n])+str(counts[n]) for n in range(len(counts))])
            if name not in self._clusters:
                self._clusters[name] = []
            self._clusters[name].append(group)

    def has_substructure_match(self, smarts: str) -> bool:
        """Check if there is a substructure match.

        Parameters
        ----------
        smarts : str
            SMARTS string.

        Returns
        -------
        bool
            True if the there is a substructure match.
        """
        return self.rdkit_mol.HasSubstructMatch(Chem.MolFromSmarts(smarts))

    def get_substructure_matches(
        self, smarts: str, maxmatches: int = 1000000
    ) -> set[int]:
        """Get the indexes which match the smarts string. Note that
        the default bond type in MDANSE is
        Chem.rdchem.BondType.UNSPECIFIED.

        Parameters
        ----------
        smarts : str
            SMARTS string.
        maxmatches : int
            Maximum number of matches used in the GetSubstructMatches
            rdkit method.

        Returns
        -------
        set[int]
            An set of matched atom indices.
        """
        substruct_set = set()
        matches = self.rdkit_mol.GetSubstructMatches(
            Chem.MolFromSmarts(smarts), maxMatches=maxmatches
        )
        for match in matches:
            substruct_set.update(match)
        return substruct_set

    @property
    def atom_list(self) -> list[str]:
        """List of all non-ghost atoms in the ChemicalSystem."""
        return self._atom_types

    @property
    def atoms(self) -> list[Atom]:
        """A list of all non-ghost atoms in the ChemicalSystem, sorted by their index."""
        return self._atom_types

    @property
    def chemical_entities(self) -> Dict[str,List[int]]:
        """
        A list of all chemical entities in this ChemicalSystem. Only the entities registered in the ChemicalSystem will
        be returned; their children will not. I.e., if a Molecule object is added to this ChemicalSystem via the
        add_chemical_entity() method, this property will only return the Molecule object, not its consituent Atom
        objects.
        """
        atoms_in_clusters = np.concatenate(self._clusters.values())
        


    @property
    def configuration(self) -> _Configuration:
        """The Configuration that this ChemicalSystem is associated with."""
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: _Configuration):
        if configuration.chemical_system != self:
            raise InconsistentChemicalSystemError("Mismatch between chemical systems")

        self._configuration = configuration

    def copy(self) -> "ChemicalSystem":
        """
        Copies the instance of ChemicalSystem into a new, identical instance.

        :return: Copy of the ChemicalSystem instance
        :rtype: MDANSE.Chemistry.ChemicalEntity.ChemicalSystem
        """
        cs = ChemicalSystem(self._name)

        cs._parent = self._parent

        cs._chemical_entities = [ce.copy() for ce in self._chemical_entities]
        # for ce in self._chemical_entities:
        #     cs.add_chemical_entity(ce.copy())

        new_atoms = {atom.name: atom for atom in cs.atoms}

        for atom in cs.atoms:
            cs._bonds += atom.restore_bonds(new_atoms)

        cs._number_of_atoms = self._number_of_atoms

        cs._total_number_of_atoms = self._total_number_of_atoms

        for ce in cs._chemical_entities:
            ce._parent = cs

        if self._configuration is not None:
            conf = self._configuration.clone(cs)

            cs._configuration = conf

        if self._atoms is not None:
            _ = cs.atoms
        else:
            cs._atoms = None

        return cs

    def rebuild(self, cluster_list: List[Tuple[int]], selection: List[int] = None):
        """
        Copies the instance of ChemicalSystem into a new, identical instance.

        :param cluster_list: list of tuples of atom indices, one per cluster
        :type List[Tuple[int]]: each element is a tuple of atom indices (int)
        """

        atom_names_before = [atom.name for atom in self.atoms]
        clusters = []

        if selection is not None:
            for cluster_number, index_list in enumerate(cluster_list):
                temp = AtomCluster(
                    "cluster_" + str(cluster_number + 1),
                    [
                        self.atom_list[index]
                        for index in index_list
                        if index in selection
                    ],
                )
                if temp.number_of_atoms > 0:
                    clusters.append(temp)
        else:
            for cluster_number, index_list in enumerate(cluster_list):
                temp = AtomCluster(
                    "cluster_" + str(cluster_number + 1),
                    [self.atom_list[index] for index in index_list],
                )
                clusters.append(temp)

        self._chemical_entities = []

        self._number_of_atoms = 0

        self._total_number_of_atoms = 0

        configuration_before = self.configuration
        atom_names_after = [atom.name for atom in self.atoms]

        for cluster in clusters[::-1]:
            self.add_chemical_entity(cluster)

        if atom_names_before == atom_names_after:
            self._configuration = configuration_before
        else:
            LOG.error(f"Atoms before: {atom_names_before}")
            LOG.error(f"Atoms after: {atom_names_after}")
            raise RuntimeError(
                "ChemicalSystem.rebuild() changed the order of atoms. This needs to be handled!"
            )

    def unique_molecules(self) -> List[str]:
        """Returns the list of unique names in the chemical system"""
        result = np.unique(
            [ce.name for ce in self.chemical_entities if ce.number_of_atoms > 1]
        )
        return list(result)

    def number_of_molecules(self, molecule_name: str) -> int:
        """Returns the number of molecules with the given name in the system"""
        result = [1 for ce in self.chemical_entities if ce.name == molecule_name]
        return len(result)

    def from_element_list(self, elements: List[str]):
        for element in elements:
            ce = Atom(element)
            self.add_chemical_entity(ce)

    def load(self, h5_filename: Union[str, h5py.File]) -> None:
        """
        Loads a ChemicalSystem from an HDF5 file. The HDF5 file must be organised in such a way that it can parsed by
        MDANSE.

        :param h5_filename: The HDF5 file that contains the information about the chemical system. This is usually the
            trajectory file.
        :type h5_filename: str or h5py.File

        :return: None
        """
        h5_classes = {
            "atoms": Atom,
            "atom_clusters": AtomCluster,
            "molecules": Molecule,
            "nucleotides": Nucleotide,
            "nucleotide_chain": NucleotideChain,
            "residue": Residue,
            "peptide_chains": PeptideChain,
            "proteins": Protein,
        }

        try:
            h5_file = h5py.File(h5_filename, "r", libver="latest")
            close_file = True
        except TypeError:
            h5_file = h5_filename
            close_file = False

        grp = h5_file["/chemical_system"]
        self._chemical_entities = []

        skeleton = h5_file["/chemical_system/contents"][:]

        try:
            bonds = np.unique(h5_file["/chemical_system/bonds"], axis=0)
        except KeyError:
            bonds = []

        self._name = grp.attrs["name"]

        h5_contents = {}
        for entity_type, v in grp.items():
            if entity_type == "contents" or entity_type == "bonds":
                continue
            h5_contents[entity_type] = v[:]

        for i, (entity_type, entity_index) in enumerate(skeleton):
            entity_index = int(entity_index)
            entity_type = entity_type.decode("utf-8")

            try:
                entity_class = h5_classes[entity_type]
            except KeyError:
                h5_file.close()
                raise CorruptedFileError(
                    f"Could not create a chemical entity of type {entity_type}. The entity listed"
                    f" in the chemical system contents (located at /chemical_system/contents in "
                    f"the HDF5 file) at index {i} is not recognised as valid entity; {entity_type}"
                    f" should be one of: atoms, atom_clusters, molecules, nucleotides, nucleotide"
                    f"_chains, residues, residue_chains, or proteins."
                )
            try:
                arguments = [
                    literal_eval(arg.decode("utf8"))
                    for arg in h5_contents[entity_type][entity_index]
                ]
            # except AttributeError:
            #     print(f"Wrong argument in entity_type: {entity_type} at index {entity_index} ")
            #     print(f"The entry was: {h5_contents[entity_type][entity_index]}")
            except KeyError:
                raise CorruptedFileError(
                    f"Could not find chemical entity {entity_type}, listed in chemical system "
                    f"contents (/chemical_system/contents) at index {i}, in the chemical system "
                    f"itself (/chemical_system). The chemical_system group in the HDF5 file "
                    f"contains only the following datasets: {h5_contents.keys()}."
                )
            except IndexError:
                raise CorruptedFileError(
                    f"The chemical entity {entity_type}, listed in chemical system contents "
                    f"(/chemical_system/contents) at index {i}, could not be found in the "
                    f"{entity_type} dataset (/chemical_system/{entity_type}) because the "
                    f"index registered in contents, {entity_index}, is out of range of the dataset"
                    f", which contains only {len(h5_contents[entity_type])} elements."
                )
            except (ValueError, SyntaxError, RuntimeError) as e:
                raise CorruptedFileError(
                    f"The data used for reconstructing the chemical system could not be parsed "
                    f"from the HDF5 file. The data located at /chemical_system/{entity_type}["
                    f"{entity_index}] is corrupted.\nThe provided data is: "
                    f"{h5_contents[entity_type][entity_index]}\nThe original error is: {e}"
                )
            finally:
                h5_file.close()

            try:
                ce = entity_class.build(h5_contents, *arguments)
            except CorruptedFileError as e:
                raise CorruptedFileError(str(e).replace("INDEX", str(entity_index)))
            except InconsistentAtomNamesError as e:
                raise CorruptedFileError(
                    f"Could not reconstruct {entity_class} from the HDF5 Trajectory because its "
                    f"constituent atoms recorded in the trajectory are different from those "
                    f"expected of this entity with this code ({arguments[1]}). The entity that "
                    f"raised this error is located in the trajectory at /chemical_system/"
                    f"{entity_type} at index {entity_index} while its constituent atoms are at "
                    f"/chemical_system/atoms at indices {arguments[0]}.\nThe full data in the "
                    f"trajectory of the entity that raised the error is "
                    f"{h5_contents[entity_type][entity_index]}\nThe original error is {e}"
                )
            except IndexError as e:
                raise CorruptedFileError(
                    f"Could not reconstruct {entity_class} from the HDF5 Trajectory because one "
                    f"or more of its constituent atoms are missing from the trajectory. This "
                    f"entity is located in the trajectory at /chemical_system/{entity_type} at "
                    f"index {entity_index}. Its constituent atoms are located at /chemical_system/"
                    f'atoms at indices {arguments[0]}, but only {len(h5_contents["atoms"])} atoms '
                    f"are present in the trajectory.\nThe full data in the trajectory of the "
                    f"entity that raised the error is {h5_contents[entity_type][entity_index]}"
                    f"\nThe original error is {e}"
                )
            except KeyError as e:
                raise CorruptedFileError(
                    f"Could not reconstruct {entity_class} from the HDF5 Trajectory because one of"
                    f" its constituent parts could not be found in the trajectory. This entity is "
                    f"located in the trajectory at /chemical_system/{entity_type} at index "
                    f"{entity_index}.\nThe full data in the trajectory of the entity that raised "
                    f"the error is {h5_contents[entity_type][entity_index]}\nThe original error is"
                    f" {e}"
                )
            except TypeError as e:
                raise CorruptedFileError(
                    f"Could not reconstruct {entity_class} from the HDF5 Trajectory because the "
                    f"data associated with it does not match the expected arguments required for"
                    f" reconstructing this entity. This entity is located in the trajectory at "
                    f"/chemical_system/{entity_type} at index {entity_index} and the associated "
                    f"data is {arguments}.\nThe full data in the trajectory of the entity that "
                    f"raised the error is {h5_contents[entity_type][entity_index]}\nThe original "
                    f"error is {e}."
                )
            except ValueError as e:
                raise CorruptedFileError(
                    f"Could not reconstruct {entity_class} from the HDF5 Trajectory because the "
                    f"data associated with it is in an incorrect format. This entity is located "
                    f"in the trajectory at /chemical_system/{entity_type} at index {entity_index} "
                    f"and the associated data is {arguments}.\nThe full data in the trajectory of "
                    f"the entity that raised the error is {h5_contents[entity_type][entity_index]}"
                    f"\nThe original error is {e}."
                )

            self.add_chemical_entity(ce)

        self._bonds = list(bonds)

        if close_file:
            h5_file.close()

        self._h5_file = None

    @property
    def number_of_atoms(self) -> int:
        """The number of non-ghost atoms in the ChemicalSystem."""
        return self._number_of_atoms

    @property
    def total_number_of_atoms(self) -> int:
        """The number of all atoms in the ChemicalSystem, including ghost ones."""
        return self._total_number_of_atoms

    def serialize(self, h5_file: h5py.File) -> None:
        """
        Serializes the contents of the ChemicalSystem object and stores all the data necessary to reconstruct it into
        the provided HDF5 file.

        :param h5_file: The file into which the ChemicalSystem is saved
        :type h5_file: h5py.File

        :return: None
        """
        string_dt = h5py.special_dtype(vlen=str)

        grp = h5_file.create_group("/chemical_system")

        grp.attrs["name"] = self._name

        h5_contents = {}

        contents = []
        for ce in self._chemical_entities:
            entity_type, entity_index = ce.serialize(h5_contents)
            contents.append((entity_type, str(entity_index)))

        for k, v in h5_contents.items():
            grp.create_dataset(k, data=v, dtype=string_dt)
        grp.create_dataset("contents", data=contents, dtype=string_dt)

        h5_bonds = np.array(self._bonds).astype(np.int32)
        grp.create_dataset("bonds", data=h5_bonds, dtype=np.int32)
