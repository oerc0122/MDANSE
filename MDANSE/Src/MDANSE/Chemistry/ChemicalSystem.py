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
from typing import TYPE_CHECKING, List, Tuple
import copy

import h5py
import numpy as np
from rdkit import Chem
from MDANSE.MLogging import LOG

if TYPE_CHECKING:
    from MDANSE.MolecularDynamics.Configuration import _Configuration
    from MDANSE.MolecularDynamics.Trajectory import Trajectory


class BasicCluster:

    def __init__(self, index_list: List[int], **kwargs):
        self._name = kwargs.get("name", "unknown")
        self._atoms = index_list


class ChemicalSystem:

    def __init__(self, name: str = "", trajectory=None):
        """

        :param name: The name of the ChemicalSystem
        :type name: str
        """

        self._name = name
        self._trajectory = trajectory

        self._atom_types = []
        self._atom_indices = []
        self._labels = {}  # arbitrary tag attached to atoms (e.g. residue name)

        self._bonds = []

        self._clusters = {}

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
        return f"ChemicalSystem {self._name} consisting of {len(self._atom_types)} atoms in {len(self._clusters)} molecules"

    def initialise_atoms(self, element_list: List[str]):
        self._atom_indices = [
            self.add_atom(self._trajectory.get_atom_property(symbol, "atomic_number"))
            for symbol in element_list
        ]
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
            name = " ".join(
                [str(unique_atoms[n]) + str(counts[n]) for n in range(len(counts))]
            )
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
    def configuration(self) -> _Configuration:
        """The Configuration that this ChemicalSystem is associated with."""
        return self._configuration

    @configuration.setter
    def configuration(self, configuration: _Configuration):
        self._configuration = configuration

    def copy(self) -> "ChemicalSystem":
        """
        Copies the instance of ChemicalSystem into a new, identical instance.

        :return: Copy of the ChemicalSystem instance
        :rtype: MDANSE.Chemistry.ChemicalSystem.ChemicalSystem
        """
        cs = ChemicalSystem(self._name)

        for attribute_name, attribute_value in self.__dict__.items():
            if attribute_name in ["rdkit_mol", "_configuration"]:
                continue
            setattr(cs, attribute_name, copy.deepcopy(attribute_value))

        if self._configuration is not None:
            conf = self._configuration.clone(cs)

            cs._configuration = conf

        return cs

    def unique_molecules(self) -> List[str]:
        """Returns the list of unique names in the chemical system"""
        return list[self._clusters.keys()]

    def number_of_molecules(self, molecule_name: str) -> int:
        """Returns the number of molecules with the given name in the system"""
        return len(self._clusters[molecule_name])

    @property
    def number_of_atoms(self) -> int:
        """The number of non-ghost atoms in the ChemicalSystem."""
        return self._total_number_of_atoms

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

        grp = h5_file.create_group("/composition")
        grp.attrs["name"] = self._name

        grp.create_dataset("atom_types", data=self._atom_types, dtype=string_dt)
        grp.create_dataset("atom_indices", data=self._atom_indices)

        grp.create_dataset("bonds", data=np.array(self._bonds))

        label_group = grp.create_group("labels")
        for key, value in self._labels.items():
            label_group.create_dataset(key, data=value)
        clusters_group = grp.create_group("clusters")
        for key, value in self._clusters.items():
            clusters_group.create_dataset(key, data=value)

    def load(self, trajectory_instance: "Trajectory"):

        source = trajectory_instance.file
        self._trajectory = trajectory_instance
        self.rdkit_mol = Chem.RWMol()

        grp = source["/composition"]
        self._name = grp.attrs["name"]

        atom_types = [binary.decode("utf-8") for binary in grp["atom_types"][:]]
        self.initialise_atoms(atom_types)
        old_indices = [int(tag) for tag in grp["atom_indices"][:]]
        if not np.allclose(old_indices, self._atom_indices):
            LOG.error("Atoms got re-indexed on loading the trajectory")

        self.add_bonds([[int(pair[0]), int(pair[1])] for pair in grp["bonds"][:]])

        self._labels = {}
        for label in grp["labels"].keys():
            self._labels[str(label)] = [int(tag) for tag in grp[f"labels/{str(label)}"]]

        for cluster in grp["clusters"].keys():
            self._clusters[str(cluster)] = [
                [int(x) for x in line] for line in grp[f"clusters/{cluster}"]
            ]
