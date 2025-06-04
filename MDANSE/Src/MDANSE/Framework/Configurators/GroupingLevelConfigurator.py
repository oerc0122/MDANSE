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
from typing import Callable, Optional
import itertools as it

import numpy as np

from MDANSE.Framework.Configurators.SingleChoiceConfigurator import (
    SingleChoiceConfigurator,
)
from MDANSE.Mathematics.Arithmetic import weighted_sum
from MDANSE.Framework.OutputVariables.IOutputVariable import OutputData


class GroupingLevelConfigurator(SingleChoiceConfigurator):
    """
    This configurator allows to choose the level of granularity in the atom selection.

    When reading the trajectory, the level of granularity will be applied by grouping the atoms of the selection
    to a single dummy-atoms located on the center of gravity of those atoms.

    The level of granularity currently supported are:

    * 'atom': no grouping will be performed
    * 'group': the atoms that belongs to an AtomCluster object will be grouped as a single atom per object while the ones that belongs to a Molecule, NucleotideChain, PeptideChain and Protein object will be grouped according to the chemical group they belong to (e.g. peptide group, methyl group ...)
    * 'residue': the atoms that belongs to anAtomCluster or Molecule object will be grouped as a single atom per object while the ones thta belongs to a NucleotideChain, PeptideChain or Protein object will be grouped according to the residue to which they belong to (e.g. Histidine, Cytosyl ...)
    * 'chain': the atoms that belongs to an AtomCluster or Molecule object will be grouped as a single atom per object while the ones that belongs to a NucleotideChain, PeptideChain or Protein object will be grouped according to the chain they belong to
    * 'molecule': the atoms that belongs to any chemical entity will be grouped as a single atom per object
    """

    _default = "atom"

    def __init__(self, name, choices=None, **kwargs):
        """
        Initializes the configurator.

        :param name: the name of the configurator as it will appear in the configuration
        :type name: str
        :param choices: the level of granularities allowed for the input value. If None all levels are allowed.
        :type choices: one of ['atom','group','residue','chain','molecule'] or None
        """
        usual_choices = ["atom", "molecule"]

        if choices is None:
            choices = usual_choices
        else:
            choices += [x for x in usual_choices if x not in choices]

        SingleChoiceConfigurator.__init__(self, name, choices=choices, **kwargs)

    def configure(self, value: str):
        """
        Parameters
        ----------
        value : str
            The level of granularity at which the atoms should be grouped
        """
        self._original_input = value

        if value is None:
            value = "atom"

        value = str(value)

        SingleChoiceConfigurator.configure(self, value)

        self["level"] = value

        if value == "atom":
            return

        trajConfig = self._configurable[self._dependencies["trajectory"]]
        atomSelectionConfig = self._configurable[self._dependencies["atom_selection"]]
        chemical_system = trajConfig["instance"].chemical_system
        indices = []
        flatten_indices = []
        elements = []
        names = []
        masses = []
        group_names = []
        group_elements = {}
        group_n_atms = {}
        mass_lookup = chemical_system.atom_property("atomic_weight")

        if value == "molecule":
            if len(trajConfig["instance"].chemical_system.unique_molecules()) == 0:
                self.error_status = "The trajectory does not contain molecules."
                return

            for mol_name in chemical_system._clusters.keys():
                n_atms = 0
                group_name_elements = []
                mol_selected = False
                for mol_number, cluster in enumerate(
                    chemical_system._clusters[mol_name]
                ):
                    for x in cluster:
                        if x not in atomSelectionConfig["flatten_indices"]:
                            continue
                        mol_selected = True
                        indices.append([x])
                        flatten_indices.append(x)
                        elements.append([chemical_system.atom_list[x]])
                        group_name_elements.append(chemical_system.atom_list[x])
                        names.append(f"[{mol_name}]_{chemical_system.atom_list[x]}")
                        masses.append([mass_lookup[x]])
                        n_atms += 1
                if mol_selected:
                    group_names.append(mol_name)
                    group_elements[mol_name] = group_name_elements
                    group_n_atms[mol_name] = n_atms

        atomSelectionConfig["indices"] = indices
        atomSelectionConfig["flatten_indices"] = flatten_indices
        atomSelectionConfig["elements"] = elements
        atomSelectionConfig["masses"] = masses
        atomSelectionConfig["names"] = names
        atomSelectionConfig["selection_length"] = len(names)
        atomSelectionConfig["unique_names"] = sorted(set(names))

        self["group_names"] = sorted(set(group_names))
        self["group_elements"] = group_elements
        self["group_n_atms"] = group_n_atms
        self["name_to_element"] = {}
        for name, element in zip(names, elements):
            self["name_to_element"][name] = element[0]
        if atomSelectionConfig["selection_length"] == 0:
            self.error_status = "This option resulted in nothing being selected in the current trajectory"

    def get_element_from_label(self, label: str) -> str:
        """Returns the element for a given label.

        Parameters
        ----------
        label : str
            The label of the element e.g. [H2_O1]_H

        Returns
        -------
        str
            The element of the inputted label.
        """
        if self["level"] == "atom":
            return label
        return self["name_to_element"][label]

    def add_grouped_totals(
        self,
        output_data: dict[str, np.ndarray],
        result_name: str,
        data_type: str,
        dim: int = 1,
        conc_exp: float = 1.0,
        intra: bool = False,
        **kwargs,
    ):
        """Add the grouped totals to the output data.

        Parameters
        ----------
        output_data : dict[str, np.ndarray]
            Dictionary of data arrays containing analysis results.
        result_name : str
            The name of the results.
        data_type : str
            The plotting type of the data.
        dim : int
            Number of repeats of the elements.
        conc_exp : float
            The exponent the at the product of the concentrations are taken
            to (e.g. (c_i * c_j)**0.5 which is used for DCSF jobs).
        intra: bool
            Add total results for intra results.
        """
        tot_n_atms = self._configurable[self._dependencies["atom_selection"]][
            "selection_length"
        ]

        if self["level"] == "atom":
            return

        if dim == 1:
            for grp in self["group_names"]:
                grp_ele = sorted(set(self["group_elements"][grp]))
                conc = self["group_n_atms"][grp] / tot_n_atms
                labels = [((grp, ele), "") for ele in grp_ele]
                results = (
                    weighted_sum(output_data, result_name + "_[%s]_%s", labels) / conc
                )
                output_data.add(
                    f"{result_name}_[{grp}]_total",
                    data_type,
                    results.shape,
                    **kwargs,
                )
                output_data[f"{result_name}_[{grp}]_total"][...] = results
                output_data[f"{result_name}_[{grp}]_total"].scaling_factor = conc
        elif dim == 2:
            if intra:
                for grp in self["group_names"]:
                    eles = sorted(set(self["group_elements"][grp]))
                    conc = (self["group_n_atms"][grp] / tot_n_atms) ** conc_exp
                    labels = [
                        ((grp, *pair), "")
                        for pair in it.combinations_with_replacement(eles, 2)
                    ]

                    results = (
                        weighted_sum(output_data, result_name + "_[%s]_%s%s", labels)
                        / conc
                    )

                    output_data.add(
                        f"{result_name}_[{grp}]_total",
                        data_type,
                        results.shape,
                        **kwargs,
                    )
                    output_data[f"{result_name}_[{grp}]_total"][...] = results
                    output_data[f"{result_name}_[{grp}]_total"].scaling_factor = conc
                return

            for grp_i, grp_j in it.combinations_with_replacement(
                self["group_names"], 2
            ):
                eles_i = sorted(set(self["group_elements"][grp_i]))
                eles_j = sorted(set(self["group_elements"][grp_j]))
                conc_i = self["group_n_atms"][grp_i] / tot_n_atms
                conc_j = self["group_n_atms"][grp_j] / tot_n_atms
                conc = (conc_i * conc_j) ** conc_exp

                if grp_i == grp_j:
                    iterable = it.combinations_with_replacement(eles_i, 2)
                else:
                    iterable = it.product(eles_i, eles_j)
                labels = [((grp_i, grp_j, *pair), "") for pair in iterable]

                results = (
                    weighted_sum(output_data, result_name + "_[%s][%s]_%s%s", labels)
                    / conc
                )

                output_data.add(
                    f"{result_name}_[{grp_i}][{grp_j}]_total",
                    data_type,
                    results.shape,
                    **kwargs,
                )
                output_data[f"{result_name}_[{grp_i}][{grp_j}]_total"][...] = results
                output_data[
                    f"{result_name}_[{grp_i}][{grp_j}]_total"
                ].scaling_factor = conc
        else:
            raise NotImplementedError("Grouped total for dim > 2 not implemented.")

    def pair_labels(self, intra=False) -> list[tuple[str, tuple[str, str]]]:
        """
        Parameters
        ----------
        intra : bool
            Returns the intra label data if true.

        Returns
        -------
        list[tuple[str, tuple[str, str]]]
            The labels of the results and the labels of the individual
            atoms in a tuple.
        """
        labels = []

        if self["level"] == "atom":
            atom_selection = self._configurable[self._dependencies["atom_selection"]]
            selected_elements = atom_selection["unique_names"]
            element_pairs = sorted(
                it.combinations_with_replacement(selected_elements, 2),
            )
            for ele_i, ele_j in element_pairs:
                labels.append((f"{ele_i}{ele_j}", (ele_i, ele_j)))
            return labels

        if intra:
            for grp in self["group_names"]:
                eles = sorted(set(self["group_elements"][grp]))
                for ele_i, ele_j in it.combinations_with_replacement(eles, 2):
                    pair_label = f"[{grp}]_{ele_i}{ele_j}"
                    label_i = f"[{grp}]_{ele_i}"
                    label_j = f"[{grp}]_{ele_j}"
                    labels.append((pair_label, (label_i, label_j)))
            return labels

        for grp_i, grp_j in it.combinations_with_replacement(self["group_names"], 2):
            eles_i = sorted(set(self["group_elements"][grp_i]))
            eles_j = sorted(set(self["group_elements"][grp_j]))
            if grp_i == grp_j:
                iterable = it.combinations_with_replacement(eles_i, 2)
            else:
                iterable = it.product(eles_i, eles_j)
            for ele_i, ele_j in iterable:
                pair_label = f"[{grp_i}][{grp_j}]_{ele_i}{ele_j}"
                label_i = f"[{grp_i}]_{ele_i}"
                label_j = f"[{grp_j}]_{ele_j}"
                labels.append((pair_label, (label_i, label_j)))
        return labels

    def update_pair_results(
        self,
        calc_func: Callable[
            [str, str], tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]
        ],
        output_data: OutputData,
        result_name: str,
    ):
        """Updates the output data with pair results.

        Parameters
        ----------
        calc_func : Callable[[str, str], tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]]
            A function which calculates the total, inter and intra
            molecular results, given two atom labels.
        output_data : OutputData
            The output data object to write the results to.
        result_name : str
            The name of the results.
        """
        if self["level"] == "atom":
            atom_selection = self._configurable[self._dependencies["atom_selection"]]
            selected_elements = atom_selection["unique_names"]
            element_pairs = sorted(
                it.combinations_with_replacement(selected_elements, 2),
            )
            for ele_i, ele_j in element_pairs:
                total, inter, intra = calc_func(ele_i, ele_j)
                output_data[f"{result_name}_{ele_i}{ele_j}"][...] = total
                if intra is not None and inter is not None:
                    output_data[f"{result_name}_inter_{ele_i}{ele_j}"][...] = inter
                    output_data[f"{result_name}_intra_{ele_i}{ele_j}"][...] = intra
            return

        for grp_i, grp_j in it.combinations_with_replacement(self["group_names"], 2):
            eles_i = sorted(set(self["group_elements"][grp_i]))
            eles_j = sorted(set(self["group_elements"][grp_j]))
            if grp_i == grp_j:
                iterable = it.combinations_with_replacement(eles_i, 2)
            else:
                iterable = it.product(eles_i, eles_j)
            for ele_i, ele_j in iterable:
                label_i = f"[{grp_i}]_{ele_i}"
                label_j = f"[{grp_j}]_{ele_j}"
                total, inter, intra = calc_func(label_i, label_j)
                if inter is None or intra is None:
                    raise ValueError(
                        f"Grouping level {self['level']} was used, so atom "
                        f"groups exist but there were no inter and intra "
                        f"results."
                    )
                output_data[f"{result_name}_[{grp_i}][{grp_j}]_{ele_i}{ele_j}"][...] = (
                    total
                )
                output_data[f"{result_name}_inter_[{grp_i}][{grp_j}]_{ele_i}{ele_j}"][
                    ...
                ] = inter
                if grp_i == grp_j:
                    output_data[f"{result_name}_intra_[{grp_i}]_{ele_i}{ele_j}"][
                        ...
                    ] = intra
