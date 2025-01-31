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
import collections
import itertools
from typing import List

import numpy as np

from MDANSE.Core.Error import Error
from MDANSE.Framework.Jobs.IJob import IJob


class NeutronDynamicTotalStructureFactorError(Error):
    pass


def create_fake_b(
    b_coh: List[float], b_inc: List[float], atom_counts: List[int]
) -> List[float]:
    """Creates random values of b (scattering length) for atoms in the
    simulation. By definition, b_coh is the average of b,
    and b_inc is the standard deviation of b.
    Therefore, a random b value can be assigned to each atom
    which reproduces the b_coh and b_inc of individual atoms,
    at the same time allowing to calculate the average
    b_inc for the entire system

    Parameters
    ----------
    b_coh : List[float]
        list of coherent scattering length values per element
    b_inc : List[float]
        list of incoherent scattering length values per element
    atom_counts : List[int]
        total number of atoms of each element

    Returns
    -------
    List[float]
        a list of random scattering lengths b for each atom
    """
    new_b = []
    for n in range(len(atom_counts)):
        new_b += list(np.random.normal(b_coh[n], b_inc[n], atom_counts[n]))
    return new_b


class NeutronDynamicTotalStructureFactor(IJob):
    """
    Computes the dynamic total structure factor for a set of atoms as the sum of the incoherent and coherent structure factors
    """

    enabled = True

    label = "Neutron Dynamic Total Structure Factor"

    category = (
        "Analysis",
        "Scattering",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    settings = collections.OrderedDict()
    settings["trajectory"] = ("HDFTrajectoryConfigurator", {})
    settings["dcsf_input_file"] = (
        "HDFInputFileConfigurator",
        {"label": "MDANSE Coherent Structure Factor", "default": "dcsf.mda"},
    )
    settings["disf_input_file"] = (
        "HDFInputFileConfigurator",
        {"label": "MDANSE Incoherent Structure Factor", "default": "disf.mda"},
    )
    settings["atom_selection"] = (
        "AtomSelectionConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
    )
    settings["atom_transmutation"] = (
        "AtomTransmutationConfigurator",
        {
            "dependencies": {
                "trajectory": "trajectory",
                "atom_selection": "atom_selection",
            }
        },
    )
    settings["output_files"] = (
        "OutputFilesConfigurator",
        {"formats": ["MDAFormat", "TextFormat"]},
    )
    settings["running_mode"] = ("RunningModeConfigurator", {})

    def initialize(self):
        """
        Initialize the input parameters and analysis self variables
        """
        super().initialize()

        self.numberOfSteps = 1

        # Check time consistency
        if "time" not in self.configuration["dcsf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No time found in dcsf input file"
            )
        if "time" not in self.configuration["disf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No time found in disf input file"
            )

        dcsf_time = self.configuration["dcsf_input_file"]["instance"]["time"][:]
        disf_time = self.configuration["disf_input_file"]["instance"]["time"][:]

        if not np.all(dcsf_time == disf_time):
            raise NeutronDynamicTotalStructureFactorError(
                "Inconsistent times between dcsf and disf input files"
            )

        self._outputData.add("time", "LineOutputVariable", dcsf_time, units="ps")

        # Check time window consistency
        if "time_window" not in self.configuration["dcsf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No time window found in dcsf input file"
            )
        if "time_window" not in self.configuration["disf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No time window found in disf input file"
            )

        dcsf_time_window = self.configuration["dcsf_input_file"]["instance"][
            "time_window"
        ][:]
        disf_time_window = self.configuration["disf_input_file"]["instance"][
            "time_window"
        ][:]

        if not np.all(dcsf_time_window == disf_time_window):
            raise NeutronDynamicTotalStructureFactorError(
                "Inconsistent time windows between dcsf and disf input files"
            )

        self._outputData.add(
            "time_window", "LineOutputVariable", dcsf_time_window, units="au"
        )

        # Check q values consistency
        if "q" not in self.configuration["dcsf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No q values found in dcsf input file"
            )
        if "q" not in self.configuration["disf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No q values found in disf input file"
            )

        dcsf_q = self.configuration["dcsf_input_file"]["instance"]["q"][:]
        disf_q = self.configuration["disf_input_file"]["instance"]["q"][:]

        if not np.all(dcsf_q == disf_q):
            raise NeutronDynamicTotalStructureFactorError(
                "Inconsistent q values between dcsf and disf input files"
            )

        self._outputData.add("q", "LineOutputVariable", dcsf_q, units="1/nm")

        # Check omega consistency
        if "omega" not in self.configuration["dcsf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No omega found in dcsf input file"
            )
        if "omega" not in self.configuration["disf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No omega found in disf input file"
            )

        dcsf_omegas = self.configuration["dcsf_input_file"]["instance"]["omega"][:]
        disf_omegas = self.configuration["disf_input_file"]["instance"]["omega"][:]

        if not np.all(dcsf_omegas == disf_omegas):
            raise NeutronDynamicTotalStructureFactorError(
                "Inconsistent omegas between dcsf and disf input files"
            )

        self._outputData.add("omega", "LineOutputVariable", dcsf_omegas, units="rad/ps")

        # Check omega window consistency
        if "omega_window" not in self.configuration["dcsf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No omega window found in dcsf input file"
            )
        if "omega_window" not in self.configuration["disf_input_file"]["instance"]:
            raise NeutronDynamicTotalStructureFactorError(
                "No omega window found in disf input file"
            )

        dcsf_omega_window = self.configuration["dcsf_input_file"]["instance"][
            "omega_window"
        ][:]
        disf_omega_window = self.configuration["disf_input_file"]["instance"][
            "omega_window"
        ][:]

        if not np.all(dcsf_omega_window == disf_omega_window):
            raise NeutronDynamicTotalStructureFactorError(
                "Inconsistent omega windows between dcsf and disf input files"
            )

        self._outputData.add(
            "omega_window", "LineOutputVariable", dcsf_omegas, units="au"
        )

        # Check f(q,t) and s(q,f) for dcsf
        self._elementsPairs = sorted(
            itertools.combinations_with_replacement(
                self.configuration["atom_selection"]["unique_names"], 2
            )
        )
        for pair in self._elementsPairs:
            if (
                "f(q,t)_{}{}".format(*pair)
                not in self.configuration["dcsf_input_file"]["instance"]
            ):
                raise NeutronDynamicTotalStructureFactorError(
                    "Missing f(q,t) in dcsf input file"
                )
            if (
                "s(q,f)_{}{}".format(*pair)
                not in self.configuration["dcsf_input_file"]["instance"]
            ):
                raise NeutronDynamicTotalStructureFactorError(
                    "Missing s(q,f) in dcsf input file"
                )
            if (
                "scaling_factor"
                not in self.configuration["dcsf_input_file"]["instance"][
                    "s(q,f)_{}{}".format(*pair)
                ].attrs.keys()
            ):
                raise NeutronDynamicTotalStructureFactorError(
                    "This DCSF file was created before the new scaling scheme. Please calculate it again."
                )

        for element in self.configuration["atom_selection"]["unique_names"]:
            if (
                "f(q,t)_{}".format(element)
                not in self.configuration["disf_input_file"]["instance"]
            ):
                raise NeutronDynamicTotalStructureFactorError(
                    "Missing f(q,t) in disf input file"
                )
            if (
                "s(q,f)_{}".format(element)
                not in self.configuration["disf_input_file"]["instance"]
            ):
                raise NeutronDynamicTotalStructureFactorError(
                    "Missing s(q,f) in disf input file"
                )
            if (
                "scaling_factor"
                not in self.configuration["disf_input_file"]["instance"][
                    "s(q,f)_{}".format(element)
                ].attrs.keys()
            ):
                raise NeutronDynamicTotalStructureFactorError(
                    "This DISF file was created before the new scaling scheme. Please calculate it again."
                )

        for element in self.configuration["atom_selection"]["unique_names"]:
            fqt = self.configuration["disf_input_file"]["instance"][
                "f(q,t)_{}".format(element)
            ]
            sqf = self.configuration["disf_input_file"]["instance"][
                "s(q,f)_{}".format(element)
            ]
            self._outputData.add(
                "f(q,t)_inc_%s" % element,
                "SurfaceOutputVariable",
                fqt,
                axis="q|time",
                units="au",
            )
            self._outputData.add(
                "s(q,f)_inc_%s" % element,
                "SurfaceOutputVariable",
                sqf,
                axis="q|omega",
                units="nm2/ps",
            )
            self._outputData.add(
                "f(q,t)_inc_weighted_%s" % element,
                "SurfaceOutputVariable",
                fqt.shape,
                axis="q|time",
                units="au",
            )
            self._outputData.add(
                "s(q,f)_inc_weighted_%s" % element,
                "SurfaceOutputVariable",
                sqf.shape,
                axis="q|omega",
                units="nm2/ps",
            )

        for pair in self._elementsPairs:
            fqt = self.configuration["dcsf_input_file"]["instance"][
                "f(q,t)_{}{}".format(*pair)
            ]
            sqf = self.configuration["dcsf_input_file"]["instance"][
                "s(q,f)_{}{}".format(*pair)
            ]
            self._outputData.add(
                "f(q,t)_coh_%s%s" % pair,
                "SurfaceOutputVariable",
                fqt,
                axis="q|time",
                units="au",
            )
            self._outputData.add(
                "s(q,f)_coh_%s%s" % pair,
                "SurfaceOutputVariable",
                sqf,
                axis="q|omega",
                units="nm2/ps",
            )
            self._outputData.add(
                "f(q,t)_coh_weighted_%s%s" % pair,
                "SurfaceOutputVariable",
                fqt.shape,
                axis="q|time",
                units="au",
            )
            self._outputData.add(
                "s(q,f)_coh_weighted_%s%s" % pair,
                "SurfaceOutputVariable",
                sqf.shape,
                axis="q|omega",
                units="nm2/ps",
            )

        nQValues = len(dcsf_q)
        nTimes = len(dcsf_time)
        nOmegas = len(dcsf_omegas)

        self._outputData.add(
            "f(q,t)_coh_total",
            "SurfaceOutputVariable",
            (nQValues, nTimes),
            axis="q|time",
            units="au",
        )
        self._outputData.add(
            "f(q,t)_inc_total",
            "SurfaceOutputVariable",
            (nQValues, nTimes),
            axis="q|time",
            units="au",
        )
        self._outputData.add(
            "f(q,t)_total",
            "SurfaceOutputVariable",
            (nQValues, nTimes),
            axis="q|time",
            units="au",
        )

        self._outputData.add(
            "s(q,f)_coh_total",
            "SurfaceOutputVariable",
            (nQValues, nOmegas),
            axis="q|omega",
            units="nm2/ps",
            main_result=True,
            partial_result=True,
        )
        self._outputData.add(
            "s(q,f)_inc_total",
            "SurfaceOutputVariable",
            (nQValues, nOmegas),
            axis="q|omega",
            units="nm2/ps",
            main_result=True,
            partial_result=True,
        )
        self._outputData.add(
            "s(q,f)_total",
            "SurfaceOutputVariable",
            (nQValues, nOmegas),
            axis="q|omega",
            units="nm2/ps",
            main_result=True,
        )
        self._input_disf_weight = (
            self.configuration["disf_input_file"]["instance"][
                "metadata/inputs/weights"
            ][0]
            .decode()
            .strip('"')
        )
        self._input_dcsf_weight = (
            self.configuration["dcsf_input_file"]["instance"][
                "metadata/inputs/weights"
            ][0]
            .decode()
            .strip('"')
        )

    def run_step(self, index):
        """
        Runs a single step of the job.\n

        :Parameters:
            #. index (int): The index of the step.
        :Returns:
            #. index (int): The index of the step.
            #. rho (np.array): The exponential part of I(k,t)
        """

        return index, None

    def combine(self, index, x):
        """
        Combines returned results of run_step.\n
        :Parameters:
            #. index (int): The index of the step.\n
            #. x (any): The returned result(s) of run_step
        """

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """

        nAtomsPerElement = self.configuration["atom_selection"].get_natoms()

        # Compute concentrations
        nTotalAtoms = 0
        for val in list(nAtomsPerElement.values()):
            nTotalAtoms += val

        b_coh_list, b_inc_list, atom_count_list = [], [], []
        for element, ni in nAtomsPerElement.items():
            bcoh = self.configuration["trajectory"]["instance"].get_atom_property(
                element, "b_coherent"
            )
            binc = self.configuration["trajectory"]["instance"].get_atom_property(
                element, "b_incoherent"
            )
            b_coh_list.append(bcoh)
            b_inc_list.append(binc)
            atom_count_list.append(ni)

        fake_b_values = create_fake_b(b_coh_list, b_inc_list, atom_count_list)
        effective_b_coh = np.mean(fake_b_values)
        effective_b_inc = np.std(fake_b_values)

        norm_natoms = 1.0 / nTotalAtoms
        # Compute coherent functions and structure factor
        for pair in self._elementsPairs:
            bi = self.configuration["trajectory"]["instance"].get_atom_property(
                pair[0], "b_coherent"
            )
            bj = self.configuration["trajectory"]["instance"].get_atom_property(
                pair[1], "b_coherent"
            )

            self._outputData["f(q,t)_coh_weighted_%s%s" % pair][:] = (
                self._outputData["f(q,t)_coh_%s%s" % pair][:] * bi * bj
            )
            self._outputData["s(q,f)_coh_weighted_%s%s" % pair][:] = (
                self._outputData["s(q,f)_coh_%s%s" % pair][:] * bi * bj
            )
            if pair[0] == pair[1]:  # Add a factor 2 if the two elements are different
                self._outputData["f(q,t)_coh_total"][:] += self._outputData[
                    "f(q,t)_coh_weighted_%s%s" % pair
                ][:]
                self._outputData["s(q,f)_coh_total"][:] += self._outputData[
                    "s(q,f)_coh_weighted_%s%s" % pair
                ][:]
            else:
                self._outputData["f(q,t)_coh_total"][:] += (
                    2 * self._outputData["f(q,t)_coh_weighted_%s%s" % pair][:]
                )
                self._outputData["s(q,f)_coh_total"][:] += (
                    2 * self._outputData["s(q,f)_coh_weighted_%s%s" % pair][:]
                )

        # Compute incoherent functions and structure factor
        for element, ni in nAtomsPerElement.items():

            self._outputData["f(q,t)_inc_weighted_%s" % element][:] = (
                self._outputData["f(q,t)_inc_%s" % element][:] * effective_b_inc**2
            )
            self._outputData["s(q,f)_inc_weighted_%s" % element][:] = (
                self._outputData["s(q,f)_inc_%s" % element][:] * effective_b_inc**2
            )

            self._outputData["f(q,t)_inc_total"][:] += self._outputData[
                "f(q,t)_inc_weighted_%s" % element
            ][:]
            self._outputData["s(q,f)_inc_total"][:] += self._outputData[
                "s(q,f)_inc_weighted_%s" % element
            ][:]

        sigma_coh = 4 * np.pi * effective_b_coh**2
        sigma_inc = 4 * np.pi * effective_b_inc**2

        # Compute total F(Q,t) = inc + coh
        self._outputData["f(q,t)_total"][:] = (
            self._outputData["f(q,t)_coh_total"][:] * norm_natoms
            + self._outputData["f(q,t)_inc_total"][:] * norm_natoms
        )
        self._outputData["s(q,f)_coh_total"][:] *= (
            norm_natoms * sigma_coh / (sigma_coh + sigma_inc)
        )
        self._outputData["s(q,f)_inc_total"][:] *= (
            norm_natoms * sigma_inc / (sigma_coh + sigma_inc)
        )

        self._outputData["s(q,f)_total"][:] = (
            self._outputData["s(q,f)_coh_total"][:]
            + self._outputData["s(q,f)_inc_total"][:]
        )

        self._outputData.write(
            self.configuration["output_files"]["root"],
            self.configuration["output_files"]["formats"],
            self._info,
            self,
        )
        self.configuration["trajectory"]["instance"].close()
        self.configuration["disf_input_file"]["instance"].close()
        self.configuration["dcsf_input_file"]["instance"].close()
        super().finalize()
