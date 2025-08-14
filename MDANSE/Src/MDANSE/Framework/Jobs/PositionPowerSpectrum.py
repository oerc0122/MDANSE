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

import collections

import numpy as np
from scipy.signal import correlate

from MDANSE.Framework.AtomGrouping.grouping import (
    add_grouped_totals,
)
from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum
from MDANSE.Mathematics.Signal import get_spectrum
from MDANSE.MLogging import LOG


class PositionPowerSpectrum(IJob):
    """Calculates the position power spectrum.

    Power spectrum (using Fast Fourier Transform) of atomic trajectories calculated
    from the Positional Autocorrelation Function (PACF). The calculation result
    is similar to the density of states, which is calculated from the velocity
    autocorrelation function. In fact, PPS and DOS may show the same features
    in systems where position and velocity are strongly correlated (i.e. solids).

    This calculation is used by TrajectoryFilter to create a preview of the spectrum,
    so that a desired range of atomic vibrational modes can be isolated.
    """

    label = "PositionPowerSpectrum"

    category = (
        "Analysis",
        "Dynamics",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    settings = collections.OrderedDict()
    settings["trajectory"] = ("HDFTrajectoryConfigurator", {})
    settings["frames"] = (
        "CorrelationFramesConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
    )
    settings["instrument_resolution"] = (
        "InstrumentResolutionConfigurator",
        {"dependencies": {"trajectory": "trajectory", "frames": "frames"}},
    )
    settings["projection"] = (
        "ProjectionConfigurator",
        {"label": "project coordinates"},
    )
    settings["grouping_level"] = (
        "GroupingLevelConfigurator",
        {
            "dependencies": {
                "trajectory": "trajectory",
            }
        },
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
            }
        },
    )
    settings["weights"] = (
        "WeightsConfigurator",
        {
            "default": "atomic_weight",
            "dependencies": {
                "trajectory": "trajectory",
                "atom_selection": "atom_selection",
                "atom_transmutation": "atom_transmutation",
            },
        },
    )
    settings["output_files"] = ("OutputFilesConfigurator", {})
    settings["running_mode"] = ("RunningModeConfigurator", {})

    def initialize(self):
        """
        Initialize the input parameters and analysis self variables
        """
        super().initialize()

        self.n_steps = len(self.trajectory.atom_indices)

        instr_resolution = self.configuration["instrument_resolution"]

        self.add_ideal_results = (
            self.configuration["instrument_resolution"]["kernel"].lower() != "ideal"
        )

        self.labels = [
            (element, (element,)) for element in self.trajectory.get_natoms()
        ]

        self._output_data.add(
            "pps/axes/time",
            "LineOutputVariable",
            self.configuration["frames"]["duration"],
            units="ps",
        )
        self._output_data.add(
            "pacf/axes/time",
            "LineOutputVariable",
            self.configuration["frames"]["duration"],
            units="ps",
        )

        self._output_data.add(
            "pps/res/time_window",
            "LineOutputVariable",
            instr_resolution["time_window_positive"],
            axis="pps/axes/time",
            units="au",
        )

        self._output_data.add(
            "pps/axes/omega",
            "LineOutputVariable",
            instr_resolution["omega"],
            units="rad/ps",
        )
        self._output_data.add(
            "pps/axes/romega",
            "LineOutputVariable",
            instr_resolution["romega"],
            units="rad/ps",
        )
        self._output_data.add(
            "pps/res/omega_window",
            "LineOutputVariable",
            instr_resolution["omega_window"],
            axis="pps/axes/omega",
            units="au",
        )

        for element in self.trajectory.unique_names:
            self._output_data.add(
                f"pacf/{element}",
                "LineOutputVariable",
                (self.configuration["frames"]["n_frames"],),
                axis="pacf/axes/time",
                units="nm2",
            )
            self._output_data.add(
                f"pps/{element}",
                "LineOutputVariable",
                (instr_resolution["n_romegas"],),
                axis="pps/axes/romega",
                units="au",
                main_result=True,
                partial_result=True,
            )
            if self.add_ideal_results:
                self._output_data.add(
                    f"pps/ideal/{element}",
                    "LineOutputVariable",
                    (instr_resolution["n_romegas"],),
                    axis="pps/axes/romega",
                    units="au",
                )

        self._output_data.add(
            "pacf/total",
            "LineOutputVariable",
            (self.configuration["frames"]["n_frames"],),
            axis="pacf/axes/time",
            units="nm2",
        )
        self._output_data.add(
            "pps/total",
            "LineOutputVariable",
            (instr_resolution["n_romegas"],),
            axis="pps/axes/romega",
            units="au",
            main_result=True,
        )
        if self.add_ideal_results:
            self._output_data.add(
                "pps/ideal/total",
                "LineOutputVariable",
                (instr_resolution["n_romegas"],),
                axis="pps/axes/romega",
                units="au",
            )

        self._atoms = self.trajectory.atom_names

    def run_step(self, index):
        """
        Runs a single step of the job.\n

        :Parameters:
            #. index (int): The index of the step.
        :Returns:
            #. index (int): The index of the step.
            #. atomic_e_s (np.array): The calculated energy spectrum for atom of index=index
            #. atomic_p_a_c_f (np.array): The calculated position auto-correlation function for atom of index=index
        """
        LOG.debug(f"Running step: {index}")
        trajectory = self.trajectory

        # get atom index
        atom_index = self.trajectory.atom_indices[index]

        series = trajectory.read_atomic_trajectory(
            atom_index,
            first=self.configuration["frames"]["first"],
            last=self.configuration["frames"]["last"] + 1,
            step=self.configuration["frames"]["step"],
        )

        series = series - np.average(series, axis=0)

        series = self.configuration["projection"]["projector"](series)

        n_configs = self.configuration["frames"]["n_configs"]
        atomic_p_a_c_f = correlate(series, series[:n_configs], mode="valid") / (
            3 * n_configs
        )
        return index, atomic_p_a_c_f.T[0]

    def combine(self, index, x):
        """
        Combines returned results of run_step.\n
        :Parameters:
            #. index (int): The index of the step.\n
            #. x (any): The returned result(s) of run_step
        """

        # The symbol of the atom.
        element = self._atoms[self.trajectory.atom_indices[index]]

        self._output_data[f"pacf/{element}"] += x

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        n_atoms_per_element = self.trajectory.get_natoms()
        for element, number in n_atoms_per_element.items():
            self._output_data[f"pacf/{element}"][:] /= number
            self._output_data[f"pps/{element}"][:] = get_spectrum(
                self._output_data[f"pacf/{element}"],
                self.configuration["instrument_resolution"]["time_window"],
                self.configuration["instrument_resolution"]["time_step"],
                fft="rfft",
            )
            if self.add_ideal_results:
                self._output_data[f"pps/ideal/{element}"][:] = get_spectrum(
                    self._output_data[f"pacf/{element}"],
                    None,
                    self.configuration["instrument_resolution"]["time_step"],
                    fft="rfft",
                )

        selected_weights, all_weights = self.trajectory.get_weights(
            prop=self.configuration["weights"]["property"]
        )
        if self.configuration["weights"]["property"] in ("b_coherent", "b_incoherent"):
            for weights in selected_weights, all_weights:
                for key, value in weights.items():
                    weights[key] = value**2
        weight_dict = get_weights(
            selected_weights,
            all_weights,
            n_atoms_per_element,
            self.trajectory.get_all_natoms(),
            1,
        )
        assign_weights(self._output_data, weight_dict, "pacf/%s", self.labels)
        assign_weights(self._output_data, weight_dict, "pps/%s", self.labels)
        if self.add_ideal_results:
            assign_weights(self._output_data, weight_dict, "pps/ideal/%s", self.labels)

        n_selected = sum(n_atoms_per_element.values())
        n_total = sum(self.trajectory.get_all_natoms().values())
        fact = n_selected / n_total

        self._output_data["pacf/total"][:] = (
            weighted_sum(
                self._output_data,
                "pacf/%s",
                self.labels,
            )
            / fact
        )
        self._output_data["pacf/total"].scaling_factor = fact
        self._output_data["pps/total"][:] = (
            weighted_sum(
                self._output_data,
                "pps/%s",
                self.labels,
            )
            / fact
        )
        self._output_data["pps/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "pacf",
            "LineOutputVariable",
            axis="pacf/axes/time",
            units="nm2",
        )
        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "pps",
            "LineOutputVariable",
            axis="pps/axes/romega",
            units="au",
            main_result=True,
            partial_result=True,
        )
        if self.add_ideal_results:
            self._output_data["pps/ideal/total"][:] = (
                weighted_sum(
                    self._output_data,
                    "pps/ideal/%s",
                    self.labels,
                )
                / fact
            )
            self._output_data["pps/ideal/total"].scaling_factor = fact
            add_grouped_totals(
                self.trajectory,
                self._output_data,
                "pps/ideal",
                "LineOutputVariable",
                axis="pps/axes/romega",
                units="au",
            )

        self._output_data.write(
            self.configuration["output_files"]["root"],
            self.configuration["output_files"]["formats"],
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
