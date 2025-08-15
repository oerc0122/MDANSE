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

from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum
from MDANSE.Mathematics.Signal import correlation, normalisation_factor


class GeneralAutoCorrelationFunction(IJob):
    """Computes the (auto)correlation for any available trajectory variable."""

    enabled = False

    label = "General AutoCorrelation Function"

    category = (
        "Analysis",
        "Dynamics",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    settings = collections.OrderedDict()
    settings["trajectory"] = ("HDFTrajectoryConfigurator", {})
    settings["frames"] = (
        "FramesConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
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
    settings["trajectory_variable"] = (
        "TrajectoryVariableConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
    )
    settings["weights"] = (
        "WeightsConfigurator",
        {
            "dependencies": {
                "trajectory": "trajectory",
                "atom_selection": "atom_selection",
                "atom_transmutation": "atom_transmutation",
            }
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

        self.labels = [
            (element, (element,)) for element in self.trajectory.get_natoms()
        ]

        # Will store the time.
        self._output_data.add(
            "gacf/axes/time",
            "LineOutputVariable",
            self.configuration["frames"]["duration"],
            units="ps",
        )

        # Will store the mean square displacement evolution.
        for element in self.trajectory.unique_names:
            self._output_data.add(
                f"gacf/{element}",
                "LineOutputVariable",
                (self.configuration["frames"]["number"],),
                axis="gacf/axes/time",
                units="au",
                main_result=True,
                partial_result=True,
            )

        self._output_data.add(
            "gacf/total",
            "LineOutputVariable",
            (self.configuration["frames"]["number"],),
            axis="gacf/axes/time",
            units="au",
            main_result=True,
        )

    def run_step(self, index):
        """
        Runs a single step of the job.\n

        :Parameters:
            #. index (int): the index of the step.
        :Returns:
            #. index (int): the index of the step.
            #. atomic_gacf (np.array): the calculated auto-correlation function for the index
        """

        atom_index = self.trajectory.atom_indices[index]

        series = self.configuration["trajectory"][
            "instance"
        ].read_configuration_trajectory(
            atom_index,
            first=self.configuration["frames"]["first"],
            last=self.configuration["frames"]["last"] + 1,
            step=self.configuration["frames"]["step"],
            variable=self.configuration["trajectory_variable"]["value"],
        )

        atomic_gacf = correlation(series, axis=0, average=1)

        return index, atomic_gacf

    def combine(self, index, x):
        """
        Combines returned results of run_step.\n
        :Parameters:
            #. index (int): The index of the step.\n
            #. x (any): The returned result(s) of run_step
        """

        element = self._atoms[self.trajectory.atom_indices[index]]

        self._output_data[f"gacf/{element}"] += x

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        n_atoms_per_element = self.trajectory.get_natoms()

        for element, number in n_atoms_per_element.items():
            self._output_data[f"gacf/{element}"] /= number

        weights = self.trajectory.get_weights(
            prop=self.configuration["weights"]["property"]
        )
        weight_dict = get_weights(weights, n_atoms_per_element, 1)
        assign_weights(self._output_data, weight_dict, "gacf/%s", self.labels)
        gacf_total = weighted_sum(self._output_data, "gacf/%s", self.labels)

        self._output_data["gacf/total"][:] = gacf_total

        self._output_data.write(
            self.configuration["output_files"]["root"],
            self.configuration["output_files"]["formats"],
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
