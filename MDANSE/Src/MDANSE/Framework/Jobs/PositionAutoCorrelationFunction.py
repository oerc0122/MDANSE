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

import numpy as np
from scipy.signal import correlate

from MDANSE.Framework.AtomGrouping.grouping import add_grouped_totals
from MDANSE.Framework.ConfigDescriptors import (
    AtomSelection,
    AtomTransmutation,
    CorrelationWindow,
    FramesConfigDesc,
    GroupingLevel,
    MDANSETrajectoryFile,
    OutputFileConfigDesc,
    PartialCharge,
    ProjectionConfigDesc,
    RunningModeConfigDesc,
    Weights,
)
from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum


class PositionAutoCorrelationFunction(IJob):
    """Calculates the position autocorrelation function.

    Like the velocity autocorrelation function, but using positions instead of
    velocities.
    """

    label = "Position AutoCorrelation Function"

    category = (
        "Analysis",
        "Dynamics",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    trajectory = MDANSETrajectoryFile()
    frames = FramesConfigDesc(depends={"trajectory": "trajectory"})
    frame_window = CorrelationWindow(depends={"frames": "frames"})
    grouping_level = GroupingLevel(depends={"trajectory": "trajectory"})
    atom_selection = AtomSelection(depends={"trajectory": "trajectory"})
    atom_transmutation = AtomTransmutation(depends={"trajectory": "trajectory"})
    projection = ProjectionConfigDesc(label="Project coordinates")
    atom_charges = PartialCharge(
        depends={"trajectory": "trajectory"},
    )
    weights = Weights()
    output_files = OutputFileConfigDesc()
    running_mode = RunningModeConfigDesc()

    def initialize(self):
        """
        Initialize the input parameters and analysis self variables
        """
        super().initialize()

        self.numberOfSteps = len(self.trajectory.atom_indices)

        self.labels = [
            (element, (element,)) for element in self.trajectory.get_natoms()
        ]

        # Will store the time.
        self._outputData.add(
            "pacf/axes/time",
            "LineOutputVariable",
            self.frames.time,
            units="ps",
        )

        # Will store the mean square displacement evolution.
        for element in self.trajectory.unique_names:
            self._outputData.add(
                f"pacf/{element}",
                "LineOutputVariable",
                (len(self.frames),),
                axis="pacf/axes/time",
                units="nm2",
                main_result=True,
                partial_result=True,
            )

        self._atoms = self.trajectory.atom_names

    def run_step(self, index):
        """
        Runs a single step of the job.\n

        :Parameters:
            #. index (int): The index of the step.
        :Returns:
            #. index (int): The index of the step.
            #. atomicPACF (np.array): The calculated position auto-correlation function for atom index
        """

        # get atom index
        atom_index = self.trajectory.atom_indices[index]

        series = self.trajectory.read_atomic_trajectory(
            atom_index,
            first=self.frames.first_index,
            last=self.frames.last_index + 1,
            step=self.frames.step_index,
        )

        series = series - np.average(series, axis=0)
        series = self.projection(series)

        atomicPACF = correlate(series, series[: self.frame_window], mode="valid") / (
            3 * self.frame_window
        )
        return index, atomicPACF.T[0]

    def combine(self, index, x):
        """
        Combines returned results of run_step.\n
        :Parameters:
            #. index (int): The index of the step.\n
            #. x (any): The returned result(s) of run_step
        """

        # The symbol of the atom.
        element = self._atoms[self.trajectory.atom_indices[index]]

        # The MSD for element |symbol| is updated.
        self._outputData[f"pacf/{element}"] += x

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        nAtomsPerElement = self.trajectory.get_natoms()

        for element, number in list(nAtomsPerElement.items()):
            self._outputData[f"pacf/{element}"] /= number

        selected_weights, all_weights = self.trajectory.get_weights(prop=self.weights)
        if self.weights in ("b_coherent", "b_incoherent"):
            for weights in selected_weights, all_weights:
                for key, value in weights.items():
                    weights[key] = abs(value) ** 2

        weight_dict = get_weights(
            selected_weights,
            all_weights,
            nAtomsPerElement,
            self.trajectory.get_all_natoms(),
            1,
        )
        assign_weights(self._outputData, weight_dict, "pacf/%s", self.labels)

        n_selected = sum(nAtomsPerElement.values())
        n_total = sum(self.trajectory.get_all_natoms().values())
        fact = n_selected / n_total

        pacfTotal = weighted_sum(self._outputData, "pacf/%s", self.labels) / fact
        self._outputData.add(
            "pacf/total",
            "LineOutputVariable",
            pacfTotal,
            axis="pacf/axes/time",
            units="nm2",
            main_result=True,
        )
        self._outputData["pacf/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._outputData,
            "pacf",
            "LineOutputVariable",
            axis="pacf/axes/time",
            units="nm2",
            main_result=True,
            partial_result=True,
        )

        self._outputData.write(
            self.output_files.path,
            self.output_files.out_formats,
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
