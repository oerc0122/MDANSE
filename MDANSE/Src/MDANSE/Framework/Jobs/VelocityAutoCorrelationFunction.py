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

from scipy.signal import correlate

from MDANSE.Framework.AtomGrouping.grouping import add_grouped_totals
from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.Framework.Parameters import (
    AtomSelection,
    AtomTransmutation,
    CorrelationWindow,
    FrameSelect,
    GroupingLevel,
    InterpOrder,
    MDANSETrajectory,
    OutputFile,
    Projection,
    RunningMode,
    Weights,
)
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum
from MDANSE.Mathematics.Signal import differentiate


class VelocityAutoCorrelationFunction(IJob):
    r"""Calculates the velocity autocorrelation function of the selected atoms.

    The Velocity AutoCorrelation Function (VACF) is a property describing the dynamics
    of a molecular system. It reveals the underlying nature of the forces acting on
    the system. Its Fourier Transform gives the cartesian density of states for a set
    of atoms.
    """

    label = "Velocity AutoCorrelation Function"

    category = (
        "Analysis",
        "Dynamics",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    trajectory = MDANSETrajectory(
        selection="atom_selection",
        transmutation="atom_transmutation",
        grouping="grouping_level",
    )
    frames = FrameSelect(depends={"trajectory": "trajectory"})
    frame_window = CorrelationWindow(depends={"frames": "frames"})
    interpolation_order = InterpOrder(
        label="Velocities", depends={"trajectory": "trajectory", "frames": "frames"}
    )
    projection = Projection(label="Project coordinates")
    grouping_level = GroupingLevel(depends={"trajectory": "trajectory"})
    atom_selection = AtomSelection(depends={"trajectory": "trajectory"})
    atom_transmutation = AtomTransmutation(depends={"trajectory": "trajectory"})
    weights = Weights(
        depends={
            "trajectory": "trajectory",
        }
    )
    output_files = OutputFile()
    running_mode = RunningMode()

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
            "vacf/axes/time",
            "LineOutputVariable",
            self.frame_window.times,
            units="ps",
        )

        for element in self.trajectory.unique_names:
            self._outputData.add(
                f"vacf/{element}",
                "LineOutputVariable",
                (self.frame_window.n_frames,),
                axis="vacf/axes/time",
                units="nm2/ps2",
                main_result=True,
                partial_result=True,
            )

        self._outputData.add(
            "vacf/total",
            "LineOutputVariable",
            (self.frame_window.n_frames,),
            axis="vacf/axes/time",
            units="nm2/ps2",
            main_result=True,
        )

        self._atoms = self.trajectory.atom_names

    def run_step(self, index):
        """
        Runs a single step of the job.\n

        :Parameters:
            #. index (int): The index of the step.
        :Returns:
            #. index (int): The index of the step.
            #. atomicDOS (np.array): The calculated density of state for atom of index=index
            #. atomicVACF (np.array): The calculated velocity auto-correlation function for atom of index=index
        """

        trajectory = self.trajectory

        # get atom index
        atom_index = self.trajectory.atom_indices[index]

        if self.interpolation_order == 0:
            series = trajectory.read_configuration_trajectory(
                atom_index,
                first=self.frames.index_start,
                last=self.frames.index_stop + 1,
                step=self.frames.index_step,
                variable="velocities",
            )
        else:
            series = trajectory.read_atomic_trajectory(
                atom_index,
                first=self.frames.index_start,
                last=self.frames.index_stop + 1,
                step=self.frames.index_step,
            )

            for axis in range(3):
                series[:, axis] = differentiate(
                    series[:, axis],
                    order=self.interpolation_order,
                    dt=self.frames.time_step,
                )

        series = self.projection.projector(series)

        n_configs = self.frame_window.n_configs
        atomicVACF = correlate(series, series[:n_configs], mode="valid") / (
            3 * n_configs
        )
        return index, atomicVACF.T[0]

    def combine(self, index, x):
        """
        Combines returned results of run_step.

        Parameters
        ----------
        index : int
            The index of the step.
        x : Any
            The returned result(s) of run_step.
        """

        # The symbol of the atom.
        element = self._atoms[self.trajectory.atom_indices[index]]

        self._outputData[f"vacf/{element}"] += x

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        nAtomsPerElement = self.trajectory.get_natoms()
        for element, number in nAtomsPerElement.items():
            self._outputData[f"vacf/{element}"] /= number

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
        assign_weights(self._outputData, weight_dict, "vacf/%s", self.labels)

        n_selected = sum(nAtomsPerElement.values())
        n_total = sum(self.trajectory.get_all_natoms().values())
        fact = n_selected / n_total

        vacfTotal = weighted_sum(self._outputData, "vacf/%s", self.labels)
        self._outputData["vacf/total"][:] = vacfTotal / fact
        self._outputData["vacf/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._outputData,
            "vacf",
            "LineOutputVariable",
            axis="vacf/axes/time",
            units="nm2/ps2",
            main_result=True,
            partial_result=True,
        )

        self._outputData.write(
            self.output_files.root,
            self.output_files.out_format,
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
