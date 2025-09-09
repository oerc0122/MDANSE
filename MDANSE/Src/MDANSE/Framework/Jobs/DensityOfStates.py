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
    InstrumentResolution,
    InterpOrder,
    MDANSETrajectory,
    OutputFile,
    Projection,
    RunningMode,
    Weights,
)
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum
from MDANSE.Mathematics.Signal import differentiate, get_spectrum
from MDANSE.MLogging import LOG


class DensityOfStates(IJob):
    """Calculate the vibrational density of states of the trajectory.

    The Density Of States describes the number of vibrations per unit frequency.
    In MDANSE the DOS calculation returns the Fourier transform (FT) of the weighted
    Velocity AutoCorrelation Function (VACF). With an atomic mass weighting scheme
    the MDANSE DOS result is proportional to the actual vibrational DOS.
    The partial DOS corresponds to selected sets of atoms or molecules.
    """

    label = "Density Of States"

    category = (
        "Analysis",
        "Dynamics",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    trajectory = MDANSETrajectory(
        selection="atom_selection",
        grouping="grouping_level",
        transmutation="atom_transmutation",
    )
    frames = FrameSelect(depends={"trajectory": "trajectory"})
    frame_window = CorrelationWindow(depends={"frames": "frames"})
    grouping_level = GroupingLevel(depends={"trajectory": "trajectory"})
    atom_selection = AtomSelection(depends={"trajectory": "trajectory"})
    atom_transmutation = AtomTransmutation(depends={"trajectory": "trajectory"})
    weights = Weights(
        default="atomic_weight",
        depends={
            "trajectory": "trajectory",
        },
    )
    projection = Projection(
        label="Project coordinates",
    )
    interpolation_order = InterpOrder(
        label="Velocities",
        depends={"trajectory": "trajectory", "frames": "frames"},
    )
    instrument_resolution = InstrumentResolution(
        depends={"window": "frame_window"},
    )
    output_files = OutputFile()
    running_mode = RunningMode()

    def initialize(self):
        """
        Initialize the input parameters and analysis self variables
        """
        super().initialize()

        self.numberOfSteps = self.trajectory.get_total_natoms()

        self.add_ideal_results = self.instrument_resolution.kernel != "ideal"

        self.labels = [
            (element, (element,)) for element in self.trajectory.get_natoms()
        ]

        self._outputData.add(
            "dos/axes/time",
            "LineOutputVariable",
            self.frames.times,
            units="ps",
        )
        self._outputData.add(
            "vacf/axes/time",
            "LineOutputVariable",
            self.frames.times,
            units="ps",
        )

        self._outputData.add(
            "dos/res/time_window",
            "LineOutputVariable",
            self.instrument_resolution.time_window_positive,
            axis="dos/axes/time",
            units="au",
        )

        self._outputData.add(
            "dos/axes/omega",
            "LineOutputVariable",
            self.instrument_resolution.omega,
            units="rad/ps",
        )
        self._outputData.add(
            "dos/axes/romega",
            "LineOutputVariable",
            self.instrument_resolution.romega,
            units="rad/ps",
        )
        self._outputData.add(
            "dos/res/omega_window",
            "LineOutputVariable",
            self.instrument_resolution.omega_window,
            axis="dos/axes/omega",
            units="au",
        )

        for element in self.trajectory.unique_names:
            self._outputData.add(
                f"vacf/{element}",
                "LineOutputVariable",
                (len(self.frames),),
                axis="vacf/axes/time",
                units="nm2/ps2",
            )
            self._outputData.add(
                f"dos/{element}",
                "LineOutputVariable",
                (self.instrument_resolution.n_romegas,),
                axis="dos/axes/romega",
                units="au",
                main_result=True,
                partial_result=True,
            )
            if self.add_ideal_results:
                self._outputData.add(
                    f"dos/ideal/{element}",
                    "LineOutputVariable",
                    (self.instrument_resolution.n_romegas,),
                    axis="dos/axes/romega",
                    units="au",
                )
        self._outputData.add(
            "vacf/total",
            "LineOutputVariable",
            (len(self.frames),),
            axis="dos/axes/time",
            units="nm2/ps2",
        )
        self._outputData.add(
            "dos/total",
            "LineOutputVariable",
            (self.instrument_resolution.n_romegas,),
            axis="dos/axes/romega",
            units="au",
            main_result=True,
        )
        if self.add_ideal_results:
            self._outputData.add(
                "dos/ideal/total",
                "LineOutputVariable",
                (self.instrument_resolution.n_romegas,),
                axis="dos/axes/romega",
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
            #. atomicDOS (np.array): The calculated density of state for atom of index=index
            #. atomicVACF (np.array): The calculated velocity auto-correlation function for atom of index=index
        """
        LOG.debug(f"Running step: {index}")
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

            order = self.interpolation_order
            for axis in range(3):
                series[:, axis] = differentiate(
                    series[:, axis],
                    order=order,
                    dt=self.frames.time_step,
                )

        series = self.projection(series)

        n_configs = len(self.frames)
        atomicVACF = correlate(series, series[:n_configs], mode="valid") / (
            3 * n_configs
        )
        return index, atomicVACF.T[0]

    def combine(self, index, x):
        """
        Combines returned results of run_step.\n
        :Parameters:
            #. index (int): The index of the step.\n
            #. x (any): The returned result(s) of run_step
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
            self._outputData[f"vacf/{element}"][:] /= number
            self._outputData[f"dos/{element}"][:] = get_spectrum(
                self._outputData[f"vacf/{element}"],
                self.instrument_resolution.time_window,
                self.instrument_resolution.time_step,
                fft="rfft",
            )
            if self.add_ideal_results:
                self._outputData[f"dos/ideal/{element}"][:] = get_spectrum(
                    self._outputData[f"vacf/{element}"],
                    None,
                    self.instrument_resolution.time_step,
                    fft="rfft",
                )
        selected_weights, all_weights = self.trajectory.get_weights(prop=self.weights)
        if self.weights in {"b_coherent", "b_incoherent"}:
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
        assign_weights(self._outputData, weight_dict, "dos/%s", self.labels)
        if self.add_ideal_results:
            assign_weights(self._outputData, weight_dict, "dos/ideal/%s", self.labels)

        n_selected = sum(nAtomsPerElement.values())
        n_total = len(self.trajectory.atom_types)
        fact = n_selected / n_total

        self._outputData["vacf/total"][:] = (
            weighted_sum(self._outputData, "vacf/%s", self.labels) / fact
        )
        self._outputData["vacf/total"].scaling_factor = fact
        self._outputData["dos/total"][:] = (
            weighted_sum(self._outputData, "dos/%s", self.labels) / fact
        )
        self._outputData["dos/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._outputData,
            "vacf",
            "LineOutputVariable",
            axis="vacf/axes/time",
            units="nm2/ps2",
        )
        add_grouped_totals(
            self.trajectory,
            self._outputData,
            "dos",
            "LineOutputVariable",
            axis="dos/axes/romega",
            units="au",
            main_result=True,
            partial_result=True,
        )

        if self.add_ideal_results:
            self._outputData["dos/ideal/total"][:] = (
                weighted_sum(self._outputData, "dos/ideal/%s", self.labels) / fact
            )
            self._outputData["dos/ideal/total"].scaling_factor = fact
            add_grouped_totals(
                self.trajectory,
                self._outputData,
                "dos/ideal",
                "LineOutputVariable",
                axis="dos/axes/romega",
                units="au",
            )

        self._outputData.write(
            self.output_files.path,
            self.output_files.out_format,
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
