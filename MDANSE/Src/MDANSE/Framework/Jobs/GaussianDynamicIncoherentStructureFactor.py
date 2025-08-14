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

from MDANSE.Framework.AtomGrouping.grouping import (
    add_grouped_totals,
)
from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum
from MDANSE.Mathematics.Signal import get_spectrum
from MDANSE.MolecularDynamics.Analysis import mean_square_displacement


class GaussianDynamicIncoherentStructureFactor(IJob):
    r"""Computes the dynamic incoherent structure factor in the Gaussian approximation.

    Gaussian approximation is exact for a system of free particles and a system of
    particles undergoing brownian motion. The results of this analysis will be close
    to the Dynamic Incoherent Structure Factor analysis in the limits of very
    short :math:`\mathbf{q}` and very long :math:`\mathbf{q}`, and will differ from
    it for intermediate :math:`\mathbf{q}` values.
    """

    label = "Gaussian Dynamic Incoherent Structure Factor"

    category = (
        "Analysis",
        "Scattering",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    settings = collections.OrderedDict()
    settings["trajectory"] = ("HDFTrajectoryConfigurator", {})
    settings["frames"] = (
        "CorrelationFramesConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
    )
    settings["q_shells"] = (
        "RangeConfigurator",
        {"value_type": float, "include_last": True, "mini": 0.0},
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
            "default": "b_incoherent",
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

        self._n_q_shells = self.configuration["q_shells"]["number"]

        self._n_frames = self.configuration["frames"]["n_frames"]

        self._instr_resolution = self.configuration["instrument_resolution"]

        self._n_omegas = self._instr_resolution["n_omegas"]

        self._k_square = self.configuration["q_shells"]["value"] ** 2

        self.add_ideal_results = (
            self.configuration["instrument_resolution"]["kernel"].lower() != "ideal"
        )

        self.labels = [
            (element, (element,)) for element in self.trajectory.get_natoms()
        ]

        self._output_data.add(
            "gdisf/axes/q",
            "LineOutputVariable",
            self.configuration["q_shells"]["value"],
            units="1/nm",
        )

        self._output_data.add(
            "gdisf/axes/q2", "LineOutputVariable", self._k_square, units="1/nm2"
        )

        self._output_data.add(
            "gdisf/axes/time",
            "LineOutputVariable",
            self.configuration["frames"]["duration"],
            units="ps",
        )
        self._output_data.add(
            "gdisf/res/time_window",
            "LineOutputVariable",
            self._instr_resolution["time_window"],
            units="au",
        )

        self._output_data.add(
            "gdisf/axes/omega",
            "LineOutputVariable",
            self.configuration["instrument_resolution"]["omega"],
            units="rad/ps",
        )
        self._output_data.add(
            "gdisf/res/omega_window",
            "LineOutputVariable",
            self._instr_resolution["omega_window"],
            axis="gdisf/axes/omega",
            units="au",
        )

        for element in self.trajectory.unique_names:
            self._output_data.add(
                f"gdisf/f(q,t)/{element}",
                "SurfaceOutputVariable",
                (self._n_q_shells, self._n_frames),
                axis="gdisf/axes/q|gdisf/axes/time",
                units="au",
            )
            self._output_data.add(
                f"gdisf/s(q,f)/{element}",
                "SurfaceOutputVariable",
                (self._n_q_shells, self._n_omegas),
                axis="gdisf/axes/q|gdisf/axes/omega",
                units="au",
                main_result=True,
                partial_result=True,
            )
            self._output_data.add(
                f"msd/{element}",
                "LineOutputVariable",
                (self._n_frames,),
                axis="gdisf/axes/time",
                units="nm2",
            )
            if self.add_ideal_results:
                self._output_data.add(
                    f"gdisf/s(q,f)/ideal/{element}",
                    "SurfaceOutputVariable",
                    (self._n_q_shells, self._n_omegas),
                    axis="gdisf/axes/q|gdisf/axes/omega",
                    units="au",
                )

        self._output_data.add(
            "gdisf/f(q,t)/total",
            "SurfaceOutputVariable",
            (self._n_q_shells, self._n_frames),
            axis="gdisf/axes/q|gdisf/axes/time",
            units="au",
        )
        self._output_data.add(
            "gdisf/s(q,f)/total",
            "SurfaceOutputVariable",
            (self._n_q_shells, self._n_omegas),
            axis="gdisf/axes/q|gdisf/axes/omega",
            units="au",
            main_result=True,
        )
        self._output_data.add(
            "msd/total",
            "LineOutputVariable",
            (self._n_frames,),
            axis="gdisf/axes/time",
            units="nm2",
        )
        if self.add_ideal_results:
            self._output_data.add(
                "gdisf/s(q,f)/ideal/total",
                "SurfaceOutputVariable",
                (self._n_q_shells, self._n_omegas),
                axis="gdisf/axes/q|gdisf/axes/omega",
                units="au",
            )

        self._atoms = self.trajectory.atom_names

    def run_step(self, index: int):
        """Calculates the GDISF and MSD of an atom.

        Parameters
        ----------
        index : int
            The index of the atom that the calculation will be run over.

        Returns
        -------
        tuple[int, tuple[np.ndarray, np.ndarray]]
            A tuple which contains the job index and a tuple of the
            GDISF and MSD of an atom.
        """

        # get atom index
        atom_index = self.trajectory.atom_indices[index]

        series = self.trajectory.read_atomic_trajectory(
            atom_index,
            first=self.configuration["frames"]["first"],
            last=self.configuration["frames"]["last"] + 1,
            step=self.configuration["frames"]["step"],
        )

        series = self.configuration["projection"]["projector"](series)

        atomic_s_f = np.zeros((self._n_q_shells, self._n_frames), dtype=np.float64)

        msd = mean_square_displacement(
            series, self.configuration["frames"]["n_configs"]
        )

        for i, q2 in enumerate(self._k_square):
            gaussian = np.exp(-msd * q2 / 6.0)
            atomic_s_f[i, :] += gaussian

        return index, (atomic_s_f, msd)

    def combine(self, index: int, x: tuple[np.ndarray, np.ndarray]):
        """Add the results to the output files.

        Parameters
        ----------
        index : int
            The atom index that the calculation was run over.
        x : tuple[np.ndarray, np.ndarray]
            A tuple of the GDISF and MSD of an atom.
        """
        element = self._atoms[self.trajectory.atom_indices[index]]
        atomic_s_f, msd = x
        self._output_data[f"gdisf/f(q,t)/{element}"] += atomic_s_f
        self._output_data[f"msd/{element}"] += msd

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...)
        """
        n_atoms_per_element = self.trajectory.get_natoms()
        for element, number in n_atoms_per_element.items():
            self._output_data[f"gdisf/f(q,t)/{element}"][:] /= number
            self._output_data[f"gdisf/s(q,f)/{element}"][:] = get_spectrum(
                self._output_data[f"gdisf/f(q,t)/{element}"],
                self.configuration["instrument_resolution"]["time_window"],
                self.configuration["instrument_resolution"]["time_step"],
                axis=1,
            )
            self._output_data[f"msd/{element}"][:] /= number
            if self.add_ideal_results:
                self._output_data[f"gdisf/s(q,f)/ideal/{element}"][:] = get_spectrum(
                    self._output_data[f"gdisf/f(q,t)/{element}"],
                    None,
                    self.configuration["instrument_resolution"]["time_step"],
                    axis=1,
                )

        selected_weights, all_weights = self.trajectory.get_weights(
            prop=self.configuration["weights"]["property"]
        )
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
        assign_weights(self._output_data, weight_dict, "gdisf/f(q,t)/%s", self.labels)
        assign_weights(self._output_data, weight_dict, "gdisf/s(q,f)/%s", self.labels)
        if self.add_ideal_results:
            assign_weights(
                self._output_data, weight_dict, "gdisf/s(q,f)/ideal/%s", self.labels
            )

        n_selected = sum(n_atoms_per_element.values())
        n_total = sum(self.trajectory.get_all_natoms().values())
        fact = n_selected / n_total

        self._output_data["gdisf/f(q,t)/total"][:] = (
            weighted_sum(self._output_data, "gdisf/f(q,t)/%s", self.labels)
        ) / fact
        self._output_data["gdisf/s(q,f)/total"][:] = (
            weighted_sum(self._output_data, "gdisf/s(q,f)/%s", self.labels)
        ) / fact
        self._output_data["gdisf/f(q,t)/total"].scaling_factor = fact
        self._output_data["gdisf/s(q,f)/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "gdisf/f(q,t)",
            "SurfaceOutputVariable",
            axis="gdisf/axes/q|gdisf/axes/time",
            units="au",
        )
        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "gdisf/s(q,f)",
            "SurfaceOutputVariable",
            axis="gdisf/axes/q|gdisf/axes/omega",
            units="au",
            main_result=True,
            partial_result=True,
        )

        if self.add_ideal_results:
            self._output_data["gdisf/s(q,f)/ideal/total"][:] = (
                weighted_sum(self._output_data, "gdisf/s(q,f)/ideal/%s", self.labels)
                / fact
            )
            self._output_data["gdisf/s(q,f)/ideal/total"].scaling_factor = fact

            add_grouped_totals(
                self.trajectory,
                self._output_data,
                "gdisf/s(q,f)/ideal",
                "SurfaceOutputVariable",
                axis="gdisf/axes/q|gdisf/axes/omega",
                units="au",
            )

        # since GDISF ~ exp(-msd * q2 / 6.0) the MSD isn't weighted in
        # the exp lets save the MSD with equal weights
        selected_weights, all_weights = self.trajectory.get_weights(prop="equal")
        weight_dict = get_weights(
            selected_weights,
            all_weights,
            n_atoms_per_element,
            self.trajectory.get_all_natoms(),
            1,
        )

        assign_weights(self._output_data, weight_dict, "msd/%s", self.labels)
        self._output_data["msd/total"][:] = (
            weighted_sum(self._output_data, "msd/%s", self.labels) / fact
        )

        self._output_data["msd/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "msd",
            "LineOutputVariable",
            axis="gdisf/axes/time",
            units="nm2",
        )

        self._output_data.write(
            self.configuration["output_files"]["root"],
            self.configuration["output_files"]["formats"],
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
