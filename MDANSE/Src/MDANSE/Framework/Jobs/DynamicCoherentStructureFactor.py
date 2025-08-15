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
from math import sqrt

import numpy as np
from scipy.signal import correlate

from MDANSE.Core.Error import Error
from MDANSE.Framework.AtomGrouping.grouping import (
    add_grouped_totals,
    pair_labels,
)
from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.Framework.QVectors.IQVectors import IQVectors
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum
from MDANSE.Mathematics.Signal import get_spectrum
from MDANSE.MLogging import LOG
from MDANSE.MolecularDynamics.UnitCell import UnitCell


class DynamicCoherentStructureFactorError(Error):
    pass


class DynamicCoherentStructureFactor(IJob):
    r"""Computes the dynamic coherent structure factor :math:`S_{\text{coh}}(\mathbf{q}, \omega)` for a set of atoms.

    It can be compared to experimental data e.g. the energy-integrated, static structure
    factor :math:`S_{\text{coh}}(q)` or the dispersion and intensity of phonons.

    The coherent part is derived from correlations between pairs of atoms.
    This analysis requires the :math:`\mathbf{q}`-vectors to be commensurate
    with the reciprocal lattice of the simulation box.
    """

    label = "Dynamic Coherent Structure Factor"

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
    settings["instrument_resolution"] = (
        "InstrumentResolutionConfigurator",
        {"dependencies": {"trajectory": "trajectory", "frames": "frames"}},
    )
    settings["q_vectors"] = (
        "QVectorsConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
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
            "default": "b_coherent",
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
        """Initialize the input parameters and analysis self variables."""
        super().initialize()

        self.n_steps = self.configuration["q_vectors"]["n_shells"]

        n_q_shells = self.configuration["q_vectors"]["n_shells"]
        if not isinstance(
            self.configuration["q_vectors"]["generator"],
            LatticeQVectors,
        ):
            LOG.warning(
                f"This task should be used with a lattice-based Q vector generator. You have picked {self.configuration['q_vectors']['generator'].__class__}. The results are likely to be incorrect."
            )

        self._n_frames = self.configuration["frames"]["n_frames"]

        self._instr_resolution = self.configuration["instrument_resolution"]

        self._n_omegas = self._instr_resolution["n_omegas"]

        self._output_data.add(
            "dcsf/axes/q",
            "LineOutputVariable",
            self.configuration["q_vectors"]["shells"],
            units="1/nm",
        )

        self._output_data.add(
            "dcsf/axes/time",
            "LineOutputVariable",
            self.configuration["frames"]["duration"],
            units="ps",
        )
        self._output_data.add(
            "dcsf/res/time_window",
            "LineOutputVariable",
            self._instr_resolution["time_window"],
            units="au",
        )

        self._output_data.add(
            "dcsf/axes/omega",
            "LineOutputVariable",
            self._instr_resolution["omega"],
            units="rad/ps",
        )
        self._output_data.add(
            "dcsf/res/omega_window",
            "LineOutputVariable",
            self._instr_resolution["omega_window"],
            axis="dcsf/axes/omega",
            units="au",
        )

        self._indices_per_element = self.trajectory.get_indices()
        self.add_ideal_results = (
            self.configuration["instrument_resolution"]["kernel"].lower() != "ideal"
        )

        self.labels = pair_labels(
            self.trajectory,
        )

        for pair_str, _ in self.labels:
            self._output_data.add(
                f"dcsf/f(q,t)/{pair_str}",
                "SurfaceOutputVariable",
                (n_q_shells, self._n_frames),
                axis="dcsf/axes/q|dcsf/axes/time",
                units="au",
            )
            self._output_data.add(
                f"dcsf/s(q,f)/{pair_str}",
                "SurfaceOutputVariable",
                (n_q_shells, self._n_omegas),
                axis="dcsf/axes/q|dcsf/axes/omega",
                units="au",
                main_result=True,
                partial_result=True,
            )
            if self.add_ideal_results:
                self._output_data.add(
                    f"dcsf/s(q,f)/ideal/{pair_str}",
                    "SurfaceOutputVariable",
                    (n_q_shells, self._n_omegas),
                    axis="dcsf/axes/q|dcsf/axes/omega",
                    units="au",
                )

        self._output_data.add(
            "dcsf/f(q,t)/total",
            "SurfaceOutputVariable",
            (n_q_shells, self._n_frames),
            axis="dcsf/axes/q|dcsf/axes/time",
            units="au",
        )
        self._output_data.add(
            "dcsf/s(q,f)/total",
            "SurfaceOutputVariable",
            (n_q_shells, self._n_omegas),
            axis="dcsf/axes/q|dcsf/axes/omega",
            units="au",
            main_result=True,
        )
        if self.add_ideal_results:
            self._output_data.add(
                "dcsf/s(q,f)/ideal/total",
                "SurfaceOutputVariable",
                (n_q_shells, self._n_omegas),
                axis="dcsf/axes/q|dcsf/axes/omega",
                units="au",
            )

        self._cell_std = 0.0
        try:
            all_cells = [
                self.trajectory.unit_cell(frame)._unit_cell
                for frame in self.configuration["frames"]["value"]
            ]
        except TypeError:
            self._average_unit_cell = None
        else:
            self._average_unit_cell = UnitCell(
                np.mean(
                    all_cells,
                    axis=0,
                ),
            )
            self._cell_std = UnitCell(
                np.std(
                    all_cells,
                    axis=0,
                ),
            )

    def run_step(self, index: int) -> tuple[int, dict[str, np.ndarray] | None]:
        """Run the analysis for a single Q shell.

        Parameters
        ----------
        index : int
            index of the Q vector shell

        Returns
        -------
        int, np.ndarray
            shell index, rho density array

        """
        shell = self.configuration["q_vectors"]["shells"][index]

        if shell not in self.configuration["q_vectors"]["value"]:
            return index, None

        traj = self.trajectory

        n_q_vectors = self.configuration["q_vectors"]["value"][shell][
            "q_vectors"
        ].shape[1]

        rho = {}
        for element in self.trajectory.unique_names:
            rho[element] = np.zeros(
                (self.configuration["frames"]["number"], n_q_vectors),
                dtype=np.complex64,
            )

        cell_present = True
        cell_fixed = True
        # loop over the trajectory time steps
        for i, frame in enumerate(self.configuration["frames"]["value"]):
            unit_cell = traj.unit_cell(frame)
            if unit_cell is None:
                cell_present = False
            elif not np.allclose(
                unit_cell.direct,
                self._average_unit_cell.direct,
            ):
                cell_fixed = False
            if not cell_present:
                q_vectors = self.configuration["q_vectors"]["value"][shell]["q_vectors"]
            else:
                try:
                    hkls = self.configuration["q_vectors"]["value"][shell]["hkls"]
                except KeyError:
                    q_vectors = self.configuration["q_vectors"]["value"][shell][
                        "q_vectors"
                    ]
                else:
                    if hkls is None:
                        q_vectors = self.configuration["q_vectors"]["value"][shell][
                            "q_vectors"
                        ]
                    else:
                        q_vectors = IQVectors.hkl_to_qvectors(hkls, unit_cell)

            coords = traj.configuration(frame)["coordinates"]

            for element, idxs in self._indices_per_element.items():
                selected_coordinates = np.take(coords, idxs, axis=0)
                rho[element][i, :] = np.sum(
                    np.exp(1j * np.dot(selected_coordinates, q_vectors)),
                    axis=0,
                )
        if not cell_present:
            LOG.warning(
                "You are running the DCSF calculation on a trajectory without periodic boundary conditions."
            )
        if not cell_fixed:
            LOG.warning(
                f"The unit cell is VARIABLE with the standard deviation of {self._cell_std}. This analysis should not be used with NPT runs! PLEASE CHECK YOUR RESULTS CAREFULLY."
            )

        return index, rho

    def combine(self, index: int, x: np.ndarray):
        """Add partial results to the final array."""
        if x is not None:
            n_configs = self.configuration["frames"]["n_configs"]
            for pair_str, (label_i, label_j) in self.labels:
                # F_ab(Q,t) = F_ba(Q,t) this is valid as long as
                # n_configs is sufficiently large
                corr = correlate(x[label_i], x[label_j][:n_configs], mode="valid").T[
                    0
                ] / (n_configs * x[label_i].shape[1])
                self._output_data[f"dcsf/f(q,t)/{pair_str}"][index, :] += corr.real

    def finalize(self):
        """Apply weights and write out the results."""
        self.configuration["q_vectors"]["generator"].write_vectors_to_file(
            self._output_data,
        )

        n_atoms_per_element = self.trajectory.get_natoms()

        selected_weights, all_weights = self.trajectory.get_weights(
            prop=self.configuration["weights"]["property"]
        )
        weight_dict = get_weights(
            selected_weights,
            all_weights,
            n_atoms_per_element,
            self.trajectory.get_all_natoms(),
            2,
            conc_exp=0.5,
        )
        assign_weights(self._output_data, weight_dict, "dcsf/f(q,t)/%s", self.labels)
        assign_weights(self._output_data, weight_dict, "dcsf/s(q,f)/%s", self.labels)
        if self.add_ideal_results:
            assign_weights(
                self._output_data, weight_dict, "dcsf/s(q,f)/ideal/%s", self.labels
            )
        for pair_str, (label_i, label_j) in self.labels:
            ni = n_atoms_per_element[label_i]
            nj = n_atoms_per_element[label_j]
            self._output_data[f"dcsf/f(q,t)/{pair_str}"] /= sqrt(ni * nj)
            self._output_data[f"dcsf/s(q,f)/{pair_str}"][:] = get_spectrum(
                self._output_data[f"dcsf/f(q,t)/{pair_str}"],
                self.configuration["instrument_resolution"]["time_window"],
                self.configuration["instrument_resolution"]["time_step"],
                axis=1,
            )
            if self.add_ideal_results:
                self._output_data[f"dcsf/s(q,f)/ideal/{pair_str}"][:] = get_spectrum(
                    self._output_data[f"dcsf/f(q,t)/{pair_str}"],
                    None,
                    self.configuration["instrument_resolution"]["time_step"],
                    axis=1,
                )

        n_selected = sum(n_atoms_per_element.values())
        n_total = len(self.trajectory.atom_types)
        fact = n_selected / n_total

        self._output_data["dcsf/f(q,t)/total"][:] = (
            weighted_sum(self._output_data, "dcsf/f(q,t)/%s", self.labels) / fact
        )
        self._output_data["dcsf/f(q,t)/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "dcsf/f(q,t)",
            "SurfaceOutputVariable",
            dim=2,
            conc_exp=0.5,
            axis="dcsf/axes/q|dcsf/axes/time",
            units="au",
        )

        self._output_data["dcsf/s(q,f)/total"][:] = (
            weighted_sum(self._output_data, "dcsf/s(q,f)/%s", self.labels) / fact
        )
        self._output_data["dcsf/s(q,f)/total"].scaling_factor = fact

        add_grouped_totals(
            self.trajectory,
            self._output_data,
            "dcsf/s(q,f)",
            "SurfaceOutputVariable",
            dim=2,
            conc_exp=0.5,
            axis="dcsf/axes/q|dcsf/axes/omega",
            units="au",
            main_result=True,
            partial_result=True,
        )

        if self.add_ideal_results:
            self._output_data["dcsf/s(q,f)/ideal/total"][:] = (
                weighted_sum(self._output_data, "dcsf/s(q,f)/ideal/%s", self.labels)
                / fact
            )
            self._output_data["dcsf/s(q,f)/ideal/total"].scaling_factor = fact
            add_grouped_totals(
                self.trajectory,
                self._output_data,
                "dcsf/s(q,f)/ideal",
                "SurfaceOutputVariable",
                dim=2,
                conc_exp=0.5,
                axis="dcsf/axes/q|dcsf/axes/omega",
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
