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
from collections.abc import Iterator

import numpy as np
import numpy.typing as npt

from MDANSE.Framework.AtomGrouping.grouping import (
    add_grouped_totals,
    update_pair_results,
)
from MDANSE.Framework.Jobs.DistanceHistogram import DistanceHistogram
from MDANSE.Mathematics.Arithmetic import assign_weights, get_weights, weighted_sum


def atomic_scattering_factor(element, qvalues, trajectory):
    a = np.empty((4,), dtype=np.float64)
    a[0] = trajectory.get_atom_property(element, "xray_asf_a1")
    a[1] = trajectory.get_atom_property(element, "xray_asf_a2")
    a[2] = trajectory.get_atom_property(element, "xray_asf_a3")
    a[3] = trajectory.get_atom_property(element, "xray_asf_a4")

    b = np.empty((4,), dtype=np.float64)
    b[0] = trajectory.get_atom_property(element, "xray_asf_b1")
    b[1] = trajectory.get_atom_property(element, "xray_asf_b2")
    b[2] = trajectory.get_atom_property(element, "xray_asf_b3")
    b[3] = trajectory.get_atom_property(element, "xray_asf_b4")

    c = trajectory.get_atom_property(element, "xray_asf_c")

    return c + np.sum(
        a[:, np.newaxis]
        * np.exp(-b[:, np.newaxis] * (qvalues[np.newaxis, :] / (4.0 * np.pi)) ** 2),
        axis=0,
    )


class XRayStaticStructureFactor(DistanceHistogram):
    """Computes the X-ray static structure for a set of atoms.

    Computes the X-ray static structure from the pair distribution function for a
    set of atoms, taking into account the atomic form factor for X-rays.
    """

    label = "XRay Static Structure Factor"

    enabled = True

    category = (
        "Analysis",
        "Scattering",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    settings = collections.OrderedDict()
    settings["trajectory"] = ("HDFTrajectoryConfigurator", {})
    settings["frames"] = (
        "FramesConfigurator",
        {"dependencies": {"trajectory": "trajectory"}},
    )
    settings["r_values"] = (
        "DistHistCutoffConfigurator",
        {
            "label": "r values (nm)",
            "value_type": float,
            "include_last": True,
            "mini": 0.0,
            "dependencies": {"trajectory": "trajectory"},
        },
    )
    settings["q_values"] = (
        "RangeConfigurator",
        {"value_type": float, "include_last": True, "mini": 0.0, "default": (0, 500, 1)},
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
    settings["output_files"] = ("OutputFilesConfigurator", {})
    settings["running_mode"] = ("RunningModeConfigurator", {})

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        nq = self.configuration["q_values"]["number"]

        n_frames = self.configuration["frames"]["number"]

        self.average_density /= n_frames

        density_factor = 4.0 * np.pi * self.configuration["r_values"]["mid_points"]

        shell_surfaces = density_factor * self.configuration["r_values"]["mid_points"]

        shell_volumes = shell_surfaces * self.configuration["r_values"]["step"]

        self._output_data.add(
            "xssf/axes/q",
            "LineOutputVariable",
            self.configuration["q_values"]["value"],
            units="1/nm",
        )

        q = self._output_data["xssf/axes/q"]
        r = self.configuration["r_values"]["mid_points"]

        fact1 = 4.0 * np.pi * self.average_density

        sincqr = np.sinc(np.outer(q, r) / np.pi)

        dr = self.configuration["r_values"]["step"]

        for label, _ in self.labels:
            self._output_data.add(
                f"xssf/{label}",
                "LineOutputVariable",
                (nq,),
                axis="xssf/axes/q",
                units="au",
                main_result=True,
                partial_result=True,
            )
        self._output_data.add(
            "xssf/total",
            "LineOutputVariable",
            (nq,),
            axis="xssf/axes/q",
            units="au",
            main_result=True,
        )
        if self.intra:
            for label, _ in self.labels_intra:
                self._output_data.add(
                    f"xssf/intra/{label}",
                    "LineOutputVariable",
                    (nq,),
                    axis="xssf/axes/q",
                    units="au",
                )
            for label, _ in self.labels:
                self._output_data.add(
                    f"xssf/inter/{label}",
                    "LineOutputVariable",
                    (nq,),
                    axis="xssf/axes/q",
                    units="au",
                )
            self._output_data.add(
                "xssf/intra/total",
                "LineOutputVariable",
                (nq,),
                axis="xssf/axes/q",
                units="au",
            )
            self._output_data.add(
                "xssf/inter/total",
                "LineOutputVariable",
                (nq,),
                axis="xssf/axes/q",
                units="au",
            )

        n_atoms_per_element = self.trajectory.get_natoms()

        def calc_func(
            label_i: str, label_j: str
        ) -> Iterator[tuple[str, bool, npt.NDArray]]:
            """Calculates the xray static structure factor for a given
            pair of element labels.

            Parameters
            ----------
            label_i : str
                The element label.
            label_j : str
                The element label.

            Yields
            ------
            name : str
                The results name.
            inter : bool
                Whether results are for intermolecular atom pairs.
            results : npt.NDArray
                The results.
            """
            ni = n_atoms_per_element[label_i]
            nj = n_atoms_per_element[label_j]

            idi = self.selected_elements.index(label_i)
            idj = self.selected_elements.index(label_j)

            if label_i == label_j:
                nij = ni**2 / 2.0
            else:
                nij = ni * nj
                if self.indices_intra is not None:
                    self.h_intra[idi, idj] += self.h_intra[idj, idi]
                self.h_total[idi, idj] += self.h_total[idj, idi]

            fact = 2 * nij * n_frames * shell_volumes

            pdf_total = self.h_total[idi, idj, :] / fact
            yield (
                "xssf",
                False,
                1.0 + fact1 * np.sum((r**2) * (pdf_total - 1.0) * sincqr, axis=1) * dr,
            )

            if self.intra:
                pdf_intra = self.h_intra[idi, idj, :] / fact
                pdf_inter = pdf_total - pdf_intra
                yield (
                    "xssf/inter",
                    False,
                    1.0
                    + fact1 * np.sum((r**2) * (pdf_inter - 1.0) * sincqr, axis=1) * dr,
                )
                yield (
                    "xssf/intra",
                    True,
                    fact1 * np.sum((r**2) * pdf_intra * sincqr, axis=1) * dr,
                )

        update_pair_results(self.trajectory, calc_func, self._output_data)

        asf = {
            name: atomic_scattering_factor(
                ele,
                self._output_data["xssf/axes/q"],
                self.trajectory,
            )
            for name, ele in zip(
                self.trajectory.selection_getter(self.trajectory.atom_names),
                self.trajectory.selection_getter(self.trajectory.atom_types),
            )
        }
        all_asf = {
            name: atomic_scattering_factor(
                ele,
                self._output_data["xssf/axes/q"],
                self.trajectory,
            )
            for name, ele in zip(self.trajectory.atom_names, self.trajectory.atom_types)
        }
        weight_dict = get_weights(
            asf,
            all_asf,
            n_atoms_per_element,
            self.trajectory.get_all_natoms(),
            2,
        )

        n_selected = sum(n_atoms_per_element.values())
        n_total = sum(self.trajectory.get_all_natoms().values())
        fact = (n_selected / n_total) ** 2

        if self.intra:
            assign_weights(
                self._output_data, weight_dict, "xssf/intra/%s", self.labels_intra
            )
            assign_weights(self._output_data, weight_dict, "xssf/inter/%s", self.labels)
            assign_weights(self._output_data, weight_dict, "xssf/%s", self.labels)

            xssf_intra = weighted_sum(
                self._output_data, "xssf/intra/%s", self.labels_intra
            )
            self._output_data["xssf/intra/total"][:] = xssf_intra / fact
            xssf_inter = weighted_sum(self._output_data, "xssf/inter/%s", self.labels)
            self._output_data["xssf/inter/total"][:] = xssf_inter / fact
            self._output_data["xssf/total"][:] = (xssf_intra + xssf_inter) / fact
            self._output_data["xssf/intra/total"].scaling_factor = fact
            self._output_data["xssf/inter/total"].scaling_factor = fact
            self._output_data["xssf/total"].scaling_factor = fact
            for i in ("/intra", "/inter", ""):
                add_grouped_totals(
                    self.trajectory,
                    self._output_data,
                    f"xssf{i}",
                    "LineOutputVariable",
                    dim=2,
                    intra=i == "/intra",
                    axis="xssf/axes/q",
                    units="au",
                    main_result=i == "",
                    partial_result=i == "",
                )

        else:
            assign_weights(self._output_data, weight_dict, "xssf/%s", self.labels)
            self._output_data["xssf/total"][:] = (
                weighted_sum(self._output_data, "xssf/%s", self.labels) / fact
            )
            self._output_data["xssf/total"].scaling_factor = fact

        self._output_data.write(
            self.configuration["output_files"]["root"],
            self.configuration["output_files"]["formats"],
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
