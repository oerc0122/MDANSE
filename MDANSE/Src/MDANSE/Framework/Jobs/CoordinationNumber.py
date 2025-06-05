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

import numpy as np


from MDANSE.Framework.Jobs.DistanceHistogram import DistanceHistogram


class CoordinationNumber(DistanceHistogram):
    """
    The Coordination Number is computed from the pair distribution function for a set of atoms.
    It describes the total number of neighbours, as a function of distance, from a central atom, or the centre of a group of atoms.
    """

    label = "Coordination Number"

    enabled = True

    category = (
        "Analysis",
        "Structure",
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
            "valueType": float,
            "includeLast": True,
            "mini": 0.0,
            "dependencies": {"trajectory": "trajectory"},
        },
    )
    settings["grouping_level"] = (
        "GroupingLevelConfigurator",
        {
            "dependencies": {
                "trajectory": "trajectory",
                "atom_selection": "atom_selection",
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
                "atom_selection": "atom_selection",
                "grouping_level": "grouping_level",
            }
        },
    )
    settings["output_files"] = ("OutputFilesConfigurator", {})
    settings["running_mode"] = ("RunningModeConfigurator", {})

    def finalize(self):
        """
        Finalizes the calculations (e.g. averaging the total term, output files creations ...).
        """

        npoints = len(self.configuration["r_values"]["mid_points"])

        self._outputData.add(
            "r",
            "LineOutputVariable",
            self.configuration["r_values"]["mid_points"],
            units="nm",
        )

        self.labels = self.configuration["grouping_level"].pair_labels(all_pairs=True)
        self.labels_intra = self.configuration["grouping_level"].pair_labels(
            intra=True, all_pairs=True
        )

        if self.indices_intra is not None:
            for label, _ in self.labels_intra:
                self._outputData.add(
                    f"cn_intra_{label}",
                    "LineOutputVariable",
                    (npoints,),
                    axis="r",
                    units="au",
                )
            for label, _ in self.labels:
                self._outputData.add(
                    f"cn_inter_{label}",
                    "LineOutputVariable",
                    (npoints,),
                    axis="r",
                    units="au",
                )
        for label, _ in self.labels:
            self._outputData.add(
                f"cn_{label}",
                "LineOutputVariable",
                (npoints,),
                axis="r",
                units="au",
                main_result=True,
            )

        nFrames = self.configuration["frames"]["number"]

        densityFactor = 4.0 * np.pi * self.configuration["r_values"]["mid_points"]

        shellSurfaces = densityFactor * self.configuration["r_values"]["mid_points"]

        shellVolumes = shellSurfaces * self.configuration["r_values"]["step"]

        self.averageDensity *= 4.0 * np.pi / nFrames

        r2 = self.configuration["r_values"]["mid_points"] ** 2
        dr = self.configuration["r_values"]["step"]

        for k in list(self._concentrations.keys()):
            self._concentrations[k] /= nFrames

        nAtomsPerElement = self.configuration["atom_selection"].get_natoms()

        def calc_func(label_i, label_j):
            ni = nAtomsPerElement[label_i]
            nj = nAtomsPerElement[label_j]

            idi = self.selectedElements.index(label_i)
            idj = self.selectedElements.index(label_j)

            if idi == idj:
                nij = ni**2 / 2.0
            else:
                nij = ni * nj
                if self.indices_intra is not None:
                    self.h_intra[idi, idj] += self.h_intra[idj, idi]
                self.h_total[idi, idj] += self.h_total[idj, idi]

            fact = 2 * nij * nFrames * shellVolumes

            rho_i = self.averageDensity * self._concentrations[label_i]

            if self.indices_intra is not None:
                self.h_intra[idi, idj, :] /= fact
                self.h_total[idi, idj, :] /= fact
                cnIntra = np.add.accumulate(self.h_intra[idi, idj, :] * r2) * dr
                cnTotal = np.add.accumulate(self.h_total[idi, idj, :] * r2) * dr
                cnInter = cnTotal - cnIntra
                return rho_i * cnTotal, rho_i * cnInter, rho_i * cnIntra
            else:
                self.h_total[idi, idj, :] /= fact
                cnTotal = np.add.accumulate(self.h_total[idi, idj, :] * r2) * dr
                return rho_i * cnTotal, None, None

        self.configuration["grouping_level"].update_pair_results(
            calc_func, self._outputData, "cn", all_pairs=True
        )

        self._outputData.write(
            self.configuration["output_files"]["root"],
            self.configuration["output_files"]["formats"],
            str(self),
            self,
        )

        self.configuration["trajectory"]["instance"].close()

        super().finalize()
