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

from MDANSE.Framework.ConfigDescriptors import (
    AtomSelection,
    FramesConfigDesc,
    GroupingLevel,
    MDANSETrajectoryFile,
    OutputFileConfigDesc,
    RunningModeConfigDesc,
)
from MDANSE.Framework.Jobs.IJob import IJob
from MDANSE.MolecularDynamics.Analysis import mean_square_fluctuation


class RootMeanSquareFluctuation(IJob):
    """Calculates the Root Mean Square Fluctuation of atom positions.

    The root mean square fluctuation (RMSF) for a set of atoms is similar to the
    square root of the mean square displacement (MSD), except that it is spatially
    resolved rather than time resolved. It reveals the dynamical heterogeneity
    of the molecule over the course of an MD simulation.

    As opposed to most analysis types, the result is a single number per atom index.
    """

    label = "Root Mean Square Fluctuation"

    category = (
        "Analysis",
        "Dynamics",
    )

    ancestor = ["hdf_trajectory", "molecular_viewer"]

    trajectory = MDANSETrajectoryFile()
    frames = FramesConfigDesc(depends={"trajectory": "trajectory"})
    grouping_level = GroupingLevel(depends={"trajectory": "trajectory"}, default="atom")
    atom_selection = AtomSelection(depends={"trajectory": "trajectory"})
    output_files = OutputFileConfigDesc()
    running_mode = RunningModeConfigDesc()

    def initialize(self):
        """Initialize the input parameters and analysis self variables"""
        super().initialize()
        self.numberOfSteps = len(self.trajectory.atom_indices)

        self.group_molecules = self.grouping_level != "atom"
        self.ele_idxs = {}

        self._names = self.trajectory.atom_names

        self._outputData.add(
            "rmsf/axes/indices/all",
            "LineOutputVariable",
            self.trajectory.atom_indices,
        )
        self._outputData.add(
            "rmsf/all",
            "LineOutputVariable",
            (self.numberOfSteps,),
            axis="rmsf/axes/indices/all",
            units="nm",
            main_result=True,
        )

        for names in self.trajectory.unique_names:
            idxs = [i for i in self.trajectory.atom_indices if names == self._names[i]]
            self.ele_idxs[names] = idxs
            self._outputData.add(
                f"rmsf/axes/indices/{names}",
                "LineOutputVariable",
                idxs,
            )
            self._outputData.add(
                f"rmsf/{names}",
                "LineOutputVariable",
                (len(idxs),),
                axis=f"rmsf/axes/indices/{names}",
                units="nm",
            )

    def run_step(self, index: int) -> tuple[int, float]:
        """Runs a single step of the job.

        Parameters
        ----------
        index : int
            The index of the step.

        Returns
        -------
        intex : int
            The index of the step.
        rmsf : float
            The calculated root mean square fluctuation for atom index.

        """
        # read the particle trajectory

        if not self.group_molecules:
            struct_index = self.trajectory.atom_indices[index]
        else:
            struct_index = self.cluster_lookup[index]

        series = self.trajectory.read_atomic_trajectory(
            struct_index,
            first=self.frames.first_index,
            last=self.frames.last_index + 1,
            step=self.frames.step_index,
        )

        rmsf = mean_square_fluctuation(series, root=True)
        return index, rmsf

    def combine(self, index: int, rmsf: float) -> None:
        """Combines returned results of run_step.

        Parameters
        ----------
        index : int
            The index of the step.
        rmsf : float
            The returned result(s) of run_step

        """
        self._outputData["rmsf/rmsf"][index] = rmsf

    def finalize(self):
        """Finalizes the calculations (e.g. averaging the total term, output
        files creations ...).
        """
        if self.group_molecules:
            for grp in self.trajectory.group_lookup:
                eles = self.trajectory.group_elements(grp)
                idxs = []
                for ele in eles:
                    idxs += self.ele_idxs[f"<{grp}>/{ele}"]
                idxs = sorted(idxs)
                self._outputData.add(
                    f"rmsf/axes/indices/<{grp}>/all",
                    "LineOutputVariable",
                    idxs,
                )
                self._outputData.add(
                    f"rmsf/<{grp}>/all",
                    "LineOutputVariable",
                    (len(idxs),),
                    axis=f"rmsf/axes/indices/<{grp}>/all",
                    units="nm",
                )
                for i, idx in enumerate(idxs):
                    self._outputData[f"rmsf/<{grp}>/all"][i] = self._outputData[
                        "rmsf/all"
                    ][idx]

        # Write the output variables.
        self._outputData.write(
            self.output_files.path,
            self.output_files.out_formats,
            str(self),
            self,
        )

        self.trajectory.close()
        super().finalize()
