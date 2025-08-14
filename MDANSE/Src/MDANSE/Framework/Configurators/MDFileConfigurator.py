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

import itertools
import re
from collections.abc import Iterable

import numpy as np

from MDANSE.Core.Error import Error
from MDANSE.Framework.AtomMapping import AtomLabel
from MDANSE.Framework.Units import measure
from MDANSE.MolecularDynamics.UnitCell import UnitCell

from .FileWithAtomDataConfigurator import FileWithAtomDataConfigurator

HBAR = measure(1.05457182e-34, "kg m2 / s").toval("Da nm2 / ps")
HARTREE = measure(27.2113845, "e_v").toval("Da nm2 / ps2")
BOHR = measure(5.29177210903e-11, "m").toval("nm")


class CASTEPError(Error):
    pass


class MDFileConfigurator(FileWithAtomDataConfigurator):
    """Set a CASTEP .md file as input.

    Class representing a .md file format (documentation can be found at
    https://www.tcm.phy.cam.ac.uk/castep/MD/node13.html). It is used to
    determine the structure of the file (eg. the length of each section)
    and to read the information stored in one frame of the trajectory.
    """

    def parse(self):
        self["instance"] = open(self["filename"], "rb")  # Open the provided file.

        # Skip over the header
        while True:
            line = self["instance"].readline().decode("UTF-8")
            if re.search("END", line):
                self["instance"].readline()
                break
            # If a line storing data is read, something is wrong with the header.
            elif re.match(".*<-- h$", line):
                raise CASTEPError(
                    "The provided input file is corrupted. Due to unexpected END header line, the header"
                    "length could not be determined."
                )

        self._header_size = self["instance"].tell()  # Record the length of the header

        # Prepare a variable storing information about a non-specific frame.
        self._frame_info = {"time_step": [0], "cell_data": [], "data": []}

        self["instance"].readline()  # Skip the line storing time information.
        # Save the length of the line storing time information
        self._frame_info["time_step"].append(self["instance"].tell() - self._header_size)

        while True:
            prev_pos = self["instance"].tell()

            line = self["instance"].readline().decode("UTF-8").strip()

            # If the properties of the cell data have not been determined yet and the current line documents cell data
            if not self._frame_info["cell_data"] and re.match(".*<-- h$", line):
                # Save how far (in character number) the cell data is from the start of the frame
                self._frame_info["cell_data"].append(prev_pos - self._header_size)
                # Skip the next two lines since cell data is always three lines long
                self["instance"].readline()
                self["instance"].readline()
                # Save the length of the cell data
                self._frame_info["cell_data"].append(
                    self["instance"].tell()
                    - self._frame_info["cell_data"][0]
                    - self._header_size
                )

            # If the properties of the positional data have not been stored yet and the line stores this data
            elif not self._frame_info["data"] and re.match(".*<-- R$", line):
                # Save how far (in character number) the positional data is from the start of the frame
                self._frame_info["data"].append(prev_pos - self._header_size)

            if not line:
                # Save the length of a frame minus one line of ionic data
                self._frame_info["data"].append(
                    prev_pos - self._frame_info["data"][0] - self._header_size
                )
                break

        # Save the length of the frame, including a blank line
        self._frame_size = self["instance"].tell() - self._header_size

        # Read the whole ionic data block (positions, velocities, and forces) of the first frame
        self["instance"].seek(self._header_size + self._frame_info["data"][0])
        frame = (
            self["instance"]
            .read(self._frame_info["data"][1])
            .decode("UTF-8")
            .splitlines()
        )
        self["n_atoms"] = (
            len(frame) // 3
        )  # Save the number of atoms (length of positional data)

        # Create a list storing the chemical symbol of the element described on each line of positional data
        tmp = [f.split()[0] for f in frame[: self["n_atoms"]]]
        # Save a list of tuples where each tuple consists of the symbol on the amount of those atoms in the simulation
        self["atoms"] = [
            (element, len(list(group))) for element, group in itertools.groupby(tmp)
        ]

        # Move file handle to the end of the file
        self["instance"].seek(0, 2)
        # Save the number of frames
        self["n_frames"] = (
            self["instance"].tell() - self._header_size
        ) // self._frame_size
        self["instance"].seek(0)  # Move file handle to the beginning of the file

    def read_step(self, step):
        """
        Extracts data from one frame of the trajectory

        :param step: The number of the frame to be read.
        :type step: int

        :return: The time of the chosen frame, the cell vectors, and the positions of all atoms in three different units
        :rtype: (float, tuple, np.array)-tuple
        """

        start = (
            self._header_size + step * self._frame_size
        )  # Determine where the step-th frame starts in the file

        # Move file handle to the starts of the line storing the information about time
        self["instance"].seek(start + self._frame_info["time_step"][0])

        # Read the time stored in the line and convert its units
        time_step = float(
            self["instance"].read(self._frame_info["time_step"][1]).decode("UTF-8")
        )
        time_step *= HBAR / HARTREE

        # Read and process the cell data
        self["instance"].seek(
            start + self._frame_info["cell_data"][0]
        )  # Move to the start of cell data
        unit_cell = (
            self["instance"]
            .read(self._frame_info["cell_data"][1])
            .decode("UTF-8")
            .splitlines()
        )  # Read the cell data by line
        # Generate an array of three vectors where each vector is constructed from its components stored in each line
        unit_cell = np.array(
            [[float(bb) for bb in b.strip().split()[:3]] for b in unit_cell]
        )
        unit_cell *= BOHR
        unit_cell = UnitCell(unit_cell)

        self["instance"].seek(
            start + self._frame_info["data"][0]
        )  # Move to the start of positional data
        # Create an array composed of the data points in each line of the positional data
        config = np.array(
            self["instance"].read(self._frame_info["data"][1]).decode("UTF-8").split(),
            dtype=str,
        )
        config = np.reshape(
            config, (3 * self["n_atoms"], 7)
        )  # Reshape the 1D array so that it is organised by lines
        config = config[:, 2:5].astype(np.float64)  # Extract the coordintates only

        # Convert the units of the positions
        config[0 : self["n_atoms"], :] *= BOHR
        config[self["n_atoms"] : 2 * self["n_atoms"], :] *= BOHR * HARTREE / HBAR
        config[2 * self["n_atoms"] : 3 * self["n_atoms"], :] *= HARTREE / BOHR

        return time_step, unit_cell, config

    def close(self):
        """Closes the file."""
        self["instance"].close()

    def atom_labels(self) -> Iterable[AtomLabel]:
        """
        Yields
        ------
        AtomLabel
            An atom label.
        """
        for atm_label, _ in self["atoms"]:
            yield AtomLabel(atm_label)
