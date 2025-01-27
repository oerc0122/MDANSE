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

from typing import List

import numpy as np

from MDANSE.MolecularDynamics.Configuration import (
    contiguous_coordinates_box,
    remove_jumps,
)


def centre_of_mass(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    scaled_positions = coords * masses.reshape((len(coords), 1))
    return scaled_positions.mean(axis=0) / masses.mean()


def com_single_frame(
    coords: np.ndarray,
    cell: np.ndarray,
    masses: List[float],
    clusters: List[List[int]],
    selection: List[int] = None,
):
    """Calculates the centre of mass of a single trajectory frame.
    Assumes that the input coordinates are FRACTIONAL.

    Parameters
    ----------
    coords : np.ndarray
        Atom coordinates
    cell : np.ndarray
        3x3 unit cell array
    rcell : np.ndarray
        3x3 reciprocal cell array
    masses : List[float]
        array of atom masses
    clusters : List[List[int]]
        groups of atoms to be made contiguous
    selection : List[int]
        list of atoms to include in the centre of mass
    box_coordinates : bool, optional
        use fractional input coordinates, by default False
    """
    new_coordinates = contiguous_coordinates_box(coords, cell, clusters)
    if selection is not None:
        return centre_of_mass(new_coordinates[selection], masses[selection])
    centre_of_mass(new_coordinates, masses)


def com_trajectory(
    coords_vs_time: np.ndarray,
    cells: np.ndarray,
    rcells: np.ndarray,
    masses: List[float],
    clusters,
    selection=None,
    box_coordinates=False,
):

    trajectory = np.empty((coords_vs_time.shape[0], 3), dtype=np.float64)

    for timestep in range(len(trajectory)):
        frac_coordinates = np.matmul(coords_vs_time[timestep], rcells[timestep])
        trajectory[timestep] = com_single_frame(
            frac_coordinates, cells[timestep], masses, clusters, selection
        )

    trajectory = remove_jumps(trajectory)

    if not box_coordinates:
        for timestep in range(len(trajectory)):
            trajectory[timestep] = np.matmul(
                trajectory[timestep : timestep + 1], cells[timestep]
            )

    return trajectory
