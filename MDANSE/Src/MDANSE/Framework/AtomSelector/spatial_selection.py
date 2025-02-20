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

from typing import Union, Dict, Any, Set

import numpy as np
from scipy.spatial import KDTree

from MDANSE.Chemistry.ChemicalSystem import ChemicalSystem
from MDANSE.MolecularDynamics.Trajectory import Trajectory


def select_positions(
    trajectory: Trajectory, **function_parameters: Dict[str, Any]
) -> Set[int]:
    """Selects atoms based on their positions at a specified frame number.
    Lower and upper limits of x, y and z coordinates can be given as input.

    Parameters
    ----------
    trajectory : Trajectory
        A trajectory instance to which the selection is applied
    function_parameters : Dict[str, Any]
        may include
        "frame_number" : int
        "position_minimum" : np.ndarray[float]
        "position_maximum" : np.ndarray[float]

    Returns
    -------
    Set[int]
        Set of all the atom indices
    """
    coordinates = trajectory.coordinates(function_parameters.get("frame_number", 0))
    lower_limits = np.array(function_parameters.get("position_minimum", 3 * [-np.inf]))
    upper_limits = np.array(function_parameters.get("position_maximum", 3 * [np.inf]))
    mask1 = np.all(coordinates > lower_limits.reshape((1, 3)), axis=1)
    mask2 = np.all(coordinates < upper_limits.reshape((1, 3)), axis=1)
    return set(np.where(np.logical_and(mask1, mask2))[0])


def select_sphere(
    trajectory: Trajectory, **function_parameters: Dict[str, Any]
) -> Set[int]:
    """Selects atoms within a distance from a fixed point in space,
    based on coordinates at a specific frame number.

    Parameters
    ----------
    trajectory : Trajectory
        A trajectory instance to which the selection is applied
    function_parameters : Dict[str, Any]
        may include
        "frame_number" : int
        "sphere_centre" : np.ndarray[float]
        "sphere_radius" : float

    Returns
    -------
    Set[int]
        Set of all the atom indices
    """
    coordinates = trajectory.coordinates(function_parameters.get("frame_number", 0))
    sphere_centre = np.array(function_parameters.get("sphere_centre", 3 * [0.0]))
    sphere_radius = function_parameters.get("sphere_radius", 1.0)
    kdtree = KDTree(coordinates)
    indices = kdtree.query_ball_point(sphere_centre, sphere_radius)
    return set(indices)
