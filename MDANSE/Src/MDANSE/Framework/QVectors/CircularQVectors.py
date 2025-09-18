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

import numpy as np
from more_itertools import numeric_range

from MDANSE.Framework.Parameters import Float, Integer, Range, Vector
from MDANSE.Framework.QVectors.IQVectors import IQVectors
from MDANSE.Mathematics.Geometry import random_points_on_circle


class CircularQVectors(IQVectors):
    """Generates Q vectors on a plane.

    Vectors are grouped into annuli (called 'shells')
    based on their length.
    """

    seed = Integer(
        minimum=0,
        default=0,
    )
    shells = Range[float](minimum=0.0, default=numeric_range(0.0, 5.0, 0.5))
    n_vectors = Integer(
        minimum=1,
        default=50,
    )
    width = Float(
        minimum=1e-6,
        default=1,
    )
    axis_1 = Vector(
        non_zero=True,
        default=np.array([1.0, 0.0, 0.0]),
    )
    axis_2 = Vector(
        non_zero=True,
        default=np.array([0.0, 1.0, 0.0]),
    )

    def _generate(self):
        if self.seed != 0:
            np.random.seed(self.seed)

        axis = np.cross(self.axis_1, self.axis_2)
        axis /= np.linalg.norm(axis)

        if self._status is not None:
            self._status.start(len(self.shells))

        self.q_vectors = {}

        for q in self.shells:
            fact = q * np.sign(
                np.random.uniform(-0.5, 0.5, self.n_vectors)
            ) + self.width * np.random.uniform(-0.5, 0.5, self.n_vectors)
            v = random_points_on_circle(axis, radius=1.0, nPoints=self.n_vectors)

            self.q_vectors[q] = {}
            self.q_vectors[q]["q_vectors"] = fact * v
            self.q_vectors[q]["n_q_vectors"] = self.n_vectors
            self.q_vectors[q]["q"] = q
            self.q_vectors[q]["hkls"] = None

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
