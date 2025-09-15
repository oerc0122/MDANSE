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

from MDANSE.Framework.Parameters import Float, Integer, Range
from MDANSE.Framework.QVectors.IQVectors import IQVectors
from MDANSE.Mathematics.Geometry import random_points_on_sphere


class SphericalQVectors(IQVectors):
    """Generates vectors randomly on a sphere.

    Vectors within one shell are generated within
    a tolerance limit around a central |Q| value.
    Most calculations will produce one data point
    for |Q| by averaging the results over all
    vectors in the shell.
    """

    seed = Integer(minimum=0, default=0)
    shells = Range(minimum=0.0, include_last=True, default=(0, 5.0, 0.5))
    n_vectors = Integer(minimum=1, default=50)
    width = Float(minimum=1e-9, default=1.0)

    def _generate(self):
        if self.seed != 0:
            np.random.seed(self.seed)

        self.q_vectors = {}

        if self._status is not None:
            self._status.start(len(self.shells))

        for q in self.shells:
            fact = q * np.sign(
                np.random.uniform(-0.5, 0.5, self.n_vectors)
            ) + self.width * np.random.uniform(-0.5, 0.5, self.n_vectors)

            v = random_points_on_sphere(radius=1.0, nPoints=self.n_vectors)

            self.q_vectors[q] = {
                "q_vectors": fact * v,
                "n_q_vectors": self.n_vectors,
                "q": q,
                "hkls": None,
            }

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
