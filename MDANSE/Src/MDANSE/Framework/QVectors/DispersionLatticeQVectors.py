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

from MDANSE.Framework.Parameters import Integer, Vector
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors


class DispersionLatticeQVectors(LatticeQVectors):
    """Generates Q vectors along a direction."""

    start = Vector(
        dtype=int,
        default=np.zeros(3),
    )
    direction = Vector(non_zero=True, default=np.array([1.0, 0.0, 0.0]))
    n_steps = Integer(
        minimum=1,
        default=10,
        label="Number of steps",
    )

    def _generate(self):
        hkls = np.array(self.start)[:, np.newaxis] + np.outer(
            self.direction, np.arange(0, self.n_steps)
        )

        # The k matrix (3,n_hkls)
        vects = self.hkl_to_qvectors(hkls, self._unit_cell)

        dists = np.linalg.norm(vects, axis=0)

        if self._status is not None:
            self._status.start(len(dists))

        self.q_vectors = {}

        for i, v in enumerate(dists):
            self.q_vectors[v] = {}
            self.q_vectors[v]["q_vectors"] = vects[:, i][:, np.newaxis]
            self.q_vectors[v]["n_q_vectors"] = 1
            self.q_vectors[v]["q"] = v
            self.q_vectors[v]["hkls"] = hkls[:, i][:, np.newaxis]

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
