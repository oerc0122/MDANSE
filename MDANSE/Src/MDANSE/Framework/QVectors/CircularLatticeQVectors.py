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

import random

import numpy as np
from more_itertools import numeric_range

from MDANSE.Framework.Parameters import Float, Integer, Range, Vector
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors


class CircularLatticeQVectors(LatticeQVectors):
    """Generates Q vectors on a plane.

    Vectors are grouped into annuli (called 'shells')
    based on their length.

    Only vectors commensurate with the reciprocal
    space lattice vectors will be generated.
    |Q| values for which no valid vectors can
    be found are omitted in the output.
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
        non_null=True,
        default=np.array([1.0, 0.0, 0.0]),
    )
    axis_2 = Vector(
        non_null=True,
        default=np.array([0.0, 1.0, 0.0]),
    )

    def _generate(self):
        if self.seed != 0:
            np.random.seed(self.seed)
            random.seed(self.seed)

        hkls = np.transpose(
            [
                self.axis_1,
                self.axis_2,
            ]
        )

        qVects = self.hkl_to_qvectors(hkls, self._unit_cell)

        qMax = self.shells[-1] + 0.5 * self.width

        uvMax = np.ceil([qMax / np.linalg.norm(v) for v in qVects.T]) + 1
        # Enforce integers in uvMax
        uvMax = uvMax.astype(np.int64)

        idxs = np.mgrid[-uvMax[0] : uvMax[0] + 1, -uvMax[1] : uvMax[1] + 1]

        idxs = idxs.reshape(2, (2 * uvMax[0] + 1) * (2 * uvMax[1] + 1))

        vects = np.dot(qVects, idxs)

        dists2 = np.sum(vects**2, axis=0)

        halfWidth = self.width / 2

        if self._status is not None:
            self._status.start(len(self.shells))

        self.q_vectors = {}

        for q in self.shells:
            qmin = max(0, q - halfWidth)

            q2low = qmin * qmin
            q2up = (q + halfWidth) * (q + halfWidth)

            hits = np.where((dists2 >= q2low) & (dists2 <= q2up))[0]

            nHits = len(hits)

            if nHits != 0:
                n = min(nHits, self.n_vectors)

                if nHits > self.n_vectors:
                    hits = random.sample(hits, self.n_vectors)

                self.q_vectors[q] = {}
                self.q_vectors[q]["q_vectors"] = vects[:, hits]
                self.q_vectors[q]["n_q_vectors"] = n
                self.q_vectors[q]["q"] = q
                self.q_vectors[q]["hkls"] = self.qvectors_to_hkl(
                    vects[:, hits], self._unit_cell
                )

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
