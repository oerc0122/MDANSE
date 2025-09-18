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
import random

import numpy as np
from more_itertools import numeric_range

from MDANSE.Framework.Parameters import Float, Integer, Range, Vector
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors


class LinearLatticeQVectors(LatticeQVectors):
    """Generates vectors randomly on a straight line.

    Only vectors commensurate with the reciprocal
    space lattice vectors will be generated.
    |Q| values for which no valid vectors can
    be found are omitted in the output.

    Vectors within one shell are generated within
    a tolerance limit around a central |Q| value.
    Most calculations will produce one data point
    for |Q| by averaging the results over all
    vectors in the group, which is still called
    a shell.
    """

    seed = Integer(
        minimum=0,
        default=0,
    )
    shells = Range[float](
        minimum=0.0,
        default=numeric_range(0.0, 5.0, 0.5),
    )
    n_vectors = Integer(
        minimum=1,
        default=50,
    )
    width = Float(minimum=1e-6, default=1.0)
    axis = Vector(non_zero=True, dtype=int, default=np.array([1, 0, 0], dtype=int))

    def _generate(self):
        if self.seed != 0:
            np.random.seed(self.seed)
            random.seed(self.seed)

        # The Q vector corresponding to the input hkl.
        qVect = self.hkl_to_qvectors(self.axis, self._unit_cell)

        qMax = self.shells[-1] + 0.5 * self.width

        uMax = np.ceil(qMax / np.linalg.norm(qVect) + 1)

        idxs = np.mgrid[-uMax : uMax + 1]

        vects = np.dot(qVect[:, np.newaxis], idxs[np.newaxis, :])

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
