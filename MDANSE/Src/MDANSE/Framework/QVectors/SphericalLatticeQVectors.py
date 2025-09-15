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

from MDANSE.Framework.Parameters import (
    Float,
    Integer,
    Range,
)
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors


class SphericalLatticeQVectors(LatticeQVectors):
    """Generates vectors randomly on a sphere.

    Only vectors commensurate with the reciprocal
    space lattice vectors will be generated.
    |Q| values for which no valid vectors can
    be found are omitted in the output.

    Vectors within one shell are generated within
    a tolerance limit around a central |Q| value.
    Most calculations will produce one data point
    for |Q| by averaging the results over all
    vectors in the shell.
    """

    seed = Integer(minimum=0, default=0)
    shells = Range[float](
        minimum=0.0,
        default=(0, 5.0, 0.5),
        include_last=True,
    )
    n_vectors = Integer(minimum=1, default=50)
    width = Float(
        minimum=1e-6,
        default=1.0,
    )

    def _generate(self):
        if self.seed != 0:
            np.random.seed(self.seed)
            random.seed(self.seed)
        qMax = self.shells[-1] + 0.5 * self.width

        hklMax = np.ceil(self.qvectors_to_hkl(qMax * np.eye(3), self._unit_cell)) + 1

        hkl_vects = np.mgrid[
            -hklMax[0, 0] : hklMax[0, 0] + 1,
            -hklMax[1, 1] : hklMax[1, 1] + 1,
            -hklMax[2, 2] : hklMax[2, 2] + 1,
        ]

        hkl_vects = hkl_vects.reshape(
            3,
            np.prod(2 * np.diag(hklMax) + 1, dtype=int),
        )

        vects = self.hkl_to_qvectors(hkl_vects, self._unit_cell)

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

            if nHits:
                n = min(nHits, self.n_vectors)

                if nHits > self.n_vectors:
                    hits = random.sample(sorted(hits), self.n_vectors)

                self.q_vectors[q] = {
                    "q_vectors": vects[:, hits],
                    "n_q_vectors": n,
                    "q": q,
                    "hkls": self.qvectors_to_hkl(vects[:, hits], self._unit_cell),
                }

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
