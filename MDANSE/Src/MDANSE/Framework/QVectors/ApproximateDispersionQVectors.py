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
import itertools
import operator

import numpy as np

from MDANSE.Framework.Parameters import Float, Vector
from MDANSE.Framework.Parameters.Parameters import ConfigError
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors


class ApproximateDispersionQVectors(LatticeQVectors):
    """Generates Q vectors along a path between two points."""

    q_start = Vector(
        label="Q start (nm^-1)",
        non_null=False,
        default=np.array([0, 0, 0]),
    )
    q_end = Vector(
        label="Q end (nm^-1)",
        non_null=False,
        default=np.array([1, 0, 0]),
    )
    q_step = Float(
        label="Q step (nm^-1)",
        minimum=1e-6,
        default=0.1,
    )

    def validate(self, value, deps):
        if np.isclose(np.linalg.norm(self.q_end - self.q_start), 0.0):
            raise ConfigError("Zero-length vector cannot be used here")

        return value

    def _generate(self):
        d = np.linalg.norm(self.q_end - self.q_start)
        n = (self.q_end - self.q_start).normal()

        nSteps = int(d / self.q_step) + 1

        vects = (
            np.array(self.q_start)[:, np.newaxis]
            + np.outer(n, np.arange(0, nSteps)) * self.q_step
        )

        hkls = self.qvectors_to_hkl(vects, self._unit_cell)

        dists = np.linalg.norm(vects, axis=0)

        dists = sorted(enumerate(dists), key=operator.itemgetter(1))

        qGroups = itertools.groupby(dists, key=operator.itemgetter(1))
        qGroups = {k: [item[0] for item in v] for k, v in qGroups}

        if self._status is not None:
            self._status.start(len(qGroups))

        self.q_vectors = {}

        for k, v in qGroups.items():
            self.q_vectors[k] = {}
            self.q_vectors[k]["q"] = k
            self.q_vectors[k]["q_vectors"] = vects[:, v]
            self.q_vectors[k]["n_q_vectors"] = len(v)
            self.q_vectors[k]["hkls"] = hkls[:, v]

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
