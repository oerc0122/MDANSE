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

    settings = collections.OrderedDict()
    settings["seed"] = ("IntegerConfigurator", {"mini": 0, "default": 0})
    settings["shells"] = (
        "RangeConfigurator",
        {
            "value_type": float,
            "include_last": True,
            "mini": 0.0,
            "default": (0, 5.0, 0.5),
        },
    )
    settings["n_vectors"] = ("IntegerConfigurator", {"mini": 1, "default": 50})
    settings["width"] = ("FloatConfigurator", {"mini": 1.0e-6, "default": 1.0})

    def _generate(self):
        if self._configuration["seed"]["value"] != 0:
            np.random.seed(self._configuration["seed"]["value"])
            random.seed(self._configuration["seed"]["value"])
        q_max = (
            self._configuration["shells"]["last"]
            + 0.5 * self._configuration["width"]["value"]
        )

        hkl_max = np.ceil(self.qvectors_to_hkl(q_max * np.eye(3), self._unit_cell)) + 1

        hkl_vects = np.mgrid[
            -hkl_max[0, 0] : hkl_max[0, 0] + 1,
            -hkl_max[1, 1] : hkl_max[1, 1] + 1,
            -hkl_max[2, 2] : hkl_max[2, 2] + 1,
        ]

        hkl_vects = hkl_vects.reshape(
            3,
            np.prod(2 * np.diag(hkl_max) + 1, dtype=int),
        )

        vects = self.hkl_to_qvectors(hkl_vects, self._unit_cell)

        dists2 = np.sum(vects**2, axis=0)

        half_width = self._configuration["width"]["value"] / 2

        n_vectors = self._configuration["n_vectors"]["value"]

        if self._status is not None:
            self._status.start(self._configuration["shells"]["number"])

        self._configuration["q_vectors"] = collections.OrderedDict()

        for q in self._configuration["shells"]["value"]:
            qmin = max(0, q - half_width)

            q2low = qmin * qmin
            q2up = (q + half_width) * (q + half_width)

            hits = np.where((dists2 >= q2low) & (dists2 <= q2up))[0]

            n_hits = len(hits)

            if n_hits != 0:
                n = min(n_hits, n_vectors)

                if n_hits > n_vectors:
                    hits = random.sample(sorted(hits), n_vectors)

                self._configuration["q_vectors"][q] = {}
                self._configuration["q_vectors"][q]["q_vectors"] = vects[:, hits]
                self._configuration["q_vectors"][q]["n_q_vectors"] = n
                self._configuration["q_vectors"][q]["q"] = q
                self._configuration["q_vectors"][q]["hkls"] = self.qvectors_to_hkl(
                    vects[:, hits], self._unit_cell
                )

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
