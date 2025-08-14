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
from MDANSE.Mathematics.LinearAlgebra import Vector


class CircularLatticeQVectors(LatticeQVectors):
    """Generates Q vectors on a plane.

    Vectors are grouped into annuli (called 'shells')
    based on their length.

    Only vectors commensurate with the reciprocal
    space lattice vectors will be generated.
    |Q| values for which no valid vectors can
    be found are omitted in the output.
    """

    settings = collections.OrderedDict()
    settings["seed"] = ("IntegerConfigurator", {"mini": 0, "default": 0})
    settings["shells"] = (
        "RangeConfigurator",
        {
            "value_type": float,
            "include_last": True,
            "mini": 0.0,
            "default": (0.0, 5.0, 0.5),
        },
    )
    settings["n_vectors"] = ("IntegerConfigurator", {"mini": 1, "default": 50})
    settings["width"] = ("FloatConfigurator", {"mini": 1.0e-6, "default": 1.0})
    settings["axis_1"] = (
        "VectorConfigurator",
        {"normalize": False, "not_null": True, "value_type": int, "default": [1, 0, 0]},
    )
    settings["axis_2"] = (
        "VectorConfigurator",
        {"normalize": False, "not_null": True, "value_type": int, "default": [0, 1, 0]},
    )

    def _generate(self):
        if self._configuration["seed"]["value"] != 0:
            np.random.seed(self._configuration["seed"]["value"])
            random.seed(self._configuration["seed"]["value"])

        hkls = np.transpose(
            [
                self._configuration["axis_1"]["vector"],
                self._configuration["axis_2"]["vector"],
            ]
        )

        q_vects = self.hkl_to_qvectors(hkls, self._unit_cell)

        q_max = (
            self._configuration["shells"]["last"]
            + 0.5 * self._configuration["width"]["value"]
        )

        uv_max = np.ceil([q_max / Vector(v).length() for v in q_vects.T]) + 1
        # Enforce integers in uv_max
        uv_max = uv_max.astype(np.int64)

        idxs = np.mgrid[-uv_max[0] : uv_max[0] + 1, -uv_max[1] : uv_max[1] + 1]

        idxs = idxs.reshape(2, (2 * uv_max[0] + 1) * (2 * uv_max[1] + 1))

        vects = np.dot(q_vects, idxs)

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
                    hits = random.sample(hits, n_vectors)

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
