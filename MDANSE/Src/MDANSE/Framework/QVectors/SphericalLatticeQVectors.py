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

from MDANSE.Framework.Parameters import (
    Float,
    Integer,
    Range,
)
from MDANSE.Framework.QVectors.LatticeQVectors import LatticeQVectors
from MDANSE.Framework.QVectors.SphericalQVectors import spherical_vectors


class SphericalLatticeQVectors(LatticeQVectors):
    """Generates randomly-selected lattice vectors grouped into spheres.

    Only vectors commensurate with the reciprocal space lattice will be generated.
    If more lattice vectors than requested are available in the requested range
    for a shell, a subset of the vectors will be selected with the probability
    given as a Gaussian function with FWHM of width/2 and centred on the |q|
    which is the nominal centre of the shell.

    |Q| values for which no valid vectors can be found are omitted in the output.
    Most calculations will produce one data point per |Q| by averaging the results
    over all vectors in the shell.
    """

    seed = Integer(minimum=0, default=0)
    shells = Range(minimum=0.0, default=(0, 5.0, 0.5), include_last=True)
    n_vectors = Integer(minimum=1, default=50)
    width = Float(minimum=1e-6, default=1.0)

    def _generate(self):
        if self.seed != 0:
            np.random.seed(self.seed)
            random.seed(self.seed)


        if self._status is not None:
            self._status.start(len(self.shells))

        self.q_vectors = {}

        for q in self.shells:
            q_vectors = spherical_vectors(q, self.width, self.n_vectors)
            lattice_hkl_vectors, weights = self.lattice_vectors_with_weights(
                q_vectors,
                self._unit_cell,
            )
            selection = self.vectors_within_limits(
                self.hkl_to_qvectors(lattice_hkl_vectors, self._unit_cell),
                q_min=q - 0.5 * width,
                q_max=q + 0.5 * width,
            )
            weights = weights[selection]
            lattice_hkl_vectors = lattice_hkl_vectors.T[selection].T

            if not len(weights):
                self.q_vectors[q] = None
                continue

            self.q_vectors[q] = {
                "q_vectors": self.hkl_to_qvectors(lattice_hkl_vectors, self._unit_cell),
                "n_q_vectors": np.sum(weights),
                "weights": weights,
                "q": q,
                "hkls": lattice_hkl_vectors,
            }

            if self._status is not None:
                if self._status.is_stopped():
                    return
                self._status.update()
