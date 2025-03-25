import math
from contextlib import nullcontext as success

import numpy as np
import pytest
from MDANSE.Framework.NewQVectors import (CircularQVectors,
                                          LatticeSphericalQVectors,
                                          QVectorData, SphericalQVectors)
from MDANSE.MolecularDynamics.UnitCell import UnitCell
from test_helpers.common_fixtures import cell, cubic, hexagonal, orthorhombic


@pytest.mark.parametrize("cell, params, expected_hkl", [
    ("cubic", (1.,),
     [[-4, -2, -2], [-4, -2, -1], [-4, -2, 0]],
     ),

    ("orthorhombic", (1.,),
     [[-3, -3, -1], [-3, -3, 0], [-3, -3, 1]]),

    ("hexagonal", (1.,),
                 [[-3, -3, -3], [-3, -3, -2], [-3, -3, -1]],
     ),

    pytest.param("cubic", (1., [1, 0, 0]),
                 [[-4, -2, -2], [-4, -2, -1], [-4, -2, 0]],
                 marks=pytest.mark.xfail(reason="Centre not implemented")),

    pytest.param("orthorhombic", (2., [1, 0, 0]),
                 [[0, 0, 0], [2, 2, 2], [4, 4, 4]],
                 marks=pytest.mark.xfail(reason="Centre not implemented")),

    pytest.param("hexagonal", (1., [1, 0, 0]),
                 [[-3, -3, -3], [-3, -3, -2], [-3, -3, -1]],
                 marks=pytest.mark.xfail(reason="Centre not implemented")),
], indirect=["cell"])
def test_lattice_spherical(cell, params, expected_hkl):
    qvec = LatticeSphericalQVectors(cell, *params)
    results = list(qvec[3])
    assert results == [QVectorData.from_hkl(hkl, cell) for hkl in expected_hkl]


@pytest.mark.parametrize("params", [
    ({"radius": 2.}),
    ({"radius": 1.}),
])
def test_spherical(params):
    qvec = SphericalQVectors(**params)

    for q in qvec[10]:
        print(q)
        assert q.mod_q == pytest.approx(params["radius"])

@pytest.mark.parametrize("params", [
    ({"radius": 2., "axis": [1, 0, 0]}),
    ({"radius": 1., "axis": [1, 1, 0]}),
    ({"radius": 1., "axis": [-1, 1, 1.3]}),
])
def test_circular(params):
    qvec = CircularQVectors(**params)

    for q in qvec[10]:
        assert np.dot(q.q, params["axis"]) == pytest.approx(0.)
        assert q.mod_q == pytest.approx(params["radius"])
