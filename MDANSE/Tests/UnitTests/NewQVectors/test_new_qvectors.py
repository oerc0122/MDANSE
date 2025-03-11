import pytest

from MDANSE.Framework.NewQVectors import QVectorGenerator


@pytest.mark.parametrize("instance, expected", [
    ("PathQVectors", "PathLinearQVectors"),
    ("PathLinearQVectors", "PathLinearQVectors"),
    ("PathSegmentQVectors", "PathSegmentQVectors"),
    ("DispersionQVectors", "PathSegmentQVectors"),
    ("ApproximateDispersionQVectors", "PathSegmentQVectors"),
    ("GridQVectors", "GridQVectors"),
    ("GeneratorQVectors", "GeneratorQVectors"),
    ("ListQVectors", "ListQVectors"),
    ("LatticeSphericalQVectors", "LatticeSphericalQVectors"),
    ("SphericalQVectors", "SphericalQVectors"),
])
def test_factory(instance, expected):
    assert QVectorGenerator.get(instance).__name__ == expected
