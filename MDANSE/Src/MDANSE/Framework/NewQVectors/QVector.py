import itertools
from abc import ABC, abstractmethod
from collections.abc import Generator, Iterator
from dataclasses import dataclass
from typing import Optional
from functools import singledispatchmethod

import numpy as np
from MDANSE.MolecularDynamics.UnitCell import UnitCell
from numpy.typing import NDArray


@dataclass(slots=True)
class QVectorData:
    """Data relating to a single Q-Vector

    Returned from a QVectorGenerator"""

    q: NDArray[float]
    mod_q: float
    hkl: NDArray[int]
    hkl_exact: NDArray[float]

    def __init__(self, q, lattice: Optional[UnitCell] = None):
        self.q = np.array(q)
        self.mod_q = np.linalg.norm(q)

        if lattice is not None:
            self.hkl_exact = np.dot(lattice.direct, q)
            self.hkl = np.rint(self.hkl_exact)
        else:
            self.hkl_exact = None
            self.hkl = None

    def __eq__(self, other):
        return np.allclose(self.q, other.q)

class QVectorGenerator(ABC):
    """Abstract type for generation of Q-Vectors."""

    def __init__(self, lattice: Optional[UnitCell] = None, **kwargs):
        self.lattice = lattice
        self._ind = 0

    def __iter__(self):
        return self.generate()

    @singledispatchmethod
    def __getitem__(self, s):
        """Get a subset of the Qvectors

        Parameters
        ----------
        s : Union[slice, int]
            If :class:`slice`, get that subset.
            If :class:`int`, get the first N vectors.

        Examples
        --------
        >>> from MDANSE.Framework.NewQVectors.LinearQVectors import LinearQVectors
        >>> vec = LinearQVectors([1, 0, 0])
        >>> [qvec.mod_q for qvec in vec[3]]
        [1., 2., 3.]
        >>> [qvec.mod_q for qvec in vec[4:6]]
        [4., 5.]
        """
        raise NotImplementedError(f"Cannot slice by {type(s).__name__}")

    @__getitem__.register(slice)
    def _(self, s):
        return itertools.islice(self.generate(), s.start, s.stop, s.step)

    @__getitem__.register(int)
    def _(self, n):
        return self.generate_n(n)

    def generate_n(self, n: int) -> Iterator[QVectorData]:
        """Generate a set of N Q-Vectors with default arguments.

        Parameters
        ----------
        n : int
            Number of Q-Vectors to generate.

        Examples
        --------
        FIXME: Add docs.
        """
        return itertools.islice(self.generate(), n)

    @abstractmethod
    def generate(self, n=None, lattice: UnitCell = None) -> Generator[QVectorData, int, None]:
        raise NotImplementedError

    @abstractmethod
    def reset(self, value: int = 0):
        self._ind = value

    @property
    def index(self):
        return self._ind

    @index.setter
    def index(self, value):
        self.reset(value)
