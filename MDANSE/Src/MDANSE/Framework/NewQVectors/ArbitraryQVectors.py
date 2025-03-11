from typing import Optional, Tuple

import numpy as np
from MDANSE.Framework.NewQVectors.QVector import (
    QVecGen,
    QVecGeneratorProtocol,
    QVectorGenerator,
)
from MDANSE.MolecularDynamics.UnitCell import UnitCell
from numpy.typing import ArrayLike


@QVectorGenerator.register("ListQVectors")
class ListQVectors(QVectorGenerator):
    """Return Q vectors from a provided list of vectors.

    Parameters
    ----------
    qvectors : ArrayLike
        Sequence of vectors to return.
    hkl : bool
        Whether vectors are provided in reciprocal lattice units.
    lattice : Optional[UnitCell]
        Lattice to generate within.

    Raises
    ------
    ValueError
        If hkl and no lattice provided.

    Examples
    --------
    >>> qvec = ListQVectors([[1, 2, 3], [2, 0, 5]])
    >>> for i in qvec:
    ...     print(i.q)
    [1. 2. 3.]
    [2. 0. 5.]
    """

    def __init__(
        self, qvectors: ArrayLike, *, lattice: Optional[UnitCell] = None, **kwargs
    ):
        super().__init__(lattice=lattice, **kwargs)
        self.qvectors = np.array(qvectors)

        if self.qvectors.shape != (len(self.qvectors), 3):
            raise ValueError(
                f"`qvectors` ({qvectors.shape}) must be an (N, 3) sequence."
            )

    def _generate(self, *, lattice: Optional[UnitCell] = None) -> QVecGen:
        lattice = lattice if lattice is not None else self.lattice
        constructor = self.qvec_gen

        while self._ind < len(self):
            new_ind = yield constructor(self.qvectors[self._ind], lattice)

            self._ind += 1
            if new_ind is not None:
                self.reset(new_ind)

    def __len__(self) -> int:
        return len(self.qvectors)

    def reset(self, value: int = 0):
        super().reset(value)


@QVectorGenerator.register("GeneratorQVectors")
class GeneratorQVectors(QVectorGenerator):
    """Return Q Vectors as generated from a generator function.

    Parameters
    ----------
    qvectors : QVecGen | Generator[ArrayLike, int, None]
        Generator which returns QVectorData.
    lattice : Optional[UnitCell]
        Lattice in which to generate Q Vectors.
    """

    def __init__(
        self,
        qvectors: QVecGeneratorProtocol,
        *args,
        lattice: Optional[UnitCell] = None,
        returns_qvec_data: bool = True,
        **kwargs,
    ):
        super().__init__(lattice=lattice, **kwargs)

        self.generator = qvectors
        self.args = args
        self.returns_qvec_data = returns_qvec_data

    def _generate(self, *, lattice: Optional[UnitCell] = None) -> QVecGen:
        lattice = lattice if lattice is not None else self.lattice
        constructor = self.qvec_gen

        for qvec in self.generator(*self.args, lattice):
            yield qvec if self.returns_qvec_data else constructor(qvec, lattice)

    def reset(self, val: None = None):
        raise NotImplementedError("Cannot reset this generator")
