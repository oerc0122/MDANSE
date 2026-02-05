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

from abc import ABC, abstractmethod
from collections.abc import Iterator, Sequence
from functools import cached_property
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar, overload

import h5py
import numpy as np
from more_itertools import first, last, ilen, iterate, nth, numeric_range

from MDANSE.Framework.QVectors.IQVectors import WIDTH_NONZERO_LIMIT, IQVectors

if TYPE_CHECKING:
    import numpy.typing as npt


class _QVectorContainer(ABC):
    @property
    @abstractmethod
    def source(self) -> Any: ...

    @property
    @abstractmethod
    def shells(self) -> npt.NDArray[np.floating]: ...

    @property
    @abstractmethod
    def params(self) -> dict[str, Any]:
        """Get the Q-Vector generation parameters."""

    @property
    def n_shells(self) -> int:
        """Number of Q shells.

        Returns
        -------
        int
            Number of shells.
        """
        return len(self.shells)

    @property
    def valid_shell_ind(self) -> Iterator[int]:
        return (
            shell
            for shell in range(self.n_shells)
            if self._get_shell(shell) is not None
        )

    @property
    def valid_shells(self) -> Iterator[float]:
        return (
            shell
            for ind, shell in enumerate(self.shells)
            if self._get_shell(ind) is not None
        )

    @property
    def n_valid_shells(self) -> int:
        return ilen(self.valid_shell_ind)

    def q_at_index(self, ind: int) -> npt.NDArray[np.floating]:
        """The Nth shell's Q-vectors from the list.

        Parameters
        ----------
        ind : int
            Q-vector index.

        Returns
        -------
        npt.NDArray[np.floating]
            Q-vectors in Nth shell.
        """
        if (qvec := self._get_shell(ind, "q_vectors")) is not None:
            return qvec[:]

        return np.empty((0,), dtype=float)

    def weight_at_index(self, ind: int) -> npt.NDArray[np.floating]:
        """The Nth shell's weights from the list.

        Parameters
        ----------
        ind : int
            Q-vector index.

        Returns
        -------
        npt.NDArray[np.floating]
            Weights in Nth shell.
        """
        if (weights := self._get_shell(ind, "weights")) is not None:
            return weights[:]

        return np.ones_like(self.q_at_index(ind))

    @property
    def q_vectors(self) -> Iterator[npt.NDArray[np.floating]]:
        """Set of all QVectors."""
        return (self.q_at_index(i) for i in range(self.n_shells))

    @property
    def weights(self) -> Iterator[npt.NDArray[np.floating]]:
        """Set of all weights."""
        return (self.weight_at_index(i) for i in range(self.n_shells))

    @property
    @abstractmethod
    def available_vectors(self) -> npt.NDArray[np.integer]: ...

    @property
    @abstractmethod
    def n_found(self) -> npt.NDArray[np.integer]: ...

    @property
    @abstractmethod
    def n_used(self) -> npt.NDArray[np.integer]: ...

    @property
    def filename(self) -> None:
        """Dummy value for interface."""
        return None

    @overload
    def _get_shell(
        self, ind: int, prop: str = ..., *, strict: Literal[False] = False
    ) -> npt.NDArray[np.floating] | None: ...
    @overload
    def _get_shell(
        self, ind: int, prop: str = ..., *, strict: Literal[True] = ...
    ) -> npt.NDArray[np.floating]: ...
    @abstractmethod
    def _get_shell(
        self, ind: int, prop: str = "", *, strict: bool = False
    ) -> npt.NDArray[np.floating] | None: ...

    @overload
    def _get_all_shells(
        self, prop: str = ..., *, strict: Literal[False] = False
    ) -> Iterator[npt.NDArray[np.floating] | None]: ...
    @overload
    def _get_all_shells(
        self, prop: str = ..., *, strict: Literal[True] = ...
    ) -> Iterator[npt.NDArray[np.floating]]: ...
    def _get_all_shells(
        self, prop: str = "", *, strict: bool = False
    ) -> Iterator[npt.NDArray[np.floating] | None]:
        for index in self.valid_shell_ind:
            yield self._get_shell(index, prop, strict=strict)


class _QVecFile(_QVectorContainer):
    """Wrapper class for getting QVector info from a file.


    Parameters
    ----------
    source : h5py.File
        File to load qvectors from.
    main_dataset : str
        Path to main dataset.
    """

    def __init__(self, source: h5py.File, main_dataset: str):
        self._source = source
        self._main_dataset = main_dataset

    @property
    def source(self) -> h5py.File:
        """Data source.

        Returns
        -------
        h5py.File
            Data source.
        """
        return self._source

    @property
    def shells(self) -> npt.NDArray[np.floating]:
        """Q shells.

        Returns
        -------
        npt.NDArray[np.floating]
            |q| of Q Shells.
        """
        return self._data["q"][:]

    @property
    def params(self) -> dict[str, Any]:
        return self.source["inputs/q_vectors"][1]

    @property
    def n_found(self) -> npt.NDArray[np.integer]:
        if "n_q_found" in self._data:
            return self._data["n_q_found"]

        return np.array(
            [sum(dat) for dat in self._get_all_shells("weights", strict=True)]
        )

    @property
    def n_used(self) -> npt.NDArray[np.integer]:
        if "n_q_vectors" in self._data:
            return self._data["n_q_vectors"]

        return np.array(
            [len(dat) for dat in self._get_all_shells("weights", strict=True)]
        )

    @overload
    def _get_shell(
        self, ind: int, prop: str = "", *, strict: Literal[False] = False
    ) -> h5py.Group | h5py.Dataset | None: ...
    @overload
    def _get_shell(
        self, ind: int, prop: str = "", *, strict: Literal[True] = ...
    ) -> h5py.Group | h5py.Dataset: ...
    def _get_shell(self, ind, prop="", *, strict=False):
        """Get a shell property.

        Parameters
        ----------
        ind : int
            Shell to get property from.
        prop : str, optional
            Dataset to get, if not provided return the Group.
        strict : bool
            Whether to raise KeyError if `prop` not present.

        Returns
        -------
        h5py.Group | h5py.Dataset | None
            Requested data or `None` if not present.
        """
        if strict:
            return self._data[f"shell_{ind}/{prop}"]

        return self._data.get(f"shell_{ind}/{prop}")

    @property
    def available_vectors(self) -> npt.NDArray[np.integer]:
        """List of number of vectors in each Q-shell.

        Returns
        -------
        npt.NDArray[int]
            Number of Q-vectors in shell.
        """
        return np.array([len(self.q_at_index(n)) for n in range(self.n_shells)])

    @property
    def filename(self) -> str:
        """Filename of source.

        Returns
        -------
        str
            Source filename.
        """
        return self.source.filename

    @property
    def _data(self) -> h5py.Dataset:
        """Main dataset.

        Returns
        -------
        h5py.Dataset
            Main dataset.
        """
        return self.source[self._main_dataset]


class _QVecObj(_QVectorContainer):
    """Wrapper class for getting QVector info from a IQVectors.

    Parameters
    ----------
    source : IQVectors
        FIXME: Add docs.
    """

    def __init__(self, source: IQVectors):
        self._source = source

    @property
    def source(self) -> IQVectors:
        """Data source.

        Returns
        -------
        IQVectors
            Data source.
        """
        return self._source

    @property
    def shells(self) -> npt.NDArray[np.floating]:
        """Shells in source.

        Returns
        -------
        npt.NDArray[np.floating]
            Sequnce of shells.
        """
        return np.array(list(self.source["q_vectors"].keys()))

    @property
    def params(self) -> dict[str, Any]:
        return self.source.configuration

    @property
    def available_vectors(self) -> npt.NDArray[np.integer]:
        """Number of vectors in each shell.

        Returns
        -------
        npt.NDArray[int]
            Number of vectors in each shell.
        """
        return np.array(
            [
                shell["n_q_vectors"] if shell is not None else 0
                for shell in self.source["q_vectors"].values()
            ]
        )

    @property
    def n_shells(self) -> int:
        """Number of shells.

        Returns
        -------
        int
            Number of shells.
        """
        return len(self.shells)

    @property
    def n_found(self) -> npt.NDArray[np.integer]:
        return np.array(self._get_all_shells("n_q_found"))

    @property
    def n_used(self) -> npt.NDArray[np.integer]:
        return np.array(self._get_all_shells("n_q_vectors"))

    @overload
    def _get_shell(
        self, ind: int, prop: str = "", *, strict: Literal[False] = False
    ) -> npt.NDArray[np.floating] | None: ...
    @overload
    def _get_shell(
        self, ind: int, prop: str = "", *, strict: Literal[True] = ...
    ) -> npt.NDArray[np.floating]: ...
    def _get_shell(self, ind, prop="", *, strict=False):
        """Get a shell property.

        Parameters
        ----------
        ind : int
            Shell to get property from.
        prop : str, optional
            Dataset to get, if not provided return the Group.
        strict : bool
            Whether to raise KeyError if `prop` not present.

        Returns
        -------
        h5py.Group | h5py.Dataset | None
            Requested data or `None` if not present.
        """
        shell = nth(self.source["q_vectors"].values(), ind)
        if strict:
            if shell is None:
                raise KeyError(f"Shell {ind} does not exist")

            return shell[prop]

        shell = shell or {}
        return shell.get(prop)


T = TypeVar("T", IQVectors, h5py.File)


class QVectorStats(Generic[T]):
    """Get statistics on a Q-vector set.

    Parameters
    ----------
    source : h5py.File | IQVectors
        Source for Q-vectors.
    main_dataset : str
        Location within file where Q-vectors are stored.

    Raises
    ------
    ValueError
        Invalid source provided.
    """

    MAX_BINS_PER_SHELL = 20
    MIN_BINS_PER_SHELL = 1
    MAX_BINS_PER_PLOT = 180

    def __init__(
        self,
        source: T,
        main_dataset: str = "vector_generator",
    ):

        match source:
            case h5py.File():
                self.data = _QVecFile(source, main_dataset)
            case IQVectors():
                self.data = _QVecObj(source)
            case _:
                raise ValueError(
                    f"Do not know how to process {type(source).__name__} as {type(self).__name__}."
                )

    @property
    def data(self) -> _QVectorContainer:
        return self._data

    @data.setter
    def data(self, value: _QVectorContainer) -> None:
        self._data = value
        if hasattr(self, "mod_q_per_shell"):
            del self.mod_q_per_shell

    @property
    def qvector_params(self) -> dict[str, Any]:
        """Get the Q-Vector generation parameters."""
        return self.data.params

    @property
    def weights(self) -> Iterator[npt.NDArray[np.floating]]:
        """Set of all weights."""
        return self.data.weights

    @property
    def shells(self) -> npt.NDArray[np.floating]:
        """Requested shells."""
        return self.data.shells

    @property
    def valid_shells(self) -> Iterator[float]:
        """Valid shells."""
        return self.data.valid_shells

    @property
    def valid_shell_ind(self) -> Iterator[int]:
        """Valid shells."""
        return self.data.valid_shell_ind

    @property
    def filename(self) -> str | None:
        """Data filename if h5 file, else None."""
        return self.data.filename

    @property
    def n_shells(self) -> int:
        """Number of shells."""
        return self.data.n_shells

    @cached_property
    def mod_q_per_shell(self) -> list[npt.NDArray[np.floating]]:
        """Array of q-vector lengths."""
        return [
            np.linalg.norm(self.data.q_at_index(n), axis=0)
            for n in range(self.n_shells)
        ]

    @property
    def available_vectors(self) -> npt.NDArray[np.integer]:
        """Number of vectors in each Q shell."""
        return self.data.available_vectors

    @property
    def mean_q(self) -> npt.NDArray[np.floating]:
        """Mean length in each Q shell."""
        return np.array(
            [np.mean(qvecs) for qvecs in self.mod_q_per_shell if qvecs.size]
        )

    @property
    def mean_q_yerr(self) -> npt.NDArray[np.floating]:
        """Standard deviation in |q| of each Q shell."""
        return np.array([np.std(qvecs) for qvecs in self.mod_q_per_shell if qvecs.size])

    @property
    def q_diff(self) -> list[npt.NDArray[np.floating]]:
        """Distance between each vector and its shell mean."""
        return [
            mod_q - shell
            for shell, mod_q in zip(self.shells, self.mod_q_per_shell, strict=True)
        ]

    @property
    def qmin(self) -> np.floating:
        """Smallest Q length."""
        ind = first(self.valid_shell_ind)
        return np.min(self.mod_q_per_shell[ind])

    @property
    def qmax(self) -> np.floating:
        """Largest Q length.

        Returns
        -------
        float
            Largest |q|.
        """
        ind = last(self.valid_shell_ind)
        return np.max(self.mod_q_per_shell[ind])

    @property
    def q_step(self) -> float:
        """Average shell width or 0.1 if only one shell."""
        return np.mean(np.abs(np.diff(self.shells))) if self.n_shells > 1 else 0.1

    def _calc_bin_step(self) -> float:
        """Compute bin_step if not provided.

        Returns
        -------
        float
            Binning step.
        """
        bin_step = 0.4 * np.min(self.mean_q_yerr)

        return (
            0.2 * self.q_step
            if abs(bin_step) < 1e-09
            else max(bin_step, 0.05 * self.q_step)
        )

    @staticmethod
    def _get_bin_width(
        n_segments: int, start: float, end: float, step_size: float, width: float | None
    ) -> tuple[float, float]:
        """Return the bin width based on shell width and shell separation."""
        match width:
            case None:
                return step_size / n_segments, step_size
            case width if width <= step_size / 2:
                return step_size / n_segments, step_size
            case width if np.isclose(start, end):
                return width / n_segments, width
            case width:
                return step_size / n_segments, width

    @classmethod
    def _find_binning(
        cls, start: float, end: float, bin_width: float, peak_width: float
    ) -> tuple[float, float, float]:
        """Find ideal binning for data.

        Parameters
        ----------
        start : float
            Lower limit of the |q| shell range.
        end : float
            Upper limit of the |q| shell range.
        bin_width : float
            Initial guess at bin width.
        peak_width : float
            Peak width.

        Returns
        -------
        tuple[float, float, float]
            Best binning for data.
        """

        def _get_first_last_values(bin_width: float) -> tuple[float, float]:
            """Return the limits of the binning range based on the shell and bin sizes."""
            bins_per_shell: float = np.clip(
                peak_width // bin_width,
                cls.MIN_BINS_PER_SHELL,
                cls.MAX_BINS_PER_SHELL,
            )
            bin_width = peak_width / bins_per_shell

            first_value = start - 0.5 * (bins_per_shell + 1) * bin_width
            last_value = end + 0.5 * (bins_per_shell + 1.01) * bin_width
            return first_value, last_value

        def _trial_iter(bin_width: float) -> Iterator[tuple[float, float, float]]:
            for width in iterate(lambda x: x * 2, bin_width):
                first, last = _get_first_last_values(bin_width)
                yield first, last, width

        prev = prev_n = None
        for trial in _trial_iter(bin_width):
            n = len(numeric_range(*trial))

            if n > cls.MAX_BINS_PER_SHELL:
                return prev or trial
            if n == prev_n:
                return trial

            prev = trial
            prev_n = n

        return trial

    def qvector_binning(self, n_segments: int = 10):
        return self.qvector_binning_from_params(self.qvector_params, n_segments)

    @classmethod
    def qvector_binning_from_params(
        cls,
        qvector_params: dict[str, Any],
        n_segments: int = 10,
    ) -> npt.NDArray[np.floating] | None:
        """Calculate the range of |q| bins from the vector generator parameters.

        Parameters
        ----------
        qvector_params : dict[str, Any]
            Input parameters of any subclass of IQVectors given as a dictionary.
        n_segments : int, optional
            Starting value of the number of histogram bins per vector shell, by default 10

        Returns
        -------
        npt.NDArray[np.floating] | None
            A 1D array of |q| bin limits.
        """

        print(qvector_params)

        if (step_params := qvector_params.get("shells")) is None:
            return None

        if isinstance(step_params, Sequence):
            start, end, step_size = step_params
        else:
            start, end, step_size = (
                step_params[prop] for prop in ("first", "last", "step")
            )

        width = qvector_params.get("width", {}).get("value")
        return cls.qvector_binning_general(start, end, step_size, width, n_segments)

    @classmethod
    def qvector_binning_general(
        cls,
        start: float,
        end: float,
        step_size: float,
        width: float | None,
        n_segments: int,
    ) -> npt.NDArray[np.floating]:
        """Calculate the |q| bin limits based of vector shell limits and steps.

        Parameters
        ----------
        start : float
            Lower limit of the |q| shell range.
        end : float
            Upper limit of the |q| shell range.
        step_size : float
            Step size of the |q| shell range.
        width : float | None
            Width of a single shell in |q| units (1/nm).
        n_segments : int
            Starting value of number of histogram bins per shell.

        Returns
        -------
        npt.NDArray[np.floating]
            A 1D array of |q| histogram bin limits.
        """
        start, end = sorted((start, end))
        if np.isclose(start, end) and (
            width is None or np.isclose(width, 0, atol=WIDTH_NONZERO_LIMIT)
        ):
            return np.array([start - 0.15, start - 0.05, start + 0.05, start + 0.15])

        width = abs(width) if width else width
        step_size = abs(step_size)
        bin_width, peak_width = cls._get_bin_width(
            n_segments, start, end, step_size, width
        )

        first, last, bin_width = cls._find_binning(start, end, bin_width, peak_width)

        common_binning = np.arange(first, last, bin_width)
        offset = np.min(np.abs(start - common_binning)) / bin_width
        common_binning -= (offset - 0.5) * bin_width

        return common_binning

    def histogram(
        self, binning: npt.NDArray[np.floating] | Sequence[float] | float | None = None
    ) -> tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]:
        """Get a histogram of Q-vectors at bin-lengths.

        Parameters
        ----------
        binning : npt.NDArray[np.floating] or float, optional
            Binning for histogram:
                - if float: considered bin_step.
                - if array: considered bin boundaries.
                - if None: compute from shells.

        Returns
        -------
        x-axis : npt.NDArray[np.floating]
            Bins-centres in shells.
        histogram : npt.NDArray[np.floating]
            Histogram of available q-vectors.
        """
        if isinstance(binning, (Sequence, np.ndarray)):
            common_bins = binning
        else:
            bin_step = float(binning or self._calc_bin_step())

            minim = float(max(0.0, min(self.qmin, self.shells[0])))
            maxim = float(max(self.qmax, self.shells[-1]) + 1.1 * bin_step)
            common_bins = np.arange(minim, maxim, bin_step)

        stacked_histograms = np.vstack([
            np.histogram(qmods, common_bins)[0] for qmods in self.mod_q_per_shell
        ])
        xvals = common_bins[1:] - np.diff(common_bins) / 2

        print("Boop")
        return xvals, stacked_histograms

    def vector_angular_datasets(
        self,
        shell_key: int,
    ) -> Iterator[tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]]:
        """Return a specific q-vector shell as spherical coordinate angle datasets.

        Parameters
        ----------
        source : QVectorsConfigurator
            An instance of the QVectorsConfigurator in the GUI.
        shell_key : float
            The |q| of the vector shell.

        Returns
        -------
        tuple[SingleDataset, SingleDataset]
            Datasets of polar and azimuthal angles.
        """
        q_array = self.data.q_at_index(shell_key)
        q_weights = self.data.weight_at_index(shell_key)
        yield from self.angular_datasets_from_qarray(
            q_array,
            q_weights,
        )

    @staticmethod
    def angular_datasets_from_qarray(
        q_array: npt.NDArray[np.floating],
        weight_array: npt.NDArray[np.floating],
    ) -> Iterator[tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]]:
        """Convert an array of q vectors into two datasets of spherical angles.

        Parameters
        ----------
        q_array : npt.NDArray[np.floating]
            A (3,N) array of reciprocal space vectors.
        weight_array: npt.NDArray[np.floating],
            Array of weights per vector, based the multiplicity of each vector.

        Yields
        ------
        tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]
            Datasets of polar and azimuthal angles.
        """
        inplane_r = np.linalg.norm(q_array[:2, :], axis=0)
        polar_angles = np.arctan2(inplane_r, q_array[2, :])
        azimuthal_angles = np.arctan2(q_array[1, :], q_array[0, :])

        rounding_precision = 3
        for input_angles, normalise, hist_range in (
            (polar_angles, True, (0, np.pi)),
            (azimuthal_angles, False, (-np.pi, np.pi)),
        ):
            angles = np.round(input_angles, rounding_precision)
            counts, bins = np.histogram(angles, weights=weight_array, range=hist_range)
            unique_counts, _ = np.histogram(angles, bins=bins)
            mean_angles = (bins[1:] + bins[:-1]) / 2
            if normalise:
                counts = counts / np.sin(mean_angles)

            yield counts, unique_counts
