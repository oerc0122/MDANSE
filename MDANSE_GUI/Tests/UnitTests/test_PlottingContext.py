from pathlib import Path

import h5py
import numpy as np
import pytest
from MDANSE_GUI.Tabs.Models.PlottingContext import (PlottingContext,
                                                    SingleDataset)
from qtpy import QtCore, QtGui, QtWidgets

file_1d_name = Path(__file__).parent / "super_dos.mda"
file_2d_name = Path(__file__).parent / "disf_argon.mda"


@pytest.fixture
def file_1d():
    data_file = h5py.File(file_1d_name)
    yield data_file
    data_file.close()

@pytest.fixture
def file_2d():
    data_file = h5py.File(file_2d_name)
    yield data_file
    data_file.close()

@pytest.fixture
def data(request, file_1d, file_2d):
    if request.param == "1d":
        return SingleDataset("dos_total", file_1d)
    if request.param == "2d":
        return SingleDataset("f(q,t)_total", file_2d)

    raise FileNotFoundError(f"Unrecognised file {request.param}")


@pytest.mark.parametrize("data, expected", [
    ("1d", (501,)),
    ("2d", (10, 501))], indirect=["data"])
def test_single_dataset(data, expected):
    assert data._data.shape == expected

@pytest.mark.parametrize("data, axes, units, longest_axis", [
    ("1d", ["omega"], ["rad/ps"], ("rad/ps", "omega")),
    ("2d", ["q", "time"], ["1/nm", "ps"], ("ps", "time"))
], indirect=["data"])
def test_available_x_axes(data, axes, units, longest_axis):
    assert data.available_x_axes() == axes
    assert [data._axes_units[axis] for axis in axes] == units
    assert data.longest_axis() == longest_axis

@pytest.mark.parametrize("data, axis, expected", [
    ("1d", "rad/ps", (1, ((), 501))),
    ("2d", "ps", (10, ((0,), 501))),
    ("2d", "1/nm", (501, ((0,), 10))),
], indirect=["data"], ids=["1d", "2d_long_axis", "2d_short_axis"])
def test_curves_vs_axis(data, axis, expected):
    curves = data.curves_vs_axis(axis)
    first_axis, (ind, second_axis) = expected
    assert len(curves) == first_axis
    assert len(curves[ind]) == second_axis
