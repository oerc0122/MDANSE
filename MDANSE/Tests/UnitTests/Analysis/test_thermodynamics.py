import sys
import tempfile
import os
from os import path
import pytest

import numpy as np
import h5py
from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from MDANSE.Framework.Jobs.IJob import IJob


sys.setrecursionlimit(100000)
short_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "Data",
    "short_trajectory_after_changes.mdt",
)
com_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "Data",
    "com_trajectory.mdt",
)


@pytest.mark.parametrize(
    "trajectory_name,interp_order",
    [(short_traj, 1), (short_traj, 3), (com_traj, 1), (com_traj, 3)],
)
def test_temperature(trajectory_name, interp_order):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["frames"] = (0, 10, 1)
    parameters["interpolation_order"] = interp_order
    parameters["output_files"] = (temp_name, ("MDAFormat",), "INFO")
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    temp = IJob.create("Temperature")
    temp.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


@pytest.mark.parametrize(
    "trajectory_name,interp_order",
    [(short_traj, 1), (short_traj, 3), (com_traj, 1), (com_traj, 3)],
)
def test_temperature_nonzero(trajectory_name, interp_order):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["frames"] = (0, 10, 1)
    parameters["interpolation_order"] = interp_order
    parameters["output_files"] = (temp_name, ("MDAFormat",), "INFO")
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    temp = IJob.create("Temperature")
    temp.run(parameters, status=True)
    with h5py.File(temp_name + ".mda") as results:
        print(results.keys())
        temperature = np.array(results["/temperature"])
    os.remove(temp_name + ".mda")
    assert np.all(temperature > 0.0)


@pytest.mark.parametrize(
    "trajectory_name,output_format",
    [
        (short_traj, "MDAFormat"),
        (com_traj, "MDAFormat"),
        (short_traj, "TextFormat"),
        (com_traj, "TextFormat"),
    ],
)
def test_density(trajectory_name, output_format):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["frames"] = (0, 10, 1)
    parameters["output_files"] = (temp_name, (output_format,), "INFO")
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    den = IJob.create("Density")
    den.run(parameters, status=True)
    if output_format == "MDAFormat":
        assert path.exists(temp_name + ".mda")
        assert path.isfile(temp_name + ".mda")
        os.remove(temp_name + ".mda")
    elif output_format == "TextFormat":
        assert path.exists(temp_name + "_text.tar")
        assert path.isfile(temp_name + "_text.tar")
        os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")
