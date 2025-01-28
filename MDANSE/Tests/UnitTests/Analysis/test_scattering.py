import sys
import tempfile
import os
from os import path

import numpy as np
import h5py
import pytest

from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from MDANSE.Framework.Jobs.IJob import IJob


sys.setrecursionlimit(100000)
short_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "Data",
    "short_trajectory_after_changes.mdt",
)
result_dir = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "Results",
)


@pytest.fixture(scope="module")
def trajectory():
    trajectory = HDFTrajectoryInputData(short_traj)
    yield trajectory


@pytest.fixture(scope="module")
def qvector_spherical_lattice(trajectory):
    return (
        "SphericalLatticeQVectors",
        {"seed": 0, "shells": (5.0, 36, 10.0), "n_vectors": 10, "width": 9.0},
    )


@pytest.fixture(scope="module")
def dcsf():
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat",), "INFO")
    parameters["q_vectors"] = (
        "SphericalLatticeQVectors",
        {"seed": 0, "shells": (5.0, 36, 10.0), "n_vectors": 10, "width": 9.0},
    )
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_coherent"
    dcsf = IJob.create("DynamicCoherentStructureFactor")
    dcsf.run(parameters, status=True)
    yield temp_name + ".mda"
    os.remove(temp_name + ".mda")


@pytest.fixture(scope="module")
def disf():
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat",), "INFO")
    parameters["q_vectors"] = (
        "SphericalLatticeQVectors",
        {"seed": 0, "shells": (5.0, 36, 10.0), "n_vectors": 10, "width": 9.0},
    )
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_incoherent2"
    disf = IJob.create("DynamicIncoherentStructureFactor")
    disf.run(parameters, status=True)
    yield temp_name + ".mda"
    os.remove(temp_name + ".mda")


def test_dcsf(trajectory, qvector_spherical_lattice):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    parameters["q_vectors"] = qvector_spherical_lattice
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_coherent"
    dcsf = IJob.create("DynamicCoherentStructureFactor")
    dcsf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


def test_ccf(trajectory, qvector_spherical_lattice):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    parameters["q_vectors"] = qvector_spherical_lattice
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "equal"
    dcsf = IJob.create("CurrentCorrelationFunction")
    dcsf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


def test_output_axis_preview(qvector_spherical_lattice):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    parameters["q_vectors"] = qvector_spherical_lattice
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_coherent"
    dcsf = IJob.create("DynamicCoherentStructureFactor")
    dcsf.setup(parameters)
    axes = dcsf.preview_output_axis()
    print(axes)
    assert len(axes) == 3  # two configurators return valid arrays


def test_disf(trajectory, qvector_spherical_lattice):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    parameters["q_vectors"] = qvector_spherical_lattice
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_incoherent2"
    disf = IJob.create("DynamicIncoherentStructureFactor")
    disf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


def test_eisf(trajectory, qvector_spherical_lattice):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1)
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    parameters["q_vectors"] = qvector_spherical_lattice
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_incoherent"
    eisf = IJob.create("ElasticIncoherentStructureFactor")
    eisf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


def test_gdisf(trajectory):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    parameters["q_shells"] = (2.0, 12.2, 2.0)
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_incoherent2"
    gdisf = IJob.create("GaussianDynamicIncoherentStructureFactor")
    gdisf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    result_file = os.path.join(result_dir, f"gdisf.mda")

    with h5py.File(temp_name + ".mda") as actual,  h5py.File(result_file) as desired:
        np.testing.assert_array_almost_equal(actual["/f(q,t)_Cu"], desired["/f(q,t)_Cu"])
        np.testing.assert_array_almost_equal(actual["/f(q,t)_S"], desired["/f(q,t)_S"])
        np.testing.assert_array_almost_equal(actual["/f(q,t)_Sb"], desired["/f(q,t)_Sb"])
        np.testing.assert_array_almost_equal(actual["/f(q,t)_total"], desired["/f(q,t)_total"])
        np.testing.assert_array_almost_equal(actual["/s(q,f)_Cu"], desired["/s(q,f)_Cu"])
        np.testing.assert_array_almost_equal(actual["/s(q,f)_S"], desired["/s(q,f)_S"])
        np.testing.assert_array_almost_equal(actual["/s(q,f)_Sb"], desired["/s(q,f)_Sb"])
        np.testing.assert_array_almost_equal(actual["/s(q,f)_total"], desired["/s(q,f)_total"])

    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


def test_ndtsf(disf, dcsf, qvector_spherical_lattice):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["disf_input_file"] = disf
    parameters["dcsf_input_file"] = dcsf
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    ndtsf = IJob.create("NeutronDynamicTotalStructureFactor")
    ndtsf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")


def test_ssfsf(disf):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["sample_inc"] = disf
    parameters["running_mode"] = ("single-core",)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat", "TextFormat"), "INFO")
    ndtsf = IJob.create("StructureFactorFromScatteringFunction")
    ndtsf.run(parameters, status=True)
    assert path.exists(temp_name + ".mda")
    assert path.isfile(temp_name + ".mda")
    os.remove(temp_name + ".mda")
    assert path.exists(temp_name + "_text.tar")
    assert path.isfile(temp_name + "_text.tar")
    os.remove(temp_name + "_text.tar")
    assert path.exists(temp_name + ".log")
    assert path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")
