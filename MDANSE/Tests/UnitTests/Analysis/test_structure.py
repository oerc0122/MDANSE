from pathlib import Path

import numpy as np
import h5py
import pytest

from MDANSE.Framework.Jobs.IJob import IJob

BASE_PATH = Path(__file__).parent.parent.resolve()
short_traj = BASE_PATH / "Converted" / "short_trajectory_after_changes.mdt"
mdmc_traj = BASE_PATH / "Converted" / "Ar_mdmc_h5md.h5"
com_traj = BASE_PATH / "Converted" / "com_trajectory.mdt"
result_dir = BASE_PATH / "Results"


################################################################
# Job parameters                                               #
################################################################
@pytest.fixture(scope="function")
def parameters():
    parameters = {}
    # parameters['atom_selection'] = None
    # parameters['atom_transmutation'] = None
    # parameters['frames'] = (0, 1000, 1)
    # Set and overridden in function
    # parameters["trajectory"] = short_traj
    # parameters["running_mode"] = ("multicore", -4)
    parameters["q_vectors"] = (
        "SphericalLatticeQVectors",
        {
            "seed": 0,
            "shells": (0, 5.0, 0.5),
            "n_vectors": 100,
            "width": 0.5,
        },
    )
    parameters["q_values"] = (0.0, 10.0, 0.1)
    parameters["r_values"] = (0.0, 0.9, 0.1)
    parameters["per_axis"] = False
    parameters["reference_direction"] = (0, 0, 1)
    parameters["instrument_resolution"] = ("Gaussian", {"sigma": 1.0, "mu": 0.0})
    parameters["interpolation_order"] = "3rd order"
    parameters["projection"] = None
    parameters["normalize"] = True
    parameters["grouping_level"] = "atom"
    parameters["weights"] = "equal"
    return parameters


@pytest.mark.parametrize("traj_info", (("short_traj", short_traj),
                                       ("mdmc_traj", mdmc_traj),
                                       ("com_traj", com_traj)))
@pytest.mark.parametrize(
    "job_info",
    (
        ("RadiusOfGyration", {"rog"}),
        ("DensityProfile", {"dp"}),
        ("MolecularTrace", {
            "molecular_trace",
            "x_position",
            "y_position",
            "z_position"
        }),
        ("Eccentricity", {"eccentricity"}),

        pytest.param(("SolventAccessibleSurface", {"sas"}), marks=pytest.mark.long),
        pytest.param(("RootMeanSquareDeviation", {"rmsd"}), marks=pytest.mark.long),
        pytest.param(("RootMeanSquareFluctuation", {"rmsf"}), marks=pytest.mark.long),
        pytest.param(("Voronoi", {"mean_volume",
                                  "neighbourhood_histogram"}), marks=pytest.mark.long),
        pytest.param(("CoordinationNumber", {"cn"}), marks=pytest.mark.long),
        pytest.param(("PairDistributionFunction", {"pdf",
                                                   "rdf",
                                                   "tcf"}), marks=pytest.mark.long),
        pytest.param(("StaticStructureFactor", {"ssf"}), marks=pytest.mark.long),
        pytest.param(("XRayStaticStructureFactor", {"xssf"}), marks=pytest.mark.long),

    )
)
@pytest.mark.parametrize("running_mode", (("single-core", 1), ("multicore", -4)))
@pytest.mark.parametrize("output_format", ("MDAFormat", "TextFormat"))
def test_structure_analysis(
        request, tmp_path, parameters, traj_info, job_info, running_mode, output_format
):

    if (
            any(mark.name == "long" for mark in request.node.own_markers) and
            (running_mode[0] == "multicore" or output_format == "TextFormat")
    ):
        pytest.skip("Test is long.")

    temp_name = tmp_path / "temp"
    log_file = tmp_path / "temp.log"
    parameters["trajectory"] = str(traj_info[1])
    parameters["running_mode"] = running_mode
    parameters["output_files"] = (str(temp_name), (output_format,), "INFO")
    job = IJob.create(job_info[0])
    job.run(parameters, status=True)

    if output_format == "MDAFormat":
        out_file = tmp_path / "temp.mda"
        assert out_file.exists() and out_file.is_file()

        result_file = result_dir / f"structure_analysis_{traj_info[0]}_{job_info[0]}.mda"

        with h5py.File(out_file) as actual, h5py.File(result_file) as desired:
            keys = [key for key in desired if any(key.startswith(typ)
                                                  for typ in job_info[1])]
            for key in keys:
                np.testing.assert_array_almost_equal(actual[f"/{key}"], desired[f"/{key}"])

    elif output_format == "TextFormat":
        out_file = tmp_path / "temp_text.tar"
        assert out_file.exists() and out_file.is_file()

    assert log_file.exists() and log_file.is_file()
