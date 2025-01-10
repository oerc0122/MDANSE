import sys
import tempfile
import os
from os import path
import pytest

import numpy as np

from MDANSE.Framework.InputData.HDFTrajectoryInputData import HDFTrajectoryInputData
from MDANSE.Framework.QVectors.IQVectors import IQVectors
from MDANSE.Framework.Jobs.IJob import IJob


sys.setrecursionlimit(100000)
short_traj = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    "..",
    "Data",
    "short_trajectory_after_changes.mdt",
)


@pytest.fixture(scope="module")
def trajectory():
    trajectory = HDFTrajectoryInputData(short_traj)
    yield trajectory


def test_qvector_to_hkl_conversion(trajectory):
    for qvector_generator in IQVectors.indirect_subclasses():
        instance = IQVectors.create(qvector_generator, trajectory.chemical_system)
        instance.setup({"shells": (5.0, 50.0, 10.0)})
        unit_cell = trajectory.trajectory.unit_cell(0)
        instance.generate()
        try:
            instance._configuration["shells"]
        except KeyError:
            print(f"{qvector_generator} has no shells")
            continue
        else:
            print(f"Calculating q vectors for {qvector_generator}")
        for q in instance._configuration["shells"]["value"][:2]:
            try:
                original_qvectors = instance._configuration["q_vectors"][q]["q_vectors"]
            except KeyError:
                continue
            print("q_vectors", original_qvectors)
            if len(original_qvectors) == 0:
                continue
            hkls = instance.qvectors_to_hkl(original_qvectors, unit_cell)
            print("hkls", hkls)
            recalculated_qvectors = instance.hkl_to_qvectors(hkls, unit_cell)
            assert np.allclose(original_qvectors, recalculated_qvectors)



def test_disf(trajectory):
    temp_name = tempfile.mktemp()
    parameters = {}
    parameters["atom_selection"] = None
    parameters["atom_transmutation"] = None
    parameters["frames"] = (0, 10, 1, 5)
    parameters["instrument_resolution"] = ("Ideal", {})
    parameters["output_files"] = (temp_name, ("MDAFormat",), "INFO")
    parameters["running_mode"] = ("single-core",)
    parameters["trajectory"] = short_traj
    parameters["weights"] = "b_incoherent2"
    for qvector_generator in IQVectors.indirect_subclasses():
        instance = IQVectors.create(qvector_generator, trajectory.chemical_system)
        qvector_defaults = {
            name: value[1]["default"] for name, value in instance.settings.items()
        }
        if len(qvector_defaults) < 1:
            continue
        print(qvector_generator)
        print(qvector_defaults)
        parameters["q_vectors"] = (qvector_generator, qvector_defaults)
        disf = IJob.create("DynamicIncoherentStructureFactor")
        disf.run(parameters, status=True)
        assert path.exists(temp_name + ".mda")
        assert path.isfile(temp_name + ".mda")
        os.remove(temp_name + ".mda")
        assert path.exists(temp_name + ".log")
        assert path.isfile(temp_name + ".log")
        os.remove(temp_name + ".log")
