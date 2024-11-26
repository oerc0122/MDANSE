import os
import tempfile
import pytest

from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Configurators.HDFTrajectoryConfigurator import (
    HDFTrajectoryConfigurator,
)


file_wd = os.path.dirname(os.path.realpath(__file__))

lammps_config = os.path.join(file_wd, "Data", "lammps_test.config")
lammps_lammps = os.path.join(file_wd, "Data", "lammps_test.lammps")
lammps_moly = os.path.join(file_wd, "Data", "structure_moly.lammps")
lammps_custom = os.path.join(file_wd, "Data", "lammps_moly_custom.txt")
lammps_xyz = os.path.join(file_wd, "Data", "lammps_moly_xyz.txt")
lammps_h5md = os.path.join(file_wd, "Data", "lammps_moly_h5md.h5")
vasp_xdatcar = os.path.join(file_wd, "Data", "XDATCAR_version5")
discover_his = os.path.join(file_wd, "Data", "sushi.his")
discover_xtd = os.path.join(file_wd, "Data", "sushi.xtd")
cp2k_pos = os.path.join(file_wd, "Data", "CO2GAS-pos-1.xyz")
cp2k_vel = os.path.join(file_wd, "Data", "CO2GAS-vel-1.xyz")
cp2k_cell = os.path.join(file_wd, "Data", "CO2GAS-1.cell")
hem_cam_pdb = os.path.join(file_wd, "Data", "hem-cam.pdb")
hem_cam_dcd = os.path.join(file_wd, "Data", "hem-cam.dcd")
ase_traj = os.path.join(file_wd, "Data", "Cu_5steps_ASEformat.traj")
xyz_traj = os.path.join(file_wd, "Data", "traj-100K-npt-1000-res.xyz")
dlp_field_v2 = os.path.join(file_wd, "Data", "FIELD_Water")
dlp_history_v2 = os.path.join(file_wd, "Data", "HISTORY_Water")
dlp_field_v4 = os.path.join(file_wd, "Data", "FIELD4")
dlp_history_v4 = os.path.join(file_wd, "Data", "HISTORY4")
apoferritin_dcd = os.path.join(file_wd, "Data", "apoferritin.dcd")
apoferritin_pdb = os.path.join(file_wd, "Data", "apoferritin.pdb")
pbanew_md = os.path.join(file_wd, "Data", "PBAnew.md")
h2o_trj = os.path.join(file_wd, "Data", "H2O.trj")
h2o_xtd = os.path.join(file_wd, "Data", "H2O.xtd")
md_pdb = os.path.join(file_wd, "Data", "md.pdb")
md_xtc = os.path.join(file_wd, "Data", "md.xtc")


@pytest.mark.parametrize("compression", ["none", "gzip", "lzf"])
def test_lammps_mdt_conversion_file_exists_and_loads_up_successfully(compression):
    temp_name = tempfile.mktemp()

    parameters = {}
    parameters["config_file"] = lammps_config
    parameters["mass_tolerance"] = 0.05
    parameters["n_steps"] = 0
    parameters["output_files"] = (temp_name, 64, 128, compression, "INFO")
    parameters["smart_mass_association"] = True
    parameters["time_step"] = 1.0
    parameters["trajectory_file"] = lammps_lammps

    lammps = Converter.create("LAMMPS")
    lammps.run(parameters, status=True)

    HDFTrajectoryConfigurator("trajectory").configure(temp_name + ".mdt")

    assert os.path.exists(temp_name + ".mdt")
    assert os.path.isfile(temp_name + ".mdt")
    os.remove(temp_name + ".mdt")
    assert os.path.exists(temp_name + ".log")
    assert os.path.isfile(temp_name + ".log")
    os.remove(temp_name + ".log")
