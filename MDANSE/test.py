from pathlib import Path
from pprint import pprint

from MDANSE.Framework.Converters import ASE, CASTEP, VASP, Forcite
from MDANSE.Framework.Jobs.Density import Density
from MDANSE.Framework.Parsers import TrjFile, XDATCARFile, XYZFile, binary_file_reader

DATA_DIR = Path("~/MDANSE/MDANSE/Tests/UnitTests/Data/").expanduser()

V = VASP()
print(V)
V.fold = True
V.trajectory_file = DATA_DIR / "XDATCAR_version5"
V.output_files.path = "./doom_file.mdt"
V.time_step = 1.0

V.run(prog_bar=True)

# X = CASTEP()

# X.fold = False
# X.trajectory_file = DATA_DIR / "PBAnew.md"
# X.output_files.path = "./doom_file.mdt"
# X.time_unit = "ps"
# X.time_step = 10.0


# X.run()

# Y = ASE()

# Y.fold = False
# Y.trajectory_file = DATA_DIR / "Cu_5steps_ASEformat.traj"
# Y.output_files.path = "./final_doom_file.mdt"
# Y.time_unit = "ps"
# Y.time_step = 10.0

# Y.run()

# Z = Density()
# Z.trajectory = X.output_files.path
# Z.frames = "all"
# Z.output_files.path = "./apm"

# Z.run()

# Q = XYZFile(DATA_DIR / "SrTiO3_MD-pos-1.xyz")
# print(Q)
# print(Q.time_step)
# for i in Q.frames:
#     print(i)

# X = Forcite()

# X.fold = False
# X.xtd_file = DATA_DIR / "H2O.xtd"
# X.trj_file = DATA_DIR / "H2O.trj"
# X.output_files.path = "./doom_file.mdt"
# X.time_unit = "ps"
# X.time_step = 10.0

# X.run()
