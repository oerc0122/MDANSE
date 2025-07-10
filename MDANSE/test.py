from pprint import pprint
from pathlib import Path

from MDANSE.Framework.Converters.ASE import ASE
from MDANSE.Framework.Converters.CASTEP import CASTEP
from MDANSE.Framework.Converters.Forcite import Forcite
from MDANSE.Framework.Jobs.Density import Density
from MDANSE.Framework.Parsers.FortranUnformat import binary_file_reader
from MDANSE.Framework.Parsers.xyz import XYZFile
from MDANSE.Framework.Parsers.trj import TrjFile

X = CASTEP()

X.fold = False
X.trajectory_file = "~/MDANSE/MDANSE/Tests/UnitTests/Data/PBAnew.md"
X.output_files.path = "./doom_file.mdt"
X.time_unit = "ps"
X.time_step = 10.0


# X.run()

# Y = ASE()

# Y.fold = False
# Y.trajectory_file = "~/MDANSE/MDANSE/Tests/UnitTests/Data/Cu_5steps_ASEformat.traj"
# Y.output_files.path = "./final_doom_file.mdt"
# Y.time_unit = "ps"
# Y.time_step = 10.0

# Y.run()

# Z = Density()
# Z.trajectory = X.output_files.path
# Z.frames = "all"
# Z.output_files.path = "./apm"

# Z.run()

# Q = XYZFile("~/MDANSE/MDANSE/Tests/UnitTests/Data/SrTiO3_MD-pos-1.xyz")
# print(Q)
# print(Q.time_step)
# for i in Q.frames:
#     print(i)

X = Forcite()

X.fold = False
X.xtd_file = "~/MDANSE/MDANSE/Tests/UnitTests/Data/H2O.xtd"
X.trj_file = "~/MDANSE/MDANSE/Tests/UnitTests/Data/H2O.trj"
X.output_files.path = "./doom_file.mdt"
X.time_unit = "ps"
X.time_step = 10.0

X.run()
