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

import collections
import struct

import numpy as np

from MDANSE.Framework.Converters.Converter import Converter
from MDANSE.Framework.Units import measure
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter

FORCE_FACTOR = measure(1.0, "kcal_per_mole/ang", equivalent=True).toval("Da nm/ps2 mol")


class TrjFile(dict):
    def __init__(self, trjfilename):
        self["instance"] = open(trjfilename, "rb")

        self.parse_header()

    def parse_header(self):
        trjfile = self["instance"]

        rec = "!4x"
        rec_size = struct.calcsize(rec)
        trjfile.read(rec_size)

        # Record 1
        rec = "!4s20i8x"
        rec_size = struct.calcsize(rec)
        data = struct.unpack(rec, trjfile.read(rec_size))
        version = data[1]
        if version < 2010:
            self._fp = "f"
        else:
            self._fp = "d"

        # Diff with doc --> NTRJTI and TRJTIC not in doc
        rec = "!i"
        rec_size = struct.calcsize(rec)
        (ntrjti,) = struct.unpack(rec, trjfile.read(rec_size))
        rec = f"!{80 * ntrjti}s8x"
        rec_size = struct.calcsize(rec)
        self["title"] = struct.unpack(rec, trjfile.read(rec_size))
        self["title"] = "\n".join([t.decode("utf-8") for t in self["title"]])

        # Record 2
        rec = "!i"
        rec_size = struct.calcsize(rec)
        neexti = struct.unpack(rec, trjfile.read(rec_size))[0]
        rec = f"!{80 * neexti}s8x"
        rec_size = struct.calcsize(rec)
        trjfile.read(rec_size)

        # Record 3
        rec = "!8i8x"
        rec_size = struct.calcsize(rec)
        pertype, _, lcanon, defcel, _, _, lnpecan, ltmpdamp = struct.unpack(
            rec, trjfile.read(rec_size)
        )
        self["pertype"] = pertype
        self["defcel"] = defcel

        # Record 4
        rec = "!i"
        rec_size = struct.calcsize(rec)
        nflusd = struct.unpack(rec, trjfile.read(rec_size))[0]

        rec = f"!{nflusd}i{nflusd}i{8 * nflusd}s8x"
        rec_size = struct.calcsize(rec)
        trjfile.read(rec_size)

        rec = "!i"
        rec_size = struct.calcsize(rec)
        self["totmov"] = struct.unpack(rec, trjfile.read(rec_size))[0]

        rec = f"!{self['totmov']}i8x"
        rec_size = struct.calcsize(rec)
        self["mvofst"] = (
            np.array(struct.unpack(rec, trjfile.read(rec_size)), dtype=np.int32) - 1
        )

        # Record 4a
        rec = "!i"
        rec_size = struct.calcsize(rec)
        (leexti,) = struct.unpack(rec, trjfile.read(rec_size))
        rec = f"!{leexti}s8x"
        rec_size = struct.calcsize(rec)
        trjfile.read(rec_size)

        # Record 4b
        rec = "!i"
        rec_size = struct.calcsize(rec)
        (lparti,) = struct.unpack(rec, trjfile.read(rec_size))
        rec = f"!{lparti}s8x"
        rec_size = struct.calcsize(rec)
        trjfile.read(rec_size)

        self._header_size = trjfile.tell()

        # Frame record 1
        if version == 2000:
            rec1 = f"!{self._fp}i33{self.fp}5i8x"
        elif version == 2010:
            rec1 = f"!{self._fp}i57{self._fp}6i8x"
        else:
            rec1 = f"!{self._fp}i58{self._fp}6i8x"

        rec_size = struct.calcsize(rec1)
        data = struct.unpack(rec1, trjfile.read(rec_size))

        if version < 2010:
            self["velocities_written"] = data[-3]
            self["gradients_written"] = 0
        else:
            self["velocities_written"] = data[-4]
            self["gradients_written"] = data[-3]

        # Frame record 2
        rec = f"!12{self._fp}8x"
        rec_size = struct.calcsize(rec)
        trjfile.read(rec_size)

        # Frame record 3
        if lcanon:
            rec = f"!4{self._fp}8x"
            rec_size = struct.calcsize(rec)
            trjfile.read(rec_size)

        if pertype > 0:
            # Frame record 4
            self._def_cell_rec_pos = trjfile.tell() - self._header_size
            self._def_cell_rec = f"!22{self._fp}8x"
            self._def_cell_rec_size = struct.calcsize(self._def_cell_rec)
            trjfile.read(self._def_cell_rec_size)

        if pertype > 0:
            # Frame record 5
            rec = f"!i14{self._fp}8x"
            rec_size = struct.calcsize(rec)
            trjfile.read(rec_size)

        if lnpecan:
            # Frame record 6
            rec = f"!3{self._fp}8x"
            rec_size = struct.calcsize(rec)
            trjfile.read(rec_size)

        if ltmpdamp:
            # Frame record 7
            rec = f"!{self._fp}8x"
            rec_size = struct.calcsize(rec)
            trjfile.read(rec_size)

        self._config_rec_pos = trjfile.tell() - self._header_size

        rec_count = 3

        if self["velocities_written"]:
            rec_count += 3
        if self["gradients_written"]:
            rec_count += 3

        self._config_rec = "!" + (f"{self['totmov']}{self._fp}8x" * rec_count)

        self._config_rec_size = struct.calcsize(self._config_rec)
        trjfile.read(self._config_rec_size)

        self._frame_size = trjfile.tell() - self._header_size

        trjfile.seek(0, 2)

        self["n_frames"] = (trjfile.tell() - self._header_size) // self._frame_size

    def read_step(self, index):
        """ """

        trjfile = self["instance"]

        pos = self._header_size + index * self._frame_size

        trjfile.seek(pos, 0)

        rec = f"!{self._fp}"
        rec_size = struct.calcsize(rec)
        (time_step,) = struct.unpack(rec, trjfile.read(rec_size))

        if self["defcel"]:
            trjfile.seek(pos + self._def_cell_rec_pos, 0)
            cell = np.zeros((3, 3), dtype=np.float64)
            # ax,by,cz,bz,az,ay
            cell_data = np.array(
                struct.unpack(
                    self._def_cell_rec, trjfile.read(self._def_cell_rec_size)
                ),
                dtype=np.float64,
            )[2:8] * measure(1.0, "ang").toval("nm")
            cell[0, 0] = cell_data[0]
            cell[1, 1] = cell_data[1]
            cell[2, 2] = cell_data[2]
            cell[1, 2] = cell_data[3]
            cell[0, 2] = cell_data[4]
            cell[0, 1] = cell_data[5]

        else:
            cell = None

        trjfile.seek(pos + self._config_rec_pos, 0)

        config = struct.unpack(self._config_rec, trjfile.read(self._config_rec_size))

        rows = 1 + self["velocities_written"] + self["gradients_written"]

        config = np.transpose(np.reshape(config, (rows, 3, self["totmov"])))
        xyz = config[:, :, 0] * measure(1.0, "ang").toval("nm")

        if self["velocities_written"]:
            vel = config[:, :, 1] * measure(1.0, "ang/fs").toval("nm/ps")
        else:
            vel = None

        if self["gradients_written"]:
            gradients = config[:, :, 2] * FORCE_FACTOR
        else:
            gradients = None

        return time_step, cell, xyz, vel, gradients

    def close(self):
        self["instance"].close()


class Forcite(Converter):
    """Converts a Forcite trajectory to an MDT trajectory."""

    label = "Forcite"

    settings = collections.OrderedDict()
    settings["xtd_file"] = (
        "XTDFileConfigurator",
        {
            "wildcard": "XTD files (*.xtd);;All files (*)",
            "default": "INPUT_FILENAME.xtd",
            "label": "Input XTD file",
        },
    )
    settings["trj_file"] = (
        "InputFileConfigurator",
        {
            "wildcard": "TRJ files (*.trj);;All files (*)",
            "default": "INPUT_FILENAME.trj",
            "label": "Input TRJ file",
        },
    )
    settings["atom_aliases"] = (
        "AtomMappingConfigurator",
        {
            "default": "{}",
            "label": "Atom mapping",
            "dependencies": {"input_file": "xtd_file"},
        },
    )
    settings["fold"] = (
        "BooleanConfigurator",
        {"default": False, "label": "Fold coordinates into box"},
    )
    settings["output_files"] = (
        "OutputTrajectoryConfigurator",
        {
            "formats": ["MDTFormat"],
            "root": "xtd_file",
            "label": "MDANSE trajectory (filename, format)",
        },
    )

    def initialize(self):
        """
        Initialize the job.
        """
        super().initialize()

        self._atomic_aliases = self.configuration["atom_aliases"]["value"]

        self._xtdfile = self.configuration["xtd_file"]

        self._xtdfile.build_chemical_system(self._atomic_aliases)

        self._chemical_system = self._xtdfile.chemical_system

        self._trjfile = TrjFile(self.configuration["trj_file"]["filename"])

        # The number of steps of the analysis.
        self.n_steps = self._trjfile["n_frames"]

        if self._trjfile["velocities_written"]:
            self._velocities = np.zeros(
                (self._chemical_system.number_of_atoms, 3), dtype=np.float64
            )
        else:
            self._velocities = None

        if self._trjfile["gradients_written"]:
            self._gradients = np.zeros(
                (self._chemical_system.number_of_atoms, 3), dtype=np.float64
            )
        else:
            self._gradients = None

        # A trajectory is opened for writing.
        self._trajectory = TrajectoryWriter(
            self.configuration["output_files"]["file"],
            self._chemical_system,
            self.n_steps,
            positions_dtype=self.configuration["output_files"]["dtype"],
            chunking_limit=self.configuration["output_files"]["chunk_size"],
            compression=self.configuration["output_files"]["compression"],
            initial_charges=self.configuration["xtd_file"].get_atom_charges(),
        )

    def run_step(self, index):
        """Runs a single step of the job.

        @param index: the index of the step.
        @type index: int.

        @note: the argument index is the index of the loop note the index of the frame.
        """

        # The x, y and z values of the current frame.
        time, cell, xyz, velocities, gradients = self._trjfile.read_step(index)

        # If the universe is periodic set its shape with the current dimensions of the unit cell.
        conf = self._xtdfile._configuration
        movable_atoms = self._trjfile["mvofst"]
        conf["coordinates"][movable_atoms, :] = xyz
        if conf.is_periodic:
            if self._trjfile["defcel"]:
                conf.unit_cell = cell
            if self._configuration["fold"]["value"]:
                conf.fold_coordinates()

        if self._velocities is not None:
            self._velocities[movable_atoms, :] = velocities
            conf["velocities"] = self._velocities

        if self._gradients is not None:
            self._gradients[movable_atoms, :] = gradients
            conf["gradients"] = self._gradients

        self._trajectory.dump_configuration(
            conf,
            time,
            units={
                "time": "ps",
                "unit_cell": "nm",
                "coordinates": "nm",
                "velocities": "nm/ps",
                "gradients": "Da nm/ps2",
            },
        )

        return index, None

    def combine(self, index, x):
        """
        @param index: the index of the step.
        @type index: int.

        @param x:
        @type x: any.
        """

        pass

    def finalize(self):
        """
        Finalize the job.
        """

        self._trjfile.close()

        # Close the output trajectory.
        self._trajectory.write_standard_atom_database()
        self._trajectory.close()

        super().finalize()
