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

import codecs
import io
import tarfile
import time
from importlib import metadata
from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np

from MDANSE import PLATFORM
from MDANSE.Framework.Formats.IFormat import IFormat

if TYPE_CHECKING:
    from MDANSE.Framework.Jobs.IJob import IJob


def length_stringio(input: io.BytesIO) -> int:
    result = input.getbuffer().nbytes
    return result


class TextFormat(IFormat):
    """
    This class handles the writing of output variables in Text format. Each output variable is written into separate Text files which are further
    added to a single archive file.
    """

    extension = ".dat"

    extensions = [".dat", ".txt"]

    @classmethod
    def write(
        cls,
        filename: Path | str,
        data,
        header: str = "",
        run_instance: IJob = None,
    ):
        """
        Write a set of output variables into a set of Text files.

        Each output variable will be output in a separate Text file. All the Text files will be compressed into a tar file.

        :param filename: the path to the output archive file that will contain the Text files written for each output variable.
        :type filename: str
        :param data: the data to be written out.
        :type data: dict of Framework.OutputVariables.IOutputVariable
        :param header: the header to add to the output file.
        :type header: str
        """

        filename = Path(filename)
        filename = filename.parent / (filename.stem + "_text.tar")

        PLATFORM.create_directory(filename.parent)
        tf = tarfile.open(filename, "w")

        if header:
            real_buffer = io.BytesIO()
            temp_str = codecs.getwriter("utf-8")(real_buffer)
            for line in header:
                temp_str.write(str(line))
            temp_str.write("\n\n")
            real_buffer.seek(0)
            info = tarfile.TarInfo(name="jobinfo.txt")
            info.size = length_stringio(real_buffer)
            info.mtime = time.time()
            tf.addfile(tarinfo=info, fileobj=real_buffer)

        if run_instance is not None:
            inputs = run_instance.output_configuration()
            real_buffer = io.BytesIO()
            temp_str = codecs.getwriter("utf-8")(real_buffer)
            temp_str.write(f"run type: {run_instance.__class__.__name__}\n")
            temp_str.write(f"MDANSE version: {metadata.version('MDANSE')}\n")
            for key, value in inputs.items():
                temp_str.write(f"parameters[{str(key)}] = {str(value)}\n")
            temp_str.write("\n\n")
            real_buffer.seek(0)
            info = tarfile.TarInfo(name="job_parameters.txt")
            info.size = length_stringio(real_buffer)
            info.mtime = time.time()
            tf.addfile(tarinfo=info, fileobj=real_buffer)

        for var in list(data.values()):
            real_buffer = io.BytesIO()
            temp_str = codecs.getwriter("utf-8")(real_buffer)
            temp_str.write(var.info())
            temp_str.write("\n\n")
            cls.write_data(temp_str, var, data)
            real_buffer.seek(0)

            info = tarfile.TarInfo(name=f"{var.varname}{cls.extensions[0]}")
            info.size = length_stringio(real_buffer)
            info.mtime = time.time()
            tf.addfile(tarinfo=info, fileobj=real_buffer)

        tf.close()

    @classmethod
    def write_data(cls, fileobject, data, all_data):
        """
        Write an Framework.OutputVariables.IOutputVariable into a file-like object

        :param fileobject: the file object where the output variable should be written.
        :type fileobject: python file-like object
        :param data: the output variable to write (subclass of NumPy array).
        :type data: Framework.OutputVariables.IOutputVariable
        :param all_data: the complete set of output variables
        :type all_data: dict of Framework.OutputVariables.IOutputVariable

        :attention: this is a recursive method.
        """

        if data.ndim > 2:
            fileobject.write("Can not write Text output for data of dimensionality > 2")

        elif data.ndim == 2:
            x_data, y_data = data.axis.split("|")

            if x_data == "index":
                x_values = np.arange(data.shape[0])
                fileobject.write(f"# 1st column: {x_data} (au)\n")
            else:
                x_values = all_data[x_data]
                fileobject.write(f"# 1st column: {x_values.varname} ({x_values.units})\n")

            if y_data == "index":
                y_values = np.arange(data.shape[1])
                fileobject.write(f"# 1st row: {y_data} (au)\n\n")
            else:
                y_values = all_data[y_data]
                fileobject.write(f"# 1st row: {y_values.varname} ({y_values.units})\n\n")

            if np.allclose(np.imag(data), 0.0):
                z_data = np.zeros(
                    (data.shape[0] + 1, data.shape[1] + 1), dtype=np.float64
                )
                data = np.real(data)
            else:
                z_data = np.zeros(
                    (data.shape[0] + 1, data.shape[1] + 1), dtype=np.complex128
                )
            z_data[1:, 0] = x_values
            z_data[0, 1:] = y_values
            z_data[1:, 1:] = data

            np.savetxt(fileobject, z_data)
            fileobject.write("\n")

        else:
            x_data = data.axis.split("|")[0]

            if x_data == "index":
                x_values = np.arange(data.size)
                fileobject.write(f"# 1st column: {x_data} (au)\n")
            else:
                x_values = all_data[x_data]
                fileobject.write(f"# 1st column: {x_values.varname} ({x_values.units})\n")

            fileobject.write(f"# 2nd column: {data.varname} ({data.units})\n\n")

            np.savetxt(fileobject, np.column_stack([x_values, data]))
            fileobject.write("\n")
