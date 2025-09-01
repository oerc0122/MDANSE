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

from collections.abc import Iterable
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
from ase.io.formats import ioformats

from MDANSE.Framework.Formats import OutputFormats
from MDANSE.Framework.Parameters.Parameters import Parameter
from MDANSE.MLogging import LogLevels
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter

from .BaseTypes import Integer, PathParam
from .Choices import MultipleChoice, SingleChoice
from .Parameters import ConfigError, CustomConfig


class OldOutputSettings(NamedTuple):
    path: Path | str
    format: OutputFormats | str | int | Iterable[OutputFormats | str | int]
    loglevel: LogLevels | int | str


class OldTrajectorySettings(NamedTuple):
    path: Path | str
    dtype: int
    chunk_size: int
    compression: str
    format: OutputFormats | str | int | Iterable[OutputFormats | str | int]
    loglevel: LogLevels | int | str


class OutputFormat(MultipleChoice):
    def __init__(
        self,
        allowed_formats: Iterable[OutputFormats | str | int],
        **kwargs,
    ):
        allowed_formats = map(OutputFormats, allowed_formats)
        super().__init__(
            n_choices=None,
            choices={format.name for format in allowed_formats},
            **kwargs,
        )


class OutputFile(CustomConfig):
    out_format = OutputFormat(("MDAFormat", "TextFormat"))
    log_level = SingleChoice(
        choices=LogLevels,
        default=LogLevels.NONE,
        label="Reporting log level.",
    )
    path = PathParam(mode="w")

    def __init__(
        self,
        formats: Iterable[OutputFormats | str | int] = ("MDAFormat", "TextFormat"),
        level: LogLevels = LogLevels.NONE,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.log_level = level
        self.out_format = formats

    @property
    def root(self) -> Path:
        return self.path.with_suffix("")

    @property
    def write_logs(self) -> bool:
        return self.log_level is not None

    @property
    def configuration(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "log_level": self.log_level,
            "formats": self.out_format,
        }

    @classmethod
    def from_old_tuple(cls, args: OldOutputSettings):
        out = cls()
        out.path, out.out_format, out.log_level = args
        return out


class OutputFilePlusMem(OutputFile):
    out_format = OutputFormat(("MDAFormat", "TextFormat", "FileInMemory"))


class OutputTrajectory(OutputFile):
    out_format = OutputFormat(
        ("MDTFormat",),
    )
    dtype = SingleChoice(
        choices=(np.float16, np.float32, np.float64),
        default=np.float64,
        aliases=(
            dict.fromkeys(("64", "float64", "float", float, 64), np.float64)
            | dict.fromkeys(("32", "float32", "single", 32), np.float32)
            | dict.fromkeys(("16", "float16", "half", 16), np.float16)
        ),
    )
    chunk_size = Integer(minimum=32, maximum=65536, default=128)
    compression = SingleChoice(
        choices=("none", *TrajectoryWriter.allowed_compression), default="none"
    )

    def __init__(
        self,
        level: LogLevels = LogLevels.NONE,
        dtype: int | str | None = None,
        chunk_size: int | None = None,
        compression: str = "none",
        **kwargs,
    ):
        super().__init__(formats=("MDTFormat",), level=level, **kwargs)
        if dtype is not None:
            self.dtype = str(dtype)
        if chunk_size is not None:
            self.chunk_size = chunk_size
        self.compression = self.compression

    @classmethod
    def from_old_tuple(cls, args: OldTrajectorySettings):
        """Construct object from old-style tuple."""
        out = cls()
        out.path, out.dtype, out.chunk_size, out.compression, out.log_level = args
        return out

    @property
    def dtype_size(self) -> int:
        return int(str(self.dtype)[-4:-2])

    @property
    def configuration(self) -> dict[str, Any]:
        return super().configuration() | {
            "dtype": self.dtype_size,
            "chunk_size": self.chunk_size,
            "compression": self.compression,
        }


class ASEOutputFormat(OutputFile):
    out_format = SingleChoice(
        choices=(key for key, val in ioformats.items() if val.can_write)
    )

    def __init__(self, fmt):
        self.out_format = fmt
