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

from MDANSE.Framework.ConfigDescriptors.AbsConfigDesc import Parameter
from MDANSE.Framework.Formats import OutputFormats
from MDANSE.MLogging import LogLevels
from MDANSE.MolecularDynamics.Trajectory import TrajectoryWriter

from .AbsConfigDesc import ConfigError
from .BaseTypesDescriptor import IntegerConfigDesc, PathConfigDesc
from .ChoiceConfigDesc import MultipleChoiceConfigDesc, SingleChoiceConfigDesc


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


class LogLevelConfigDesc(SingleChoiceConfigDesc):
    def __init__(self, **kwargs):
        super().__init__(choices={level.name for level in LogLevels}, **kwargs)

    def validate(self, value: LogLevels | int | str, *_deps) -> LogLevels:
        try:
            value = LogLevels(value)
        except Exception as err:
            raise ConfigError("Invalid config level.") from err


class OutputFormatConfigDesc(MultipleChoiceConfigDesc):
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


class OutputFileConfigDesc(Parameter):
    out_format = OutputFormatConfigDesc(("MDAFormat", "TextFormat"))
    log_level = LogLevelConfigDesc(default=LogLevels.NONE, label="Reporting log level.")
    path = PathConfigDesc(mode="w")

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

    def __str__(self) -> str:
        return f"{type(self).__name__}({self.path=}, {self.out_format=}, {self.log_level=})"

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


class OutputFilePlusMemConfigDesc(OutputFileConfigDesc):
    out_format = OutputFormatConfigDesc(("MDAFormat", "TextFormat", "FileInMemory"))


class OutputTrajectoryConfigDesc(OutputFileConfigDesc):
    out_format = OutputFormatConfigDesc(
        ("MDTFormat",),
    )
    dtype = SingleChoiceConfigDesc(
        choices=(np.float16, np.float32, np.float64),
        default=np.float64,
        aliases=(
            dict.fromkeys(("64", "float64", "float", float, 64), np.float64)
            | dict.fromkeys(("32", "float32", "single", 32), np.float32)
            | dict.fromkeys(("16", "float16", "half", 16), np.float16)
        ),
    )
    chunk_size = IntegerConfigDesc(minimum=32, maximum=65536, default=128)
    compression = SingleChoiceConfigDesc(
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

    def __str__(self) -> str:
        return f"{type(self).__name__}({self.path=}, {self.log_level=}, {self.dtype=}, {self.chunk_size=}, {self.compression=})"

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
