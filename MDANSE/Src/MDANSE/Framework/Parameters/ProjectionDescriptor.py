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

from enum import Enum, auto
from functools import singledispatchmethod

import numpy as np
import numpy.typing as npt

from MDANSE.Framework.Projectors.IProjector import IProjector
from MDANSE.MLogging import LOG

from .AbsConfigDesc import ConfigError, CustomConfig
from .BaseTypesDescriptor import Vector
from .ChoiceConfigDesc import SingleChoice


class ProjType(Enum):
    NULL = auto()
    AXIAL = auto()
    PLANAR = auto()

    NullProjector = NULL
    AxialProjector = AXIAL
    PlanarProjector = PLANAR


class Projection(CustomConfig):
    proj_type = SingleChoice(choices=ProjType)
    axis = Vector(
        optional=True,
        non_zero=True,
    )

    def __init__(
        self,
        proj_type: ProjType = ProjType.NULL,
        axis: npt.NDArray[float] | None = None,
        **kwargs,
    ):
        super().__init__(**kwargs)
        self.proj_type = proj_type
        self.axis = axis

    @property
    def projector(self) -> IProjector:
        if not self.check_status():
            LOG.warning(
                f"Projector ({self.proj_type.name}) requested, but no axis defined, returning Null"
            )
            return IProjector.create("null")

        proj = IProjector.create(self.proj_type.name.lower())

        if self.proj_type is not ProjType.NULL:
            proj.set_axis(self.axis)

    @projector.setter
    def projector(
        self,
        value: IProjector | tuple[ProjType | str, npt.NDArray[float] | None] | None,
    ) -> None:
        if isinstance(value, IProjector):
            self.proj_type = ProjType(type(value).__name__)

            if self.proj_type is not ProjType.NULL:
                self.axis = value._axis

        elif isinstance(value, tuple):
            name, axis = value
            self.proj_type = ProjType(name)
            self.axis = axis

        elif value is None:
            self.proj_type = ProjType.NULL
            self.axis = None

        else:
            raise ConfigError(f"Unable to parse {type(value).__name__} as IProjector")

    def check_status(self) -> bool:
        return (self.proj_type is ProjType.NULL) or (
            self.proj_type in (ProjType.AXIAL, ProjType.PLANAR)
            and self.first_axis is not None
        )

    def validate(self, _desc, _value) -> IProjector | None:
        return self.projector
