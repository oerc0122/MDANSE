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

from .AbsConfigDesc import to_class
from .AtomMappingDesc import (
    AtomMapping,
    AtomSelection,
    AtomTransmutation,
    GroupingLevel,
    PartialCharge,
    PartialChargeMapper,
)
from .BaseTypesDescriptor import (
    Array,
    Boolean,
    Float,
    Integer,
    NumericRange,
    PathParam,
    Range,
    String,
    Vector,
)
from .ChoiceConfigDesc import (
    DynamicMultiChoice,
    DynamicSingleChoice,
    MultipleChoice,
    SingleChoice,
)
from .Filter import TrajectoryFilter
from .FramesDescriptors import CorrelationWindow, Frames, FrameSelect, InterpOrder
from .InputFileDescriptors import MDANSEResult, MDANSETrajectory
from .OutputFileDescriptors import ASEOutputFormat, OutputFile, OutputTrajectory
from .ProjectionDescriptor import Projection
from .QVectorsDescriptor import QVectors, QVectorsSelect
from .RangeDescriptors import RangeCellCutoff
from .Resolution import InstrumentResolution
from .RunningModeDescriptor import RunningMode
from .Weights import Weights

__all__ = [
    "Array",
    "AtomMapping",
    "AtomSelection",
    "AtomTransmutation",
    "Boolean",
    "CorrelationWindow",
    "DynamicMultiChoice",
    "DynamicSingleChoice",
    "Filter",
    "Float",
    "FrameSelect",
    "Frames",
    "GroupingLevel",
    "InstrumentResolution",
    "Integer",
    "InterpOrder",
    "MDANSEResult",
    "MDANSETrajectory",
    "MultipleChoice",
    "NumericRange",
    "OutputFile",
    "OutputTrajectory",
    "PartialCharge",
    "PartialChargeMapper",
    "PathParam",
    "Projection",
    "QVectors",
    "QVectorsSelect",
    "Range",
    "RangeCellCutoff",
    "RunningMode",
    "SingleChoice",
    "String",
    "Vector",
    "Weights",
    "to_class",
]
