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

from .AtomMappingDesc import (
    AtomMapping,
    AtomSelection,
    AtomTransmutation,
    GroupingLevel,
    PartialCharge,
    PartialChargeMapper,
)
from .BaseTypesDescriptor import (
    ArrayConfigDesc,
    BooleanConfigDesc,
    FloatConfigDesc,
    IntegerConfigDesc,
    NumericRangeConfigDesc,
    PathConfigDesc,
    RangeConfigDesc,
    StringConfigDesc,
    VectorConfigDesc,
)
from .ChoiceConfigDesc import (
    DynamicMultiChoiceConfigDesc,
    DynamicSingleChoiceConfigDesc,
    MultipleChoiceConfigDesc,
    SingleChoiceConfigDesc,
)
from .FramesDescriptors import CorrelationWindow, FramesConfigDesc, InterpOrder
from .InputFileDescriptors import MDANSETrajectoryFile
from .OutputFileDescriptors import (
    ASEOutputFormat,
    OutputFileConfigDesc,
    OutputTrajectoryConfigDesc,
)
from .ProjectionDescriptor import ProjectionConfigDesc
from .RangeDescriptors import RangeCellCutoff
from .Resolution import InstrumentResolution
from .RunningModeDescriptor import RunningModeConfigDesc
from .Weights import Weights

__all__ = [
    "ArrayConfigDesc",
    "AtomMapping",
    "AtomSelection",
    "AtomTransmutation",
    "BooleanConfigDesc",
    "CorrelationFramesDesc",
    "CorrelationWindow",
    "DynamicMultiChoiceConfigDesc",
    "DynamicSingleChoiceConfigDesc",
    "FloatConfigDesc",
    "FramesConfigDesc",
    "GroupingLevel",
    "InstrumentResolution",
    "IntegerConfigDesc",
    "InterpOrder",
    "MDANSETrajectoryFile",
    "MultipleChoiceConfigDesc",
    "NumericRangeConfigDesc",
    "OutputFileConfigDesc",
    "OutputTrajectoryConfigDesc",
    "PartialCharge",
    "PartialChargeMapper",
    "PathConfigDesc",
    "ProjectionConfigDesc",
    "RangeCellCutoff",
    "RangeConfigDesc",
    "RunningModeConfigDesc",
    "SingleChoiceConfigDesc",
    "StringConfigDesc",
    "VectorConfigDesc",
    "Weights",
]
