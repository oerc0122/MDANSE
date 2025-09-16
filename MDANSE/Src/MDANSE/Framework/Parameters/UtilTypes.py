from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple, NewType

import numpy as np
import numpy.typing as npt
from typing_extensions import TypedDict, TypeVar

if TYPE_CHECKING:
    from collections.abc import Container, Iterable

    from MDANSE.Core.SubclassFactory import SubclassFactory
    from MDANSE.Framework.Parameters.Frames import Frames, FrameWindow
    from MDANSE.Framework.Parameters.Parameters import ConfigureDescriptor
    from MDANSE.Framework.Projectors.IProjector import IProjector
    from MDANSE.Framework.QVectors.IQVectors import IQVectors
    from MDANSE.Mathematics.Signal import Filter as SFilter
    from MDANSE.MolecularDynamics.Trajectory import Trajectory

P = TypeVar("P")
T = TypeVar("T")
CB = TypeVar("CB", default=T)
Num = TypeVar("Num", int, float)
K = TypeVar("K")
V = TypeVar("V")
CD = TypeVar("CD", bound="ConfigureDescriptor")
SC = TypeVar("SC", bound="SubclassFactory")

DescID = NewType("DescID", str)


class Depends(TypedDict, total=False, extra_items=object):
    trajectory: Trajectory
    frames: Frames
    window: FrameWindow
    selection: Container[int]
    transmutation: dict[int, str]
    generator: IQVectors
    projector: IProjector
    filter: SFilter


class PredictionResult(NamedTuple):
    """Strings used for showing the predicted output range in the GUI.

    :label: An explanation for the user of what physical property the numbers represent.
    :data: Data representing prediction if numeric otherwise attribute containing data.
    :unit: The internal unit of MDANSE of the numerical array. Used for converting to GUI units.
    """
    label: str
    data: Iterable[float] | str
    unit: str = "none"
