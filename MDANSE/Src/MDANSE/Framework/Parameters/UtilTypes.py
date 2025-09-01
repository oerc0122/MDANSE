from __future__ import annotations

from collections.abc import Container
from typing import TYPE_CHECKING, Any, NewType, TypedDict

import numpy as np
import numpy.typing as npt

if TYPE_CHECKING:
    from MDANSE.Framework.Parameters.Frames import Frames
    from MDANSE.Framework.Projectors.IProjector import IProjector
    from MDANSE.Framework.QVectors.IQVectors import IQVectors
    from MDANSE.Mathematics.Signal import Filter as SFilter
    from MDANSE.MolecularDynamics.Trajectory import Trajectory

DescID = NewType("DescID", str)


class Depends(TypedDict, total=False):
    trajectory: Trajectory
    frames: Frames
    selection: Container[int]
    transmutation: dict[int, str]
    generator: IQVectors
    projector: IProjector
    filter: SFilter
