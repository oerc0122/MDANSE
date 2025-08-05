from __future__ import annotations

from collections.abc import Sequence
from functools import cached_property
from typing import NamedTuple

from more_itertools import numeric_range

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .AbsConfigDesc import ConfigError, Parameter
from .BaseTypesDescriptor import IntegerConfigDesc, RangeConfigDesc


class Frame(NamedTuple):
    time: float
    loc_index: int
    index: int


class Frames:
    def __init__(self, traj: Trajectory, samples: range):
        self.samples = samples
        time = traj.time()
        self.times = numeric_range(
            time[samples[0]],
            time[samples[-1]],
            time[1] - time[0],
        )

    def __getitem__(self, value: int) -> Frame:
        return Frame(
            time=self.times[value],
            loc_index=value,
            index=self.index_start + (value * self.index_step),
        )

    def __len__(self) -> int:
        return len(self.samples)

    @property
    def index_start(self) -> int:
        return self.samples[0]

    @property
    def index_stop(self) -> int:
        return self.samples[-1]

    @property
    def index_step(self) -> int:
        return self.samples.step

    @property
    def time_start(self) -> int:
        return self.times[0]

    @property
    def time_stop(self) -> int:
        return self.times[-1]

    @property
    def time_step(self) -> int:
        return self.times.step


class FramesConfigDesc(RangeConfigDesc):
    def required_deps(self) -> set[str]:
        return super().required_deps() | {"trajectory"}

    def validate(self, value: range, deps: dict[str, object]) -> Frames:
        trajectory = deps["trajectory"]
        nsteps = len(trajectory)
        if value in {"all", None}:
            value = (0, nsteps, 1)

        if isinstance(value, int):
            value = range(value)
        elif isinstance(value, Sequence):
            value = range(*value)
        elif isinstance(value, dict):
            value = range(**value)
        elif not isinstance(value, range):
            raise ConfigError(
                f"Do not know how to convert {type(value).__name__} to range"
            )

        self.validate_range(value.start, (0, nsteps + 1))
        self.validate_range(value.stop, (0, nsteps + 1))

        return Frames(trajectory, value)


class CorrelationFramesDesc(FramesConfigDesc, Parameter):
    frames = FramesConfigDesc(
        depends={"trajectory": "trajectory"},
        label="Frames to include in correlation.",
        tooltip="Select which frames are to be included.",
    )
    window = IntegerConfigDesc(
        minimum=1,
        label="Frames window",
        tooltip="Number of frames in correlation window.",
    )
