from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal, NamedTuple

from more_itertools import numeric_range

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .Parameters import ConfigError, CustomConfig
from .BaseTypes import Integer, NumericRange, Range
from .Choices import SingleChoice
from .UtilTypes import Depends, DescID


class Frame(NamedTuple):
    time: float
    loc_index: int
    ind: int


class Frames:
    def __init__(self, traj: Trajectory, samples: range | NumericRange | Frames):
        if isinstance(samples, Frames):
            self.samples = samples.samples
        else:
            self.samples = samples
        time = traj.time()
        self.times = NumericRange(
            time[self.samples[0]],
            time[self.samples[-1]],
            time[1] - time[0],
        )

    def __getitem__(self, value: int) -> Frame:
        return Frame(
            time=self.times[value],
            loc_index=value,
            ind=self.index_start + (value * self.index_step),
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
    def time_start(self) -> float:
        return self.times[0]

    @property
    def time_stop(self) -> float:
        return self.times[-1]

    @property
    def time_step(self) -> float:
        return self.times.step

    @property
    def duration(self) -> float:
        return self.time_stop - self.time_start


class FrameSelect(Range[int]):
    default_tooltip = "Select which frames are to be included."
    default_label = "Frames to include in correlation."

    def __init__(self, *args, dtype: None = None, **kwargs):
        super().__init__(*args, dtype=int, **kwargs)

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def get_ranges(self, deps: Depends) -> tuple[int, int]:
        return (0, len(deps["trajectory"]) + 1)

    def validate(
        self,
        value: range | tuple[int, ...] | None | Literal["all"],
        deps: dict[str, Any],
        /,
    ) -> Frames:
        trajectory: Trajectory = deps["trajectory"]
        nsteps = len(trajectory)

        if value in {"all", None}:
            value = (0, nsteps, 1)

        ranges = super().validate(value, deps)

        print(ranges)

        return Frames(trajectory, ranges)


class CorrelationWindow(Integer):
    default_label = "Correlation window"
    default_tooltip = "Number of frames in correlation window."

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("frames")}

    def get_ranges(self, deps: Depends):
        return (2, len(deps["frames"]))


class InterpOrder(SingleChoice[int | str, int]):
    default_label = "Velocity interpolation order."
    default_tooltip = (
        "Select velocity interpolation order, "
        "if the file contains velocities use 0 to use those."
    )

    base_aliases = {
        "no interpolation": 0,
        "1st order": 1,
        "2nd order": 2,
        "3rd order": 3,
        "4th order": 4,
        "5th order": 5,
    }

    def __init__(self, *args, aliases: None = None, choices: None = None, **kwargs):
        if aliases is not None:
            raise ConfigError("aliases may not be non-None")
        if choices is not None:
            raise ConfigError("choices may not be non-None")

        super().__init__(*args, choices=(), aliases=self.base_aliases, **kwargs)

    @property
    def choices(self) -> list[int]:
        assert isinstance(self.last_choices, set)
        return sorted(self.last_choices)

    @choices.setter
    def choices(self, value) -> None: ...

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory"), DescID("frames")}

    def get_choices(self, deps: Depends) -> set[int]:
        maxi = min(6, len(deps["frames"]))
        mini = 0 if "velocities" in deps["trajectory"].variables() else 1
        return set(range(mini, maxi))
