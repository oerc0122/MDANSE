from __future__ import annotations

import math
from typing import Any, Literal, NamedTuple, SupportsInt

import numpy as np

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .BaseTypes import Integer, NumericRange, Range
from .Choices import SingleChoice
from .Parameters import ConfigError
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

        if any(time := traj.time()):
            step = time[1] - time[0]
            self.times = NumericRange(
                time[self.samples[0]],
                time[self.samples[-1]] + (step / 2),
                step,
            )
        else:
            self.times = self.samples

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

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.index_start}, {self.index_stop}, {self.index_step})"


FrameInput = Frames | range | tuple[int, ...] | None | Literal["all"]


class FrameWindow:
    def __init__(self, window: int, frames: Frames):
        self.window = window
        self._frames = frames

    @property
    def n_frames(self):
        return self.window

    @property
    def n_configs(self):
        return len(self._frames) - self.window + 1

    @property
    def times(self):
        return self._frames.times[: self.window]

    @property
    def duration(self):
        return [time - self.times[0] for time in self.times]

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.window})"


class FrameSelect(Range[int]):
    default_tooltip = "Select which frames are to be included."
    default_label = "Frames to include in correlation."

    def __init__(
        self, *args, default: FrameInput = "all", dtype: None = None, **kwargs
    ):
        super().__init__(*args, default=default, dtype=int, **kwargs)

    def required_deps(self) -> set[DescID]:
        return super().required_deps() | {DescID("trajectory")}

    def get_ranges(self, deps: Depends) -> tuple[int, int]:
        return (0, len(deps["trajectory"]) + 1)

    def __set__(self, owner: object, value: tuple[SupportsInt, ...]):

        # Support legacy
        corr = [x for x in self.dependents if isinstance(x, CorrelationWindow)]
        window = None
        if len(value) == 4 and corr:
            value, window = value[:3], value[3]

        super().__set__(owner, value)

        if corr and window:
            corr[0].__set__(owner, window)

    def validate(
        self,
        value: FrameInput,
        deps: dict[str, Any],
        /,
    ) -> Frames:
        if isinstance(value, Frames):
            value = value.samples

        trajectory: Trajectory = deps["trajectory"]
        nsteps = len(trajectory)

        if value is None or value == "all":
            value = (0, nsteps, 1)

        ranges = super().validate(value, deps)
        return Frames(trajectory, ranges)


class CorrelationWindow(Integer):
    default_label = "Correlation window"
    default_tooltip = "Number of frames in correlation window."

    def __init__(self, *args, default: int | None = None, **kwargs):
        super().__init__(*args, default=default, optional=True, **kwargs)

    def validate(self, value: SupportsInt, deps: Depends):
        if isinstance(value, FrameWindow):
            value = value.window
        elif value is None:
            value = math.ceil(len(deps["frames"]) / 2)

        value = super().validate(value, deps)

        return FrameWindow(value, deps["frames"])

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
        if getattr(self, "last_choices", None) is None:
            self.last_choices = self.get_choices()

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
