from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Literal, NamedTuple

from more_itertools import numeric_range

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .AbsConfigDesc import ConfigError, CustomConfig
from .BaseTypesDescriptor import IntegerConfigDesc, RangeConfigDesc
from .ChoiceConfigDesc import DynamicSingleChoiceConfigDesc


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


class FramesConfigDesc(RangeConfigDesc):
    default_tooltip = "Select which frames are to be included."
    default_label = "Frames to include in correlation."

    def required_deps(self) -> set[str]:
        return super().required_deps() | {"trajectory"}

    def get_ranges(self, _value: None, deps: dict[str, Any]) -> tuple[int, int]:
        return (0, len(deps["trajectory"]) + 1)

    def validate(
        self,
        value: range | tuple[int, ...] | None | Literal["all"],
        deps: dict[str, Any],
    ) -> Frames:
        trajectory: Trajectory = deps["trajectory"]
        nsteps = len(trajectory)

        if value in {"all", None}:
            value = (0, nsteps, 1)

        value = super().validate(value, deps)

        return Frames(trajectory, value)


class CorrelationWindow(IntegerConfigDesc):
    default_label="Correlation window"
    default_tooltip="Number of frames in correlation window."

    def required_deps(self) -> set[str]:
        return super().required_deps() | {"frames"}

    def get_ranges(self, deps: dict[str, Any]):
        return (2, len(deps["frames"]))


class InterpOrder(DynamicSingleChoiceConfigDesc[int]):
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

    def __init__(self, *args, aliases: None = None, **kwargs):
        if aliases is not None:
            raise ConfigError("Aliases may not be non-None")

        super().__init__(self, *args, aliases=self.base_aliases, **kwargs)

    @property
    def choices(self) -> list[int]:
        return sorted(self.last_choices)

    def required_deps(self) -> set[str]:
        return super().required_deps() | {"trajectory", "frames"}

    def get_choices(self, deps) -> set[int]:
        maxi = min(6, len(deps["frames"]))
        mini = 0 if "velocities" in deps["trajectory"].variables() else 1
        return set(range(mini, maxi))
