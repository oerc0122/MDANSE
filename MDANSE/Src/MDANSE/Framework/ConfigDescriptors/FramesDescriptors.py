from __future__ import annotations

from collections.abc import Sequence

from .AbsConfigDesc import ConfigError, GUIComponent
from .BaseTypesDescriptor import IntegerConfigDesc, RangeConfigDesc


class FramesConfigDesc(RangeConfigDesc):
    def validate(self, value: range, deps: dict[str, object]) -> range:
        if value in {"all", None}:
            value = (0, len(deps["trajectory"]), 1)

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

        self.validate_range(value.start, 0, len(deps["trajectory"]))
        self.validate_range(value.stop, 0, len(deps["trajectory"]))

        return value


class CorrelationFramesDesc(FramesConfigDesc, GUIComponent):
    frames = FramesConfigDesc(
        label="Frames to include in correlation.",
        tooltip="Select which frames are to be included.",
    )
    window = IntegerConfigDesc(
        minimum=1,
        label="Frames window",
        tooltip="Number of frames in correlation window.",
    )
