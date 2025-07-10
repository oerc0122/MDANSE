import h5py

from MDANSE.MolecularDynamics.Trajectory import Trajectory

from .BaseTypesDescriptor import PathConfigDesc


class MDANSETrajectoryFile(PathConfigDesc):
    def __init__(
        self,
        mode: None = None,
        *,
        extensions: None = None,
        directory: None = None,
        **params,
    ):
        super().__init__(mode="r", extensions=("mdt", "h5"), **params)

    def validate(self, value, *_) -> Trajectory:
        value = super().validate(value)

        return Trajectory(value)
