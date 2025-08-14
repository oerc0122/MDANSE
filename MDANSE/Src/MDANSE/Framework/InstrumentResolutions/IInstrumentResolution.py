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

import abc

import numpy as np

from MDANSE.Core.Error import Error
from MDANSE.Core.SubclassFactory import SubclassFactory
from MDANSE.Framework.Configurable import Configurable


class InstrumentResolutionError(Error):
    pass


class IInstrumentResolution(Configurable, metaclass=SubclassFactory):
    def __init__(self):
        Configurable.__init__(self)

        self._omega_window = None

        self._time_window = None

    @abc.abstractmethod
    def set_kernel(self, omegas, dt):
        pass

    def apply_fft(self, omega_window, dt):
        time_window = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(omega_window)) / dt)
        return time_window

    @property
    def omega_window(self):
        if self._omega_window is None:
            raise InstrumentResolutionError("Undefined omega window")

        return self._omega_window

    @property
    def time_window(self):
        if self._time_window is None:
            raise InstrumentResolutionError("Undefined time window")

        return self._time_window
