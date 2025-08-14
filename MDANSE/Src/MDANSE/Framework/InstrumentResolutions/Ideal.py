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

import collections

import numpy as np

from MDANSE.Framework.InstrumentResolutions.IInstrumentResolution import (
    IInstrumentResolution,
)


class Ideal(IInstrumentResolution):
    """Defines an ideal instrument resolution with a Dirac response"""

    settings = collections.OrderedDict()

    def set_kernel(self, omegas, dt):
        n_omegas = len(omegas)
        self._omega_window = np.zeros(n_omegas, dtype=np.float64)
        self._omega_window[int(n_omegas / 2)] = 1.0

        self._time_window = np.ones(n_omegas, dtype=np.float64)
