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

import multiprocessing
from warnings import warn

from MDANSE.Framework.Parameters import Integer, SingleChoice
from MDANSE.Framework.Parameters.Parameters import ConfigWarning, CustomConfig


class RunningMode(CustomConfig):
    mode = SingleChoice(
        choices=("single-core", "multicore", "remote"), default="single-core"
    )
    n_procs = Integer(minimum=1, maximum=multiprocessing.cpu_count(), default=1)

    def validate(self, desc, value):
        if self.mode == "single-core" and self.n_procs > 1:
            warn(
                "Requested more than one process for single-core job.",
                category=ConfigWarning,
            )
