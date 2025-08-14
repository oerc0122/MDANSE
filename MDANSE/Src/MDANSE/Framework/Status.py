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
import time


class Status(metaclass=abc.ABCMeta):
    """
    This class defines an interface for status objects.
    This kind of object is used to store the status a loop-based task.
    """

    def __init__(self):
        self._update_step = 1

        self._current_step = 0
        self._n_steps = None
        self._finished = False
        self._stopped = False
        self._start_time = time.time()
        self._deltas = [self._start_time, self._start_time + 1.0]
        self._elapsed_time = "N/A"
        self._last_refresh = self._start_time

    @abc.abstractmethod
    def finish_status(self):
        pass

    @abc.abstractmethod
    def start_status(self):
        pass

    @abc.abstractmethod
    def stop_status(self):
        pass

    @abc.abstractmethod
    def update_status(self):
        pass

    @abc.abstractmethod
    def fixed_status(self, current_progress: int):
        pass

    @property
    def current_step(self):
        return self._current_step

    @property
    def elapsed_time(self):
        return str(self._deltas[1] - self._deltas[0])

    def finish(self):
        self._finished = True

        self.finish_status()

    def get_current_step(self):
        return self._current_step

    def get_elapsed_time(self):
        return self.elapsed_time

    def get_number_of_steps(self):
        return self._n_steps

    def is_finished(self):
        return self._finished

    def is_running(self):
        return not self._finished and not self._stopped

    def is_stopped(self):
        return self._stopped

    @property
    def n_steps(self):
        return self._n_steps

    def start(self, n_steps, rate=None):
        if self._n_steps is not None:
            return

        self._n_steps = n_steps

        if rate is not None:
            self._update_step = max(0, int(rate * n_steps))

        self.start_status()

    def stop(self):
        self._stopped = True
        self.stop_status()

    def update(self, force=False):
        if self._update_step == 0:
            return

        self._current_step += 1

        last_update = time.time()

        self._deltas[1] = last_update

        if force or ((last_update - self._last_refresh) > 5):
            self._last_refresh = last_update

        self.update_status()
