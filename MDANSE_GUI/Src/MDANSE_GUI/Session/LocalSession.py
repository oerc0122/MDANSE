#    This file is part of MDANSE_GUI.
#
#    MDANSE_GUI is free software: you can redistribute it and/or modify
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

import json
from pathlib import Path

from qtpy.QtCore import QObject, Signal, Slot

from MDANSE import PLATFORM
from MDANSE.MLogging import LOG
from MDANSE_GUI.Tabs.Settings.LocalSettings import LocalSettings

json_encoder = json.encoder.JSONEncoder()
json_decoder = json.decoder.JSONDecoder()


class LocalSession(QObject):
    """Stores different parameters, and the state
    of the software.

    The intention is to have a single session that
    can be accessed by different tabs.

    At the moment, LocalSession is meant to be used
    for testing individual GUI components, and
    is largely a mock object.
    """

    new_units = Signal(dict)
    new_cmap = Signal(str)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._parameters = {}
        self._paths = {}
        self._units = {}
        self._colours = {}
        self._state = None
        self._filename = None
        self.populate_defaults()

    def populate_defaults(self):
        self._paths["root_directory"] = str(Path.home())
        self._units["energy"] = "meV"
        self._units["time"] = "fs"
        self._units["distance"] = "ang"
        self._units["reciprocal"] = "1/ang"
        self._colours["colormap"] = "viridis"
        self._colours["style"] = "ggplot"

    @Slot(dict)
    def update_units(self, inp: dict):
        self._units.update(inp)
        self.new_units.emit(self._units)

    @Slot(str)
    def update_cmap(self, inp: str):
        self._colours["colormap"] = input
        self.new_cmap.emit(self._colours["colormap"])

    def get_parameter(self, key: str) -> str:
        return self._parameters.get(key, None)

    def get_path(self, key: str) -> str:
        value = self._paths.get(key, str(Path.cwd()))
        return value

    def set_path(self, key: str, value: str):
        self._paths[key] = value

    def get_unit(self, key: str) -> str:
        value = self._units.get(key, "1")
        return value

    def obtain_settings(self, gui_element):
        return LocalSettings()

    def sections(self) -> list[dict[str, str]]:
        return [self._units, self._colours]

    @Slot()
    def save(self, fname: Path | str | None = None):
        all_items = {
            "paths": self._paths,
            "units": self._units,
            "colours": self._colours,
        }

        if fname is None:
            fname = (
                self._filename
                if self._filename
                else PLATFORM.application_directory() / "gui_session.json"
            )

        try:
            with open(fname, "w") as out_file:
                json.dump(all_items, out_file, cls=json_encoder)
        except Exception:
            return

        self._filename = fname

    def load(self, fname: Path | str | None = None):
        if fname is None:
            fname = PLATFORM.application_directory() / "gui_session.json"

        try:
            with open(fname, encoding="utf-8") as source:
                all_items = json.load(source, cls=json_decoder)
        except Exception:
            LOG.warning(f"Failed to read session settings from {fname}")
            return

        self._paths = all_items["paths"]
        self._units = all_items["units"]
        if "colours" in all_items:
            self._colours = all_items["colours"]
        self._filename = fname
