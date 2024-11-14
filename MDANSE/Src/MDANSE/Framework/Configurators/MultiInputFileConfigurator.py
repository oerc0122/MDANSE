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
import ast
import os

from MDANSE import PLATFORM
from MDANSE.Framework.Configurators.IConfigurator import IConfigurator


class MultiInputFileConfigurator(IConfigurator):

    _default = []

    def __init__(self, name, wildcard="All files (*)", **kwargs):
        IConfigurator.__init__(self, name, **kwargs)
        self._wildcard = wildcard

    def configure(self, values):
        """
        Parameters
        ----------
        values : str or list
            A list of file names or a string of the list which can be
            converted to a list using literal_eval.
        """
        self["values"] = self._default
        self._original_input = values

        if type(values) is str:
            if values:
                try:
                    values = ast.literal_eval(values)
                except (SyntaxError, ValueError) as e:
                    self.error_status = f"Unable to evaluate string: {e}"
                    return
                if type(values) is not list:
                    self.error_status = (
                        f"Input values should be able to be evaluated as a list"
                    )
                    return
            else:
                values = []

        if type(values) is list:
            if not all([type(value) is str for value in values]):
                self.error_status = f"Input values should be a list of str"
                return
        else:
            self.error_status = f"Input values should be able to be evaluated as a list"
            return

        values = [PLATFORM.get_path(value) for value in values]

        none_exist = []
        for value in values:
            if not os.path.isfile(value):
                none_exist.append(value)

        if none_exist:
            self.error_status = f"The files {', '.join(none_exist)} do not exist."
            return

        self["values"] = values
        self["filenames"] = values
        self.error_status = "OK"

    @property
    def wildcard(self):
        return self._wildcard

    def get_information(self) -> str:
        """
        Returns
        -------
        str
           Input file names.
        """
        try:
            filenames = self["value"]
        except KeyError:
            filenames = ["None"]
        return f"Input files: {', '.join(filenames)}.\n"
