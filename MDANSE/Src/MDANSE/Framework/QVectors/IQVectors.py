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
import abc
from typing import TYPE_CHECKING

import numpy as np

from MDANSE.Core.Error import Error
from MDANSE.Framework.Configurable import Configurable
from MDANSE.Core.SubclassFactory import SubclassFactory
from MDANSE.MLogging import LOG

if TYPE_CHECKING:
    from MDANSE.MolecularDynamics.UnitCell import UnitCell
    from MDANSE.Framework.OutputVariables.IOutputVariable import OutputData


class QVectorsError(Error):
    pass


class IQVectors(Configurable, metaclass=SubclassFactory):
    is_lattice = False

    def __init__(self, atom_configuration, status=None):
        Configurable.__init__(self)

        self._atom_configuration = atom_configuration

        self._status = status

    @abc.abstractmethod
    def _generate(self):
        pass

    def generate(self) -> bool:
        if self._configured:
            self._generate()

            if self._status is not None:
                self._status.finish()
            return True
        else:
            LOG.error(
                f"Cannot generate vectors: q vector generator is not configured correctly."
            )
            return False

    def setStatus(self, status):
        self._status = status

    @classmethod
    def qvectors_to_hkl(
        self, vector_array: np.array, unit_cell: "UnitCell"
    ) -> np.ndarray:
        """Using a unit cell definition, recalculates an array
        of q vectors to an equivalent array of HKL Miller indices.

        Parameters
        ----------
        vector_array : np.array
            a (3,N) array of scattering vectors
        unit_cell : UnitCell
            an instance of UnitCell class describing the simulation box

        Returns
        -------
        np.ndarray
            A (3,N) array of HKL values (Miller indices)
        """
        return np.dot(unit_cell.direct, vector_array) / (2 * np.pi)

    @classmethod
    def hkl_to_qvectors(self, hkls: np.array, unit_cell: "UnitCell") -> np.ndarray:
        """Converts an array of HKL values to scattering vectors
        based on the definition of a unit cell.

        Parameters
        ----------
        hkls : np.array
            A (3,N) array of HKL values (Miller indices)
        unit_cell : UnitCell
            An instance of UnitCell class describing the simulation box shape

        Returns
        -------
        np.ndarray
            a (3, N) array of Q vectors (scattering vectors)
        """
        return 2 * np.pi * np.dot(unit_cell.inverse, hkls)

    def write_vectors_to_file(self, output_data: "OutputData"):
        """Writes a summary of the generated vectors to the output
        file using an OutputData class instance.

        Parameters
        ----------
        output_data : OutputData
            An object managing the writeout to one or many output files
        """
        q_values = [float(x) for x in self._configuration["q_vectors"].keys()]
        output_data.add(
            "vector_generator_q",
            "LineOutputVariable",
            q_values,
            units="1/nm",
        )
        qvector_lengths = [
            self._configuration["q_vectors"][q]["q_vectors"].shape[1] for q in q_values
        ]
        qarray_maxlength = np.max(qvector_lengths)
        output_data.add(
            "vector_generator_qvector_array",
            "VolumeOutputVariable",
            (len(q_values), 3, qarray_maxlength),
            units="1/nm",
        )
        output_data["vector_generator_qvector_array"][:] = 0.0
        for nq, q in enumerate(q_values):
            output_data["vector_generator_qvector_array"][
                nq, :, : qvector_lengths[nq]
            ] = self._configuration["q_vectors"][q]["q_vectors"]
        try:
            hkl_arrays = [self._configuration["q_vectors"][q]["hkls"] for q in q_values]
        except KeyError:
            return
        output_data.add(
            "vector_generator_hkl_array",
            "VolumeOutputVariable",
            (len(q_values), 3, qarray_maxlength),
            units="au",
        )
        output_data["vector_generator_hkl_array"][:] = 0.0
        for nq, q in enumerate(q_values):
            output_data["vector_generator_hkl_array"][nq, :, : qvector_lengths[nq]] = (
                hkl_arrays[nq]
            )
