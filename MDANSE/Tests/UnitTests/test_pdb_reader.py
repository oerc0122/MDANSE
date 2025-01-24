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
import os
import unittest
from MDANSE.IO.MinimalPDBReader import MinimalPDBReader as PDBReader


pbd_2vb1 = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Data", "2vb1.pdb")


class TestPDBReader(unittest.TestCase):
    """
    Unittest for the geometry-related functions
    """

    def test_reader(self):
        with self.assertRaises((IOError, AttributeError)):
            reader = PDBReader("xxxxx.pdb")

        reader = PDBReader(pbd_2vb1)

        chemical_system = reader._chemical_system

        atomList = chemical_system.atom_list

        self.assertEqual(atomList[4], "C")
        print(chemical_system._labels.keys())
        self.assertTrue(10 in chemical_system._labels["LYS"])
        self.assertEqual(chemical_system.name_list[7], "HB2")
        # self.assertEqual(atomList[10].full_name, "...LYS1.HG2")
        # self.assertEqual(atomList[28].parent.name, "VAL2")
