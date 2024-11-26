import unittest

from MDANSE.Chemistry.ChemicalSystem import Molecule, Protein, ChemicalSystem
from MDANSE.MolecularDynamics.Configuration import RealConfiguration
from MDANSE.MolecularDynamics.TrajectoryUtils import *


class TestTrajectoryUtils(unittest.TestCase):
    def test_atom_index_to_molecule_index(self):
        m1 = Molecule("WAT", "w1")
        m2 = Molecule("WAT", "w2")
        cs = ChemicalSystem()
        cs.add_chemical_entity(m1)
        cs.add_chemical_entity(m2)

        result = atom_index_to_molecule_index(cs)
        self.assertDictEqual({0: 0, 1: 0, 2: 0, 3: 1, 4: 1, 5: 1}, result)

    def test_find_atoms_in_molecule(self):
        m1 = Molecule("WAT", "water")
        m2 = Molecule("WAT", "water")
        m3 = Molecule("WAT", "water")
        p = Protein("protein")
        p.set_peptide_chains([m3])

        cs = ChemicalSystem()
        for ce in [m1, m2, p]:
            cs.add_chemical_entity(ce)

        water_hydrogens = find_atoms_in_molecule(cs, "water", ["HW2", "HW1"])
        protein_oxygens = find_atoms_in_molecule(cs, "protein", ["OW"])
        empty = find_atoms_in_molecule(cs, "INVALID", ["OW"])
        empty2 = find_atoms_in_molecule(cs, "water", ["INVALID"])
        water_hydrogen_indices = find_atoms_in_molecule(
            cs, "water", ["HW2", "HW1"], True
        )

        self.assertEqual(
            [[m1["HW2"], m1["HW1"]], [m2["HW2"], m2["HW1"]]], water_hydrogens
        )
        self.assertEqual([[m3["OW"]]], protein_oxygens)
        self.assertEqual([], empty)
        self.assertEqual([[], []], empty2)
        self.assertEqual([[1, 2], [4, 5]], water_hydrogen_indices)

    def test_resolve_undefined_molecules_name(self):
        m1 = Molecule("WAT", "")
        m2 = Molecule("WAT", " water ")
        m4 = Molecule("WAT", "   ")

        m3 = Molecule("WAT", "")
        p = Protein("protein")
        p.set_peptide_chains([m3])

        cs = ChemicalSystem()
        for ce in [m1, m2, m4, p]:
            cs.add_chemical_entity(ce)

        resolve_undefined_molecules_name(cs)

        self.assertEqual("H2O1", m1.name)
        self.assertEqual(" water ", m2.name)
        self.assertEqual("H2O1", m4.name)
        self.assertEqual("", m3.name)
        self.assertEqual("protein", p.name)
