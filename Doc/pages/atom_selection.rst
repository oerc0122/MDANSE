
.. _atom-selection:

Atom Selection
==============

An analysis does not have to include all the atoms present in the system.
For large, complicated or inhomogeneous systems it may be beneficial to select
subsets of atoms, and perform analysis on those subsets.
This can be done by applying atom selection.

Reusable Selection
~~~~~~~~~~~~~~~~~~

The current design choice in MDANSE is to create atom selection in a way
that allows the same definition to be used with different trajectories.
This is achieved by specifying the selection as a sequence of operations,
each one of them producing a set of atom indices.

A simple example of a reusable atom selection is this sequence:

1. Select all water molecules,
2. Invert selection.

This can be applied in a meaningful way to a wider group of trajectories,
such as macromolecules solvated by water. For each trajectory, this sequence
will remove the solvent from the selection, allowing us to run the analysis
on the macromolecule itself.

Selection Output
----------------

The result of applying a selection to a trajectory is a set of atom indices.

Selection Input
---------------

The input of atom selection is a sequence of selection operations.
Each operation returns a set of atom indices, and user input specifies
how each set should be included in the existing selection.

Set operations
~~~~~~~~~~~~~~

Three standard set operations are implemented in the atom selection:

Union
-----

For sets A and B, their union contains all elements of A and all elements of B.

Intersection
------------

For sets A and B, their intersection contains only those elements that belong
to both A and B.

Difference
----------

For sets A and B, their difference contains only those elements that belong
to A and at the same time do not belong to B.