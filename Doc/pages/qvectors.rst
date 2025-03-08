
.. _qvector-generation:

q-Vector Generation
~~~~~~~~~~~~~~~~~~~

Reciprocal Lattice q-Vectors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`b_1`, :math:`b_2`, :math:`b_3` be the basis vectors
which span the MD cell. Any position vector in the MD cell can be
written as

.. math::
   :label: pfx86

   {{R = x^{'}}{b_{1} + y^{'}}{b_{2} + z^{'}}b_{3},}

so that it defines the position vector in the MD cell.
With :math:`x'`, :math:`y'`, :math:`z'` having
values between :math:`0` and :math:`1` if :math:`R` is in the unit cell.
The primes indicate that the coordinates are fractional coordinates. A jump due
to periodic boundary conditions can cause :math:`x'`, :math:`y'`,
:math:`z'` to jump by :math:`\pm1`. The set of dual basis
vectors :math:`b^1`, :math:`b^2`, :math:`b^3` is defined by
the relation

.. math::
   :label: pfx87

   {b_{i}{b^{j} = \delta_{i}^{j}}.}

and defines the dual basis vectors and their relation to the basis
vectors. If the :math:`q`-vectors are now chosen as

.. math::
   :label: pfx88

   {{q = 2}\pi\left( {k{b^{1} + l}{b^{2} + m}b^{3}} \right),}

so that this selection of :math:`q`-vectors produces phase changes for
handling jumps in particle trajectories. Here :math:`k`, :math:`l`, :math:`m`
are integer numbers, jumps in the particle trajectories
produce phase changes of multiples of :math:`2\pi` in the Fourier transformed
particle density, i.e. leave it unchanged.

In MDANSE, One can define a grid of
:math:`q`-shells or a grid of :math:`q`-vectors along a given direction or on a
given plane, giving in addition a *tolerance* for :math:`q`. MDANSE looks
then for :math:`q`-vectors of the form given in Eq. :math:numref:`pfx88` whose moduli
deviate within the prescribed tolerance from the equidistant :math:`q`-grid.
From these :math:`q`-vectors only a maximum number per grid-point (called
generically :math:`q`-shell also in the anisotropic case) is kept.