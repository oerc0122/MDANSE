Converged Results
=================
In MD a large number of parameters must be chosen to ensure
high quality results are obtained. Some of these include the size of the MD
box, the time step and the simulation length. The choice of these parameters
also depend on the analysis calculation, for example, the dynamics coherent
structure factor is much more difficult to obtain converged results
when compared to the dynamics incoherent structure factor. In this section
we will show how the calculations results change with different starting
parameters for a simple liquid argon system.


Simulation Box Size
~~~~~~~~~~~~~~~~~~~
Here we run a liquid argon trajectory with four different simulation box
sizes: 1.146, 2.292, 4.584 and 9.168 nm. The same atom density,
temperature, time step and simulation length is used for all cases. We
calculate the pair distribution function for all trajectories.

.. _figure-pdf:

.. figure:: ./Pictures/argon_pdf.png
   :align: center
   :width: 11.748cm
   :height: 9.393cm

   PDF calculated for a 60 ps MD simulation of liquid argon with a
   number of different MD box sizes. Blue, orange, green and red plots
   correspond to MD box sizes of 1.146, 2.292, 4.584 and 9.168 nm
   respectively.

Here we show the results for the PDF plotted to half the box size.
Our smallest box size 1.146 nm (blue in :numref:`figure-pdf`) is small
enough for the periodic image of the argon atoms to have an affect on
itself. Compared to our largest system size, the first peak is shifted
to a slightly longer distance, is too high and too broad.


Simulation Length
~~~~~~~~~~~~~~~~~
For analysis calculations which calculate a correlation function,
a balance between the length of your correlation function and the number of
configurations that you average over for each time step must be made.
Here we run the dynamic coherent structure factor calculation using
the liquid argon trajectory with the 4.584 nm box size and ideal
instrument resolution with two different correlation frames settings.

.. _figure-fqt:

.. figure:: ./Pictures/fqt_conv.png
   :align: center
   :width: 11.748cm
   :height: 9.393cm

   The coherent intermediate scattering function calculated for 120 ps
   from a 240 ps and 12 ns MD simulation of liquid argon plotted in blue
   and orange respectively.

.. _figure-sqw:

.. figure:: ./Pictures/sqw_conv.png
   :align: center
   :width: 11.748cm
   :height: 9.393cm

   The dynamics coherent structure factor calculated from a Fourier
   transform of the above coherent intermediate scattering functions
   using a 240 ps and 12 ns MD simulation of liquid argon plotted in blue
   and orange respectively.

In the first calculation (blue in :numref:`figure-fqt` and :numref:`figure-sqw`) we
use a correlation frames setting of (0, 2000, 1, 1000). So that the first
and last frames will be 0 and 2000 and the number of time
steps of the correlation function will be 1000. This will mean that for
this calculation each time step of the correlation function will be
averaged over 1001 = 2000 -- 1000 + 1 configurations. In :numref:`figure-fqt`
we can see :math:`F(q, t)` oscillate around zero, after Fourier transforming
we obtain a noisy :math:`S(q, f)` which is especially poor at zero energy.

In the second calculation (orange in :numref:`figure-fqt` and :numref:`figure-sqw`)
we use a correlation frames setting of (0, 100000, 1, 1000). So that each time
step of the correlation function will be averaged over 99001 = 100000 -- 1000 + 1
configurations but will be same length of the first calculation. Visually
we can see that :math:`F(q, t)` decays and stay close to zero, after Fourier
transforming we obtain a much smoother :math:`S(q, f)`. Clearly for dynamic
coherent structure factor calculations, to obtain high quality results longer
trajectories are needed so that a larger number of configurations are used per
time step of the correlation function.

In both calculation the correlation function time steps was set to 1000 which
corresponds to a time of 120 ps. From the :math:`F(q, t)`, we can see that
this is a sufficiently long time as the correlation function have decayed
to zero and does not change significantly beyond 30 ps or so. For other
calculations or system this might not have be the case and a more careful
choice for the correlation frames maybe required.
