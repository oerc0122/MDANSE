
.. _appendix-fca:

Correlation Functions
=====================

Most of the quantities which can be extracted from MD
simulations are time correlation functions. In MDANSE we use a correlation
window to ensure that the time averaging for each time step
is done in a consistent way. Consider two time series

.. math::
   :label: eqn-fca1

   A(k \Delta t) \quad \text{and} \quad B(k \Delta t) \qquad k = 0, \ldots, n_{\mathrm{t}}-1,

of length :math:`t_{\mathrm{tot}} = (n_{\mathrm{t}} -1) \Delta t` which are
to be correlated. In MDANSE,
correlation function are calculated by first choosing a specific
number of correlation time steps :math:`n_{\mathrm{c}}` which will define
the length of our correlation function :math:`t_{\mathrm{cor}} = (n_{\mathrm{c}} -1) \Delta t`. The correlation function of
:math:`A(k \Delta t)` and :math:`B(k \Delta t)` will be

.. math::
   :label: eqn-fca2

   C_{AB}(l \Delta t) = \frac{1}{n_{\mathrm{t}} - n_{\mathrm{c}} + 1} \sum\limits_{k=0}^{n_{\mathrm{t}} - n_{\mathrm{c}} + 1} A^{*}(k\Delta t)B([k + l]\Delta t) \qquad l = 0, \ldots, n_{c} - 1.

In case that :math:`A(k \Delta t)` and
:math:`B(k \Delta t)` are identical, the corresponding correlation function
:math:`C_{AA}(l \Delta t)` is called an *autocorrelation* function. Notice that
the prefactor is the same for all :math:`l \Delta t` time steps, this was
not the case in previous versions of MDANSE. This meant that for different
time steps a different number of configurations were used to obtain the
average correlation; leading to spuriously large correlations for some
time intervals. However in MDANSE 2 your correlation functions will be
truncated since :math:`t_{\mathrm{tot}} \geq t_{\mathrm{cor}}`.

In many cases not only is the computation of a correlation function
required, but also the computation of its Fourier spectrum. In
MDANSE the spectra can be smoothed by applying a instrument resolution
function

.. math::
   :label: eqn-fca11

   P_{AB}\left(m \Delta \omega \right) = \frac{\Delta t}{2 \pi}\sum_{l=-(n_{\mathrm{c}}-1)}^{n_{\mathrm{c}}-1}
   \exp\left[- 2 \pi i \frac{l \Delta t }{2n_{\mathrm{c}} - 1} m \Delta \omega \right] \frac{W(l \Delta t)}{W(0)} C_{AB}( \vert l \Delta t \vert )

here :math:`m = -(n_{\mathrm{c}}-1), \ldots, n_{\mathrm{c}}-1` and :math:`\Delta \omega` is the
frequency step. Notice that the we assume that the correlation is symmetric so
that :math:`C_{AB}(l \Delta t) = C_{AB}( |l \Delta t| )` which should
approximately be the case for all the correlation functions calculated
in MDANSE assuming good (equilibrated, of a sufficient length/size and
etc) MD trajectories are used. A Gaussian window is one example
instrument resolution function that can be used can be used where:

.. math::
   :label: eqn-fca12

   W(l \Delta t ) = \frac{1}{2n_{\mathrm{c}} - 1} \frac{\sqrt{2 \pi}}{ \sigma \Delta t} \sum_{m=-(n_{\mathrm{c}}-1)}^{n_{\mathrm{c}}-1} \exp\left[ 2 \pi i \frac{m \Delta \omega}{2n_{\mathrm{c}} - 1} l \Delta t \right] \exp\left[
   -\frac{1}{2}\left(\frac{ m \Delta \omega }{\sigma}\right)^2
   \right]

here :math:`l = -(n_{\mathrm{c}}-1), \ldots, n_{\mathrm{c}}-1` and :math:`\sigma` corresponds to the width of the resolution
function of the Fourier spectrum.
