
This section is dealing with specific types of analysis performed by
MDANSE. If you are not sure where these fit into the general workflow
of data analysis, please read :ref:`workflow-of-analysis`.

Analysis Theory: Scattering
===========================

This section contains background theory for following plugins:

-  :ref:`current-correlation-function`
-  :ref:`dynamic-coherent-structure-factor`
-  :ref:`dynamic-incoherent-structure-factor`
-  :ref:`elastic-incoherent-structure-factor`
-  :ref:`gaussian-dynamic-incoherent-structure-factor`
-  :ref:`neutron-dynamic-total-structure-factor`
-  :ref:`structure-factor-from-scattering-function`

This section discusses plugins used
to calculate neutron spectroscopy observables from the trajectory.
These plugins will be explored in depth in further sections, however,
before that, it is important to understand how MDANSE performs these
analyses. A part of that are the :ref:`param-q-vectors`, which
are used to perform these analyses.

.. _scattering_theory:

Background
''''''''''
**Dynamic structure factor S(q, ω)**: This is a central
concept in neutron scattering experiments. This factor characterizes how
scattering intensity changes with alterations in momentum (:math:`q`) and energy (:math:`\omega`)
during scattering events. It is instrumental in unraveling the atomic and
molecular structures of materials.

**Double differential cross-section**: The dynamic structure factor is closely related to the
double differential cross-section, which is a vital measurement in neutron
scattering. The double differential cross-section, :math:`{\mathrm{d}^{2}{\sigma/\mathit{\mathrm{d}\Omega
\mathrm{d}E}}}`, is defined as the number of
neutrons scattered per unit time into the solid angle interval
:math:`{\left\lbrack {\Omega, {\Omega + \mathrm{d}}\Omega} \right\rbrack}` with an
energy in the interval :math:`{\left\lbrack {E, {E + \mathrm{d}}E} \right\rbrack}`. To make meaningful comparisons, the double differential cross-section
is normalized by :math:`\mathrm{d}\Omega`, :math:`\mathrm{d}E`, and the flux of the incoming neutrons. The relationship
between the double differential cross-section and the dynamic structure factor
is given by:

.. math::
   :label: pfx55

   {{\frac{\mathrm{d}^{2}\sigma}{\mathrm{d}\Omega\mathrm{d}E} = N}\frac{k}{k_{0}}S\left( {q,\omega} \right).}

This equation relates the double differential cross-section, which represents
the number of neutrons scattered per unit time into specific solid angle and
energy intervals, to the dynamic structure factor, :math:`S(q, \omega)`. It includes terms
related to the number of atoms (:math:`N`) and wave numbers of scattered (:math:`k`) and
incident (:math:`k_0`) neutrons.

They are related to the corresponding neutron energies by

.. math::
   :label: pfx56
   
   {E = \hbar^{2}}k^{2}\text{/}2m

\ and

.. math::
   :label: pfx57
   
   {E_{0} = \hbar^{2}}k_{0}^{2}\text{/}2m


These equations relate the neutron energies (:math:`E` and :math:`E_0`) to their respective wave
numbers (:math:`k` and :math:`k_0`) using the mass of the neutron (:math:`m`). They are fundamental for
connecting energy and momentum in neutron scattering.

**Momentum and energy transfer**: These equations below define the
momentum (:math:`q`) and energy (:math:`\omega`) transfer in
units of the reduced Planck constant (:math:`\hbar`) based on the incident and scattered
wave numbers and energies:

.. math::
   :label: pfx58

   {{q = \frac{k_{0} - k}{\hbar}},}

.. math::
   :label: pfx59

   {{\omega = \frac{E_{0} - E}{\hbar}}.}

The modulus of the momentum transfer can be expressed in terms of a scattering
angle, energy transfer, and incident neutron energy.

.. math::
   :label: pfx60

   {{\vert q \vert = \vert k_{0} \vert \sqrt{{2 - \frac{\mathit{\hbar\omega}}{E_{0}} - 2}\cos{\theta\sqrt{1 - \frac{\mathit{\hbar\omega}}{E_{0}}}}}}.}


**Intermediate scattering function F(q, t)**:
This equation defines the dynamic structure factor :math:`S(q, \omega)` as a Fourier
transform of the intermediate scattering function :math:`F(q, t)` with respect to
time :math:`t`. It captures information about the structure and dynamics of the
scattering system [Ref16]_. It can be written as:

.. math::
   :label: pfx61

   {S{\left( {q,\omega} \right) = \frac{1}{2\pi}}{\int\limits_{- \infty}^{+ \infty}\mathrm{d}t \, }\exp\left\lbrack {{- i}\omega t} \right\rbrack F\left( {q,t} \right).}

and

.. math::
   :label: pfx62

   {\text{F}{\left( {q,t} \right) = {\sum\limits_{\alpha,\beta}{\Gamma_{\mathit{\alpha\beta}}\left\langle {\exp\left\lbrack {{- i}q\cdot\hat{R}_{\alpha}(0)} \right\rbrack\exp\left\lbrack {iq\cdot\hat{R}_{\beta}(t)} \right\rbrack} \right\rangle}}},}

.. math::
   :label: pfx63

   {{\Gamma_{\mathit{\alpha\beta}} = \frac{1}{N}}\left\lbrack {\overline{b_{\alpha}}{\overline{b_{\beta}} + \delta_{\mathit{\alpha\beta}}}\left( {\overline{b_{\alpha}^{2}} - {\overline{b_{\alpha}}}^{2}} \right)} \right\rbrack.}

The operators :math:`\hat{R}_{\alpha}(t)`
in Eq. :math:numref:`pfx62` are the position
operators of the nuclei in the sample. The brackets
:math:`\langle\ldots\rangle`
denote a quantum thermal average and the time dependence of the position
operators is defined by the Heisenberg picture. The quantities
:math:`b_{\alpha}` are the scattering lengths of the nuclei
which depend on the isotope and
the relative orientation of the spin of the neutron and the spin of the
scattering nucleus. If the spins of the nuclei and the neutron are not
prepared in a special orientation one can assume a random relative
orientation and that spin and position of the nuclei are uncorrelated.
The overline :math:`\overline{...}` appearing in :math:`{\Gamma_{\mathit{\alpha\beta}}}`
denotes an average over isotopes and relative spin orientations of
neutron and nucleus.

**Coherent and incoherent scattering**:
Usually, one splits the intermediate scattering function and the dynamic
structure factor into their *coherent* and *incoherent* parts which
describe collective and single particle motions, respectively. By defining

.. math::
   :label: pfx65

   {b_{\alpha,\mathrm{coh}}\doteq\overline{b_{\alpha}},}

.. math::
   :label: pfx66

   {b_{\alpha,\mathrm{inc}}\doteq\sqrt{\overline{b_{\alpha}^{2}} - {\overline{b_{\alpha}}}^{2}},}

the coherent and incoherent intermediate scattering functions can be
written. They are expressed as sums over pairs of nuclei, with different
treatments for coherent and incoherent scattering lengths.

.. math::
   :label: pfx67

   {\text{F}_{\text{coh}}{\left( {q,t} \right) = \frac{1}{N}}{\sum\limits_{\alpha,\beta}b_{\alpha,\mathrm{coh}}}b_{\beta,\mathrm{coh}}\left\langle {\exp\left\lbrack {{- i}q\cdot\hat{R}_{\alpha}(0)} \right\rbrack\exp\left\lbrack {iq\cdot\hat{R}_{\beta}(t)} \right\rbrack} \right\rangle,}

.. math::
   :label: pfx68

   {\text{F}_{\text{inc}}{\left( {q,t} \right) = \frac{1}{N}}{\sum\limits_{\alpha}{b_{\alpha,\mathrm{inc}}^{2}\left\langle {\exp\left\lbrack {{- i}q\cdot\hat{R}_{\alpha}(0)} \right\rbrack\exp\left\lbrack {iq\cdot\hat{R}_{\alpha}(t)} \right\rbrack} \right\rangle}}.}

MDANSE introduces the partial terms, this considers the contributions from different species :math:`(I, J)` to the scattering process.

.. math::
   :label: pfx69

   {\text{F}_{\text{coh}}{\left( {q,t} \right) = \sum\limits_{I,J\geq I}^{N_{\mathrm{species}}}}\sqrt{n_{I}n_{J}\omega_{I,\text{coh}}\omega_{J,\text{coh}}}F_{\mathit{IJ},\text{coh}}\left( {q,t} \right),}

.. math::
   :label: pfx70

   {\text{F}_{\text{inc}}{\left( {q,t} \right) = {\sum\limits_{I = 1}^{N_{\mathrm{species}}}{n_{I}\omega_{I,\text{inc}}F_{I,\text{inc}}\left( {q,t} \right)}}}}

where:

.. math::
   :label: pfx71

   {\text{F}_{\mathit{IJ},\text{coh}}{\left( {q,t} \right) = \frac{1}{\sqrt{n_{I}n_{J}}}}{\sum\limits_{\alpha}^{n_{I}}{\sum\limits_{\beta}^{n_{J}}\left\langle {\exp\left\lbrack {{- i}q\cdot\hat{R}_{\alpha}\left( 0 \right)} \right\rbrack\exp\left\lbrack {iq\cdot\hat{R}_{\beta}\left( t \right)} \right\rbrack} \right\rangle}},}

.. math::
   :label: pfx72

   {\text{F}_{I,\text{inc}}{\left( {q,t} \right) = \frac{1}{n_{I}}}{\sum\limits_{\alpha = 1}^{n_{I}}\left\langle {\exp\left\lbrack {{- i}q\cdot\hat{R}_{\alpha}\left( 0 \right)} \right\rbrack\exp\left\lbrack {iq\cdot\hat{R}_{\alpha}\left( t \right)} \right\rbrack} \right\rangle}.}

where :math:`n_I`, :math:`n_J`, :math:`n_{\mathrm{species}}`, :math:`\omega_{I,(\mathrm{coh}/\mathrm{inc})}`
and :math:`\omega_{J,(\mathrm{coh}/\mathrm{inc})}` are defined in Section :ref:`target_CN`. The corresponding dynamic structure factors are obtained by performing
the Fourier transformation defined in Eq. :math:numref:`pfx61`.


**Classical framework and corrections**:
In the classical framework the intermediate scattering functions are
interpreted as classical time correlation functions. The position
operators are replaced by time-dependent vector functions and quantum
thermal averages are replaced by classical *ensemble averages*. It is
well known that this procedure leads to a loss of the universal detailed
balance relation,

.. math::
   :label: pfx74

   {\text{S}{\left( {q,\omega} \right) = \exp}\left\lbrack {\beta\hbar\omega} \right\rbrack\text{S}\left( {{- q}{, - \omega}} \right),}

and also to a loss of all odd moments

.. math::
   :label: pfx75

   {\left\langle \omega^{2{n + 1}} \right\rangle\doteq{\int\limits_{- \infty}^{+ \infty}{\mathrm{d}\omega}}\, \omega^{2{n + 1}}S\left( {q,\omega} \right) \qquad {n = 1,2},\ldots.}

The odd moments vanish since the classical dynamic structure factor is
even in :math:`\omega`, assuming invariance of the scattering process with respect to
reflections in space. The first moment is also universal. For an atomic
liquid, containing only one sort of atoms, it reads

.. math::
   :label: pfx76

   {{\left\langle \omega \right\rangle = \frac{\hbar q^{2}}{2M}},}

where :math:`M` is the mass of the atoms.

**Recoil moment**: Formula :math:numref:`pfx76` shows that the
first moment is given by the average kinetic energy (in units of
:math:`\hbar`) of a particle which receives a momentum transfer
:math:`\hbar q`. Therefore, :math:`\langle\omega\rangle`
is called the *recoil moment*. A number of 'recipes' has been suggested
to correct classical dynamic structure factors for detailed balance and
to describe recoil effects in an approximate way. The most popular one
has been suggested by Schofield [Ref17]_

.. math::
   :label: pfx77

   {{\text{S}\left( {q,\omega} \right)\approx\exp\left\lbrack \frac{\beta\hbar\omega}{2} \right\rbrack}_{}\text{S}_{\mathrm{cl}}\left( {q,\omega} \right)}

One can easily verify that the resulting dynamic structure factor
fulfils the relation of detailed balance. Formally, the correction :math:numref:`pfx77`
is correct to first order in :math:`\hbar`. Therefore, it cannot be used
for large :math:`q`-values which correspond to large momentum transfers
:math:`\hbar q`. This is actually true for all correction
methods which have suggested
so far. For more details we refer to [Ref18]_.


**Static structure factor S(q)**: An important quantity describing structural properties of liquids is the
static structure factor. :math:`S(q)` is an integral involving the
dynamic structure factor or the coherent intermediate scattering function
at zero time delay :math:`t = 0`.

.. math::
   :label: pfx73

   {\text{S}(q)\doteq{\int\limits_{- \infty}^{+ \infty}{\mathrm{d}\omega}}\,\text{S}_{\mathrm{coh}}\left( {q,\omega} \right) = \text{F}_{\mathrm{coh}}\left( {q,0} \right).}

**Total structure factors**: MDANSE computes the partial :math:`S(Q)` as the Fourier transform of the
partial pair distribution function :math:`g(r)`, corresponding to the Faber-Ziman definition:

.. math::
   :label: pfx78
   
   {S_{\alpha\beta}(Q) = 1 + \frac{4\pi\rho_0}{Q}\int\limits_{0}^{\infty}{\mathrm{d}r \, r \sin(Qr) \left\lbrack {g_{\alpha\beta}}(r)-1 \right\rbrack}}

The total :math:`S(Q)` is computed as a weighted sum similar to the one used for
the total :math:`g(r)`. In the case of the analysis 'X-ray Static structure
factor', the :math:`Q`-dependence of the atomic form factors is taken into
account in this weighted sum.

**X-ray observable normalization**: Soper has provided experimental data (table 4 in *ISRN Physical
Chemistry*, 279463 (2013), given in file soper13_fx.dat). Here a source
of confusion is that the data can be normalized in different ways (see
Soper's paper). Using the normalization II in that reference we have
that:

.. math::
   :label: pfx79
   
    D_{x}{(Q) = \frac{\sum\limits_{\mathit{\alpha\beta}\geq\alpha}{\left( {2 - \delta_{\mathit{\alpha\beta}}} \right) c_{\alpha}c_{\beta}f_{\alpha}{(Q)}f_{\beta}{(Q)}\left\lbrack {S_{\mathit{\alpha\beta}}{(Q) - 1}} \right\rbrack}}{\sum\limits_{\alpha}{c_{\alpha}f_{\alpha}^{2}{(Q)}}} = \left\lbrack {S{(Q) - 1}} \right\rbrack}\frac{\sum\limits_{\mathit{\alpha\beta}}{c_{\alpha}c_{\beta}f_{\alpha}{(Q)}f_{\beta}{(Q)}}}{\sum\limits_{\alpha}{c_{\alpha}f_{\alpha}^{2}{(Q)}}}

Where :math:`S(Q)` would be the static structure factor (going to :math:`1` at large :math:`Q`)
computed by MDANSE. Therefore, even after using MDANSE we should
recalculate the x-ray observable using the atomic factors.

.. _current-correlation-function:

Current Correlation Function
''''''''''''''''''''''''''''

The current correlation functions :math:`C_{\alpha\beta}(q, t)` and :math:`C_{\alpha\beta}(q, \omega)`
are closely related to the intermediate scattering function :math:`F(q, t)`
and the dynamics structure factor :math:`S(q, \omega)`. The intermediate
scattering function :math:`F(q, t)` is a correlation function of the Fourier components of
particle density whereas the current correlation function :math:`C_{\alpha\beta}(q, t)`
are a correlation function of the Fourier components of the particle current:

.. math::

    C_{\alpha\beta}(q, t) = \frac{1}{N} \langle j_{\alpha}(q, t) j_{\beta}(-q, 0) \rangle \qquad\qquad j_{\alpha}(q, t) = \sum_{l} v_{l\alpha}(t) \exp(iq\cdot r_l(t))

where :math:`\alpha, \beta = x, y` or :math:`z`. The particle currents
can be projected onto longitudinal and transverse components of the
:math:`q`-vector. The longitudinal and transverse particle current are:

.. math::

    j_{\mathrm{L}}(q, t) = \sum_{l} (v_{l\alpha}(t) \cdot \hat{q})\hat{q} \, \exp(iq\cdot r_l(t))

.. math::

    j_{\mathrm{T}}(q, t) = \sum_{l} [v_{l\alpha}(t) - (v_{l\alpha}(t) \cdot \hat{q})\hat{q}] \, \exp(iq\cdot r_l(t))

where :math:`\hat{q}` are unit vectors of :math:`q`. For isotropic systems,
the longitudinal and transverse particle current are uncorrelated and the
current correlation function tensor :math:`C_{\alpha\beta}(q, t)`
will depend only on two independent components so that we can write:

.. math::

    C_{\alpha\beta}(q, t) = \hat{q}_{\alpha}\hat{q}_{\beta} C_{\mathrm{L}}(q, t) + (\delta_{\alpha\beta} - \hat{q}_{\alpha}\hat{q}_{\beta}) C_{\mathrm{T}}(q, t)

where the longitudinal and transverse current correlation functions are.

.. math::

    C_{\mathrm{L}}(q, t) = \frac{1}{N} \langle j_{\mathrm{L}}(q, t) \cdot j_{\mathrm{L}}(-q, 0) \rangle

.. math::

    C_{\mathrm{T}}(q, t) = \frac{1}{N} \langle j_{\mathrm{T}}(q, t) \cdot j_{\mathrm{T}}(-q, 0) \rangle


From the continuity equation we can obtain a relation between the
longitudinal current correlation and the dynamic structure factor.

.. math::

    C_{\mathrm{L}}(q, \omega) = \frac{\omega^2}{Q^2} S(q, \omega)

.. _dynamic-coherent-structure-factor:

Dynamic Coherent Structure Factor
'''''''''''''''''''''''''''''''''
In this analysis, MDANSE proceeds in two steps. First, it computes the partial
and total intermediate coherent scattering function using equation
:math:numref:`pfx69`. Then, the partial and total dynamic coherent structure
factors are obtained by performing the Fourier Transformation, defined in Eq.
:math:numref:`pfx61`, respectively on the total and partial intermediate
coherent scattering functions.

**Coherent intermediate scattering function**:

.. math::
   :label: pfx80
   
   {{F}_{\text{coh}}\left( {q, t} \right)\doteq{\sum\limits_{{I = 1},J\geq I}^{N_{\mathrm{species}}}\sqrt{n_{I}n_{J}\omega_{I,\text{com}}\omega_{I,\text{com}}}}{\left\langle {\rho_{I}\left( {{-q},0} \right)\rho_{J}\left( {q, t} \right)} \right\rangle},}

where :math:`N_{\mathrm{species}}` is the number of selected species, :math:`n_{I}`
is the number of atoms , :math:`\omega_{I}` is the weight and
:math:`{\rho_{I}( {q, t})}` is the Fourier components of the particle density of species
:math:`I`

.. math::
   :label: pfx83

   {\rho_{I}{\left( {q, t} \right) = \sum\limits_{\alpha}^{n_{I}}}\exp\left\lbrack {\mathit{iq}\cdot R_{\alpha}\left( t \right)} \right\rbrack.}

In MDANSE, the DCSF is averaged over :math:`q`-vectors having *approximately* the same
modulus so that

.. math::
   :label: pfx80

   F_{\text{coh}}\left( {Q_{m}, t} \right) = \overline{F_{\text{coh}}\left( {q_m, t} \right)}^{q_m}

where :math:`Q_{m} \approx \vert q_{m} \vert`.

**Reciprocal lattice q-vectors**: Let
:math:`b_1`, :math:`b_2`, :math:`b_3` be the basis vectors
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

**Negative coherent scattering lengths**: The :math:`q`-vectors can be generated isotropically, anisotropically or along
user-defined directions. The :math:`\sqrt{\omega_{I}}` may be negative
if they represent normalized coherent scattering
lengths, i.e.

.. math::
   :label: pfx89

   {{\sqrt{\omega_{I}} = \frac{b_{I,\text{coh}}}{\sqrt{\sum\limits_{I = 1}^{N_{\mathrm{species}}}{n_{I}b_{I,\text{coh}}^{2}}}}}.}

Defines the use of negative coherent scattering lengths for hydrogenous materials.
Negative coherent scattering lengths occur in hydrogenous materials
since :math:`b_{\mathrm{H},\mathrm{coh}}` is negative [Ref20]_. When the
default value of weights (``b_coherent``) is chosen for this
analysis, the result will correspond to that of the equation :math:numref:`ntdsf-eq6`
from the :ref:`neutron-dynamic-total-structure-factor`.

.. _dynamic-incoherent-structure-factor:

Dynamic Incoherent Structure Factor
'''''''''''''''''''''''''''''''''''
                      
In this analysis, MDANSE proceeds in two steps. First, it computes
the partial and total intermediate incoherent scattering function
:math:`F_{\mathrm{inc}}(q, t)` using equation :math:numref:`pfx70`. Then, the
partial and total dynamic incoherent structure factors are obtained by
performing the Fourier Transformation, defined in Eq. :math:numref:`pfx61`,
respectively on the total and partial intermediate incoherent
scattering function.

**Incoherent intermediate scattering function**:

.. math::
   :label: pfx90

   {\text{F}_{\text{inc}}\left( {q, t} \right)\doteq{\sum\limits_{I = 1}^{N_{\mathrm{species}}}{n_{I}\omega_{I,\text{inc}}}}\text{F}_{I,\text{inc}}\left( {q, t} \right)}

where :math:`N_{\mathrm{species}}`
is the number of selected species, :math:`n_I` the
number of atoms of species :math:`n_I`, :math:`\omega_{I}` the weight for species :math:`I`
(see Section :ref:`target_CN` for more details) and :math:`{F_{I,\text{inc}}\left( {q, t} \right)}`
is defined as:

.. math::
   :label: pfx92

   {\text{F}_{I,\mathrm{inc}}{\left( {q_{m}, t} \right) = \sum\limits_{\alpha = 1}^{n_{I}}}{\left\langle {\exp\left\lbrack {{-i}q\cdot R_{\alpha}(0)} \right\rbrack\exp\left\lbrack {iq\cdot R_{\alpha}(t)} \right\rbrack} \right\rangle}.}

In MDANSE, the DISF is averaged over :math:`q`-vectors having *approximately* the same
modulus so that

.. math::
   :label: pfx80

   F_{\text{inc}}\left( {Q_{m}, t} \right) = \overline{F_{\text{inc}}\left( {q_m, t} \right)}^{q_m}

where :math:`Q_{m} \approx \vert q_{m} \vert`.

.. _elastic-incoherent-structure-factor:

Elastic Incoherent Structure Factor
'''''''''''''''''''''''''''''''''''

The elastic incoherent structure
factor (EISF) appears as the amplitude of the *elastic* line in the neutron
scattering spectrum. Elastic scattering is only present for systems in which the
atomic motion is confined in space, as for solids. To understand which information
is contained in the EISF we consider for simplicity a system where only one
sort of atoms is visible to the neutrons. To a very good approximation this is
the case for all systems containing a large amount of hydrogen atoms, as biological
systems. Incoherent scattering from hydrogen dominates by far all other
contributions.

**The van Hove self-correlation function**: The EISF is defined as the
limit of the incoherent intermediate scattering function for infinite time

.. math::
   :label: pfx97

   {\mathrm{EISF}(q)\doteq\lim\limits_{t\rightarrow\infty}\text{F}_{\mathrm{inc}}\left( {q,t} \right).}

Using the above definition of the EISF one can decompose the incoherent
intermediate scattering function as follows:

.. math::
   :label: pfx98

   {\text{F}_{\text{inc}}{\left( {q,t} \right) = \mathrm{EISF}}{(q) + \text{F}_{\text{inc}}^{'}}\left( {q,t} \right),}

where :math:`F^{'}_{\mathrm{inc}}(q, t)` decays to zero for infinite time. Taking
now the Fourier transform it follows immediately that

.. math::
   :label: pfx99

   {\text{S}_{\text{inc}}{\left( {q,\omega} \right) = \mathrm{EISF}}(q)\delta{(\omega) + \text{S}_{\text{inc}}^{'}}\left( {q,\omega} \right).}

The EISF appears as the amplitude of the *elastic* line in the neutron
scattering spectrum. Elastic scattering is only present for systems in
which the atomic motion is confined in space, as for solids. To
understand which information is contained in the EISF we consider for
simplicity a system where only one sort of atoms is visible to the
neutrons. To a very good approximation this is the case for all systems
containing a large amount of hydrogen atoms, as biological systems.
Incoherent scattering from hydrogen dominates by far all other
contributions. Using the definition of the van Hove self-correlation
function :math:`G_{\mathrm{s}}(r, t)` [Ref20]_,

.. math::
   :label: pfx100

   {b_{\text{inc}}^{2}G_{\mathrm{s}}\left( {r,t} \right)\doteq\frac{1}{2\pi^{3}}{\int \mathrm{d}^{3}}q \, \exp\left\lbrack {{- i}q\cdot r} \right\rbrack\text{F}_{\mathrm{inc}}\left( {q,t} \right),}

which can be interpreted as the conditional probability to find a tagged
particle at the position :math:`r` at time :math:`t`, given it started at :math:`r = 0`,
one can write:

.. math::
   :label: pfx101

   {\mathrm{EISF}(q) = b_{\text{inc}}^{2}{\int \mathrm{d}^{3}}r \, \exp\left\lbrack {\mathit{iq}\cdot r} \right\rbrack G_{\mathrm{s}}\left( {r,{t = \infty}} \right).}

The EISF gives the sampling distribution of the points in space in the
limit of infinite time. In a real experiment this means times longer
than the time which is observable with a given instrument. The EISF
vanishes for all systems in which the particles can access an infinite
volume since :math:`G_{\mathrm{s}}(r, t)` approaches :math:`1/V` for large times. This is
the case for molecules in liquids and gases.

**EISF computation**: For computational purposes it is convenient to use the following
representation of the EISF [Ref21]_:

.. math::
   :label: pfx102

   {\mathrm{EISF}{(q) = {\sum\limits_{I = 1}^{N_{\mathrm{species}}}{n_{I}\omega_{I,\text{inc}}\mathrm{EISF}_{I}(q)}}}}

where :math:`N_{\mathrm{species}}` is the number of selected species, :math:`n_I`
the number of atoms of species :math:`I`, :math:`\omega_{I,\mathrm{inc}}` the weight for
species :math:`I` (see Section :ref:`target_CN` for more details) and for each species the
following expression for the elastic incoherent scattering function is

.. math::
   :label: pfx103

   {\mathrm{EISF}_{I}{(q) = \frac{1}{n_{I}}}{\sum\limits_{\alpha}^{n_{I}}\left\langle {|{\exp\left\lbrack {\mathit{iq}\cdot R_{\alpha}} \right\rbrack\left. {} \right|^{2}}} \right\rangle}.}

This expression is derived from definition :math:numref:`pfx97`
of the EISF and expression :math:numref:`pfx70` for the
intermediate scattering function, using that for infinite time the
relation

.. math::
   :label: pfx104
   
   \lim\limits_{t\rightarrow\infty}{\left\langle {\exp\left\lbrack {{- \mathit{iq}}\cdot R_{\alpha}(0)} \right\rbrack\exp\left\lbrack {\mathit{iq}\cdot R_{\alpha}(t)} \right\rbrack} \right\rangle = \left\langle {|{\exp\left\lbrack {\mathit{iq}\cdot R_{\alpha}} \right\rbrack\left. {} \right|^{2}}} \right\rangle}

holds. In this way the computation of the EISF is reduced to the
computation of a static thermal average. We remark at this point that
the length of the MD trajectory from which the EISF is computed
should be long enough to allow for a representative sampling of the
conformational space.

In MDANSE, the EISF is averaged over :math:`q`-vectors having *approximately* the same
modulus so that

.. math::
   :label: pfx80

   \mathrm{EISF}\left( {Q_{m}, t} \right) = \overline{\mathrm{EISF}\left( {q_m, t} \right)}^{q_m}

where :math:`Q_{m} \approx \vert q_{m} \vert`.

.. _gaussian-dynamic-incoherent-structure-factor:

Gaussian Dynamic Incoherent Structure Factor
''''''''''''''''''''''''''''''''''''''''''''

**Cumulant expansion**: The MSD can be related to the incoherent intermediate scattering
function via the cumulant expansion [Ref11]_, [Ref22]_

.. math::
   :label: pfx108

   {\text{F}_{\text{inc}}{\left( {q,t} \right) = {\sum\limits_{I = 1}^{N_{\mathrm{species}}}{n_{I}\omega_{I,\text{inc}}}}}\text{F}_{I,\text{inc}}\left( {q,t} \right)}

.. math::
   :label: pfx109

   {\text{F}_{I,\text{inc}}{\left( {q,t} \right) = \frac{1}{n_{I}}}\sum\limits_{\alpha}^{n_{I}}\exp\left\lbrack {{- q^{2}}\rho_{\alpha,1}{(t) + q^{4}}\rho_{\alpha,2}(t) + \ldots} \right\rbrack}

where :math:`N_{\mathrm{species}}` is the number of selected species, :math:`n_I`
the number of atoms of species :math:`I`, :math:`\omega_{I,\mathrm{inc}}` the weight for
species :math:`I` (see Section :ref:`target_CN` for more details). The cumulants
:math:`\rho_{\alpha,k}(t)` are identified as

.. math::

   {\rho_{\alpha,1}{(t) = \frac{1}{2!}\left\langle {d_{\alpha}^{2}\left( {t;n_{q}} \right)} \right\rangle}}

.. math::

   {\rho_{\alpha,2}{(t) = \frac{1}{4!}}\left\lbrack {{\left\langle {d_{\alpha}^{4}\left( {t;n_{q}} \right)} \right\rangle - 3}\left\langle {d_{\alpha}^{2}\left( {t;n_{q}} \right)} \right\rangle^{2}} \right\rbrack}

.. math::
   :label: pfx112
   {\vdots}

**Gaussian approximation**: The vector :math:`n_q` is the unit vector
in the direction of :math:`q`. In the Gaussian
approximation the above expansion is truncated after the
:math:`q^2`-term.

.. math::

   {\text{F}_{I,\text{inc}}\left( {q,t} \right) \approx \text{F}_{I,\text{inc}}^{\mathrm{g}}{\left( {q,t} \right) = \frac{1}{n_{I}}}\sum\limits_{\alpha}^{n_{I}}\exp\left\lbrack {{- q^{2}}\rho_{\alpha,1}(t)} \right\rbrack}

For certain model systems like the ideal gas, the
harmonic oscillator, and a particle undergoing Einstein diffusion, this
is exact. For these systems the incoherent intermediate scattering
function is completely determined by the MSD.

**Intermediate scattering function**: For an isotropic system, the
following expression for the intermediate scattering function is obtained

.. math::
   :label: pfx114

   {\text{F}_{I,\alpha,\text{inc}}^{\mathrm{g}}{\left( {Q, t} \right) = \frac{1}{n_{I}}}\sum\limits_{\alpha}^{n_{I}}\exp\left\lbrack {\frac{- Q^{2}}{6}\Delta_{\alpha}^{2}\left( { t} \right)} \right\rbrack}

.. math::
   :label: pfx115

   {\text{F}_{I,\alpha,\text{inc}}^{\mathrm{g}}{\left( {Q, t} \right) = \frac{1}{n_{I}}}\sum\limits_{\alpha}^{n_{I}}\exp\left\lbrack {\frac{- Q^{2}}{2}\Delta_{\alpha}^{2}\left( { t;n} \right)} \right\rbrack}

where the quantities :math:`\Delta_{\alpha}^{2}(t)`
and :math:`\Delta_{\alpha}^{2}\left( {t;n} \right)` are the mean-square
displacements, defined in Equations :math:numref:`pfx14`
and :math:numref:`pfx15`, respectively.
They are computed by using the algorithm described in the :ref:`analysis-msd` section.
It should be noted that the computation of the
intermediate scattering function in the Gaussian approximation is much
'cheaper' than the computation of the full intermediate scattering
function, :math:`F_{\mathrm{inc}}(q, t)`, since no averaging over different
:math:`q`-vectors needs to be performed.

.. _neutron-dynamic-total-structure-factor:

Neutron Dynamic Total Structure Factor
''''''''''''''''''''''''''''''''''''''

**Partial coherent intermediate scattering functions and dynamic structure factors**:
This is a combination of the dynamic coherent and the dynamic incoherent
structure factors. It is a fully neutron-specific analysis, where the
coherent part of the intermediate scattering function is calculated
using the atomic coherent neutron scattering lengths ``b_coherent`` and
the incoherent one is calculated using the square of the atomic
incoherent neutron scattering lengths ``b_incoherent2``. Therefore, in
this analysis the weights option is not available.

The partial coherent intermediate scattering functions
:math:`I_{\alpha\beta}^{\mathrm{coh}}(q,t)` and their corresponding Fourier
transforms giving the partial coherent dynamic structure factors,
:math:`S_{\alpha\beta}^{\mathrm{coh}}(q,\omega)` are calculated exactly in the
same way as in the DCSF analysis, i.e.:

.. math::
   :label: ntdsf-eq1
   
   I_{\alpha\beta}^{\mathrm{coh}}(q,t) = \frac{1}{\sqrt{N_{\alpha}N_{\beta}}}\sum_{i \in \alpha,j \in \beta}^{N_{\alpha},N_{\beta}}\left\langle \exp\left[- iq \cdot r_{i}(0)\right] \exp \left[iq \cdot r_{j}(t)\right] \right\rangle

where :math:`\alpha` and :math:`\beta` refer to the chemical elements,
:math:`N_{\alpha}` and :math:`N_{\beta}` are the respective number of
atoms of each type, :math:`i` and :math:`j` are two specific atoms of
type :math:`\alpha` and :math:`\beta`, respectively, and
:math:`\mathbf{r}_{i}(0)` and :math:`\mathbf{r}_{j}(t)` are their
positions at the times :math:`0` and :math:`t`, respectively. Similarly,
the partial incoherent intermediate scattering functions
:math:`I_{\alpha}^{\mathrm{inc}}(q,t)` and the partial incoherent dynamic
structure factors :math:`S_{\alpha}^{\mathrm{inc}}(q,\omega)` are obtained as in
the DISF analysis:

.. math::
   :label: ntdsf-eq2
   
   I_{\alpha}^{\mathrm{inc}}(q,t) = \frac{1}{N_{\alpha}}\sum_{i \in \alpha}^{N_{\alpha}}\left\langle \exp\left[- iq \cdot r_{i}(0)\right] \exp \left[iq \cdot r_{i}(t)\right] \right\rangle

**Combination of partial contributions**: The main difference between
this analysis and the DCSF and DISF
analyses, apart from the fact that the coherent and incoherent
contributions are calculated simultaneously, is the way the different
partial contributions are combined. In this analysis the total
incoherent, total coherent and total (coherent + incoherent) signals are
calculated as:

.. math::
   :label: ntdsf-eq3
   
   I^{\mathrm{inc}}(q,t) = \sum_{\alpha}^{N_{\alpha}}{c_{\alpha}b_{\alpha,\text{inc}}^{2}}I_{\alpha}^{\mathrm{inc}}(q,t)

.. math::
   :label: ntdsf-eq4
   
   I^{\mathrm{coh}}(q,t) = \sum_{\alpha,\beta}^{N_{\alpha},N_{\beta}}{\sqrt{c_{\alpha}c_{\beta}}b_{\alpha,\text{coh}}b_{\beta,\text{coh}}I_{\alpha\beta}^{\mathrm{coh}}(q,t)}

.. math::
   :label: ntdsf-eq5
   
   I^{\mathrm{tot}}(q,t) = I^{\mathrm{inc}}(q,t) + I^{\mathrm{coh}}(q,t)

where :math:`c_{\alpha} = N_{\alpha} / N` and
:math:`c_{\beta} =  N_{\beta} / N` are the concentration numbers
for elements :math:`\alpha` and :math:`\beta`, respectively. These
expressions correspond to the formalism and equations given in
[Ref47]_ - Chapter 1: “An introduction to neutron scattering” .

**Units conversion**:
As in the MDANSE database the coherent and incoherent neutron scattering
lengths are given in Å, the total intermediate scattering functions
above will be given in Å\ :sup:`2`/sterad/atom. Therefore, multiplying
the output from MDANSE by a factor 10\ :sup:`8` we can obtain these
neutron observables in barn/sterad/atom and compare them directly to the
experimental results (assuming the later have been properly normalized
and presented in absolute units).

On the other hand, the DISF and DCSF analyses use the standard weight
normalization procedure implemented in MDANSE (see :ref:`param-normalize`).
Therefore the total coherent intermediate scattering function
returned by the DCSF analysis is (assuming that the chosen weights are
``b_coherent``):

.. math::
   :label: ntdsf-eq6
   
   I^{\mathrm{coh}}(q,t) = \frac{\sum_{\alpha\beta}^{n}{c_{\alpha}c_{\beta}b_{\alpha,\mathrm{coh}}b_{\beta,\mathrm{coh}}I_{\alpha\beta}^{\mathrm{coh}}(q,t)}}{\sum_{\alpha\beta}^{n}{c_{\alpha}c_{\beta}b_{\alpha,\mathrm{coh}}b_{\beta,\mathrm{coh}}}}

and the incoherent intermedicate scattering function given by the DISF
analysis is (assuming that the chosen weights are ``b_incoherent2``):

.. math::
   :label: ntdsf-eq7
   
   I^{\mathrm{inc}}(q,t) = \frac{\sum_{\alpha}^{n}{c_{\alpha}b_{\alpha,\mathrm{inc}}^{2}I_{\alpha}^{\mathrm{inc}}(q,t)}}{\sum_{\alpha}^{n}{c_{\alpha}b_{\alpha,\mathrm{inc}}^{2}}}

Naturally, similar expressions apply to the dynamic structure factors,
:math:`S_{\alpha\beta}^{\mathrm{coh}}(q,\omega)` and
:math:`S_{\alpha}^{\mathrm{inc}}(q,\omega)`.

.. _structure-factor-from-scattering-function:

Structure Factor From Scattering Function
'''''''''''''''''''''''''''''''''''''''''
The structure factor from scattering function is a concept used in
scientific research, particularly in the field of neutron scattering
experiments. It relates to how particles or atoms within a material
contribute to its overall structural properties based on their scattering
behavior. This concept provides valuable insights into the material's
internal structure, dynamics, and interactions. Researchers use the Structure
Factor From Scattering Function to better understand the atomic-level details
of materials, which has applications in diverse areas, including materials
science and condensed matter physics
