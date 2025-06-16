
.. _grouping:


Results Grouping
================

Atom Grouping
^^^^^^^^^^^^^

With the default settings, MDANSE will group results by atom so that results are
summed for all atoms of each type in the system and those results are
summed using the :ref:`weighting-scheme` to obtain the total. Lets
consider a system of water and ethanol. With atom grouping,
the partial coherent intermediate scattering
functions are obtained for each unique pair of atoms in the system:
CC, CH, CO, HH, HO, and OO. For example, the C and H partial coherent
intermediate scattering functions is

.. math::
   :label: atompartial1

   F_{\text{coh},\text{CH}}{(\mathbf{q},t) =  \frac{1}{N \sqrt{c_{\text{C}}c_{\text{H}}}}}{\sum\limits_{j \in \text{C}}{\sum\limits_{k \in \text{H}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{j}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{k}\left( t \right)} \right\rbrack} \right\rangle}},

so the summations run over all carbon and hydrogen atoms in the system.
The total is a weighted sum of all partial components

.. math::
   :label: atompartial2

    F_{\text{coh}}(\mathbf{q},t) = W_{\text{CC}} F_{\text{coh}, \text{CC}}(\mathbf{q},t) + W_{\text{CH}} F_{\text{coh}, \text{CH}}(\mathbf{q},t) + \cdots

Similarly the partial incoherent intermediate
scattering functions are obtained for each unique atom in the system e.g.
C, H, and O. For exmaple, the C partial incoherent
intermediate scattering functions is

.. math::
   :label: atompartial3

   F_{\text{inc},\text{C}}{(\mathbf{q},t ) = \frac{1}{Nc_{\text{C}}}}{\sum\limits_{j \in \text{C}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{j}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle}

so the summation run over all carbon atoms and the total is a weight sum
of all partial components

.. math::
   :label: atompartial2

    F_{\text{inc}}(\mathbf{q},t) = W_{\text{C}} F_{\text{inc}, \text{C}}(\mathbf{q},t) + W_{\text{H}} F_{\text{inc}, \text{H}}(\mathbf{q},t) + \cdots,

see :ref:`weighting-scheme` for further details.

Molecule Grouping
^^^^^^^^^^^^^^^^^

MDANSE can also group results by molecule so that they are summed for
all atoms of each type on each type of molecule, again we will consider
a system of water and ethanol. With molecule grouping, the partial
coherent intermediate scattering functions are obtained for each
unique pair of molecules and their atoms types:
[EtOH][EtOH]_CC, [EtOH][EtOH]_CH, [EtOH][H2O]_HO, and etc. Where
[EtOH][H2O]_HO are all pairs of H and O atoms where H are the hydrogen
atoms in ethanol and O are the oxygen atoms in water. The [EtOH][H2O]_HO
partial coherent intermediate scattering functions is

.. math::
   :label: moleculepartial1

   F_{\text{coh},\text{HO}}^{[\text{EtOH}][\text{H2O}]}{(\mathbf{q},t) =  \frac{1}{N \sqrt{c_{\text{H}}^{\text{EtOH}}c_{\text{O}}^{\text{H2O}}}}}{\sum\limits_{j \in (\text{EtOH}\, \cap\, H)}{\sum\limits_{k \in (\text{H2O}\, \cap\, O)}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{j}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{k}\left( t \right)} \right\rbrack} \right\rangle}},

where :math:`c_{\text{H}}^{\text{EtOH}} = N_{\text{H}}^{\text{EtOH}} / N` and
:math:`c_{\text{O}}^{\text{H2O}} = N_{\text{O}}^{\text{H2O}} / N`.
:math:`N_{\text{H}}^{\text{EtOH}}` and :math:`N_{\text{O}}^{\text{H2O}}` are
the total number of atoms atoms of hydrogen in ethanol and oxygen in water respectively,
and :math:`N` is the total number of atom in the system. The molecular
coherent intermediate scattering functions between ethanol and water is
proportional to the weighted sum of the partial terms

.. math::
   :label: moleculepartial2

   \sqrt{c_{\text{EtOH}} c_{\text{H2O}}}F_{\text{coh}}^{[\text{EtOH}][\text{H2O}]}(\mathbf{q},t) =  &W_{\text{CH}}^{[\text{EtOH}][\text{H2O}]} F_{\text{coh},\text{CH}}^{[\text{EtOH}][\text{H2O}]}(\mathbf{q},t) \\\\
    &+ W_{\text{CO}}^{[\text{EtOH}][\text{H2O}]} F_{\text{coh},\text{CO}}^{[\text{EtOH}][\text{H2O}]}(\mathbf{q},t) + \cdots

where :math:`c_{\text{EtOH}} = N_{\text{EtOH}} / N` and
:math:`c_{\text{H2O}} = N_{\text{H2O}} / N`. :math:`N_{\text{EtOH}}`
and :math:`N_{\text{H2O}}` are the total number of atoms in ethanol and
water respectively and for coherent scattering lengths

.. math::
   :label: moleculepartial3

   W_{\text{HO}}^{[\text{EtOH}][\text{H2O}]} = 2\frac{\sqrt{c_{\text{H}}^{\text{EtOH}}c_{\text{O}}^{\text{H2O}}} b_{\mathrm{coh},\text{H}}b_{\mathrm{coh},\text{O}}}{c_{\text{C}}^{\text{EtOH}}c_{\text{C}}^{\text{EtOH}}  b_{\mathrm{coh},\text{C}}b_{\mathrm{coh},\text{C}} + \cdots}.

The total coherent intermediate scattering functions
is a weighted sum of all molecular terms

.. math::
   :label: moleculepartial4

   F_{\text{coh}}(\mathbf{q},t) = \sqrt{c_{\text{EtOH}} c_{\text{EtOH}}}F_{\text{coh}}^{[\text{EtOH}][\text{EtOH}]}(\mathbf{q},t) + \sqrt{c_{\text{EtOH}} c_{\text{H2O}}}F_{\text{coh}}^{[\text{EtOH}][\text{H2O}]}(\mathbf{q},t) + \cdots

where here the weights are simply the square roots of the product of the
atom concentrations.

Similarly the partial incoherent intermediate
scattering functions are obtained for each molecule and its atom types:
[EtOH]_C, [EtOH]_H, [EtOH]_O, [H2O]_H, and [H2O]_O. For example,
the partial incoherent intermediate scattering functions for ethanols
carbon atoms is

.. math::
   :label: moleculepartial5

   F_{\text{inc},\text{C}}^{\text{EtOH}}{(\mathbf{q},t ) = \frac{1}{Nc_{\text{C}}^{\text{EtOH}}}}{\sum\limits_{j \in (\text{C}\, \cap \, \text{EtOH})}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{j}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle},

the molecular incoherent intermediate scattering functions for ethanol
is

.. math::
   :label: moleculepartial6

   c_{\text{EtOH}}F_{\text{inc},\text{EtOH}}(\mathbf{q},t ) = W^{\text{EtOH}}_{\text{C}} F_{\text{inc},\text{C}}^{\text{EtOH}}(\mathbf{q},t ) + W^{\text{EtOH}}_{\text{H}}F_{\text{inc},\text{H}}^{\text{EtOH}}(\mathbf{q},t ) + \cdots.

With incoherent scattering length weights

.. math::
   :label: moleculepartial7

   W_{\text{C}}^{\text{EtOH}} = \frac{c_{\text{C}}^{\text{EtOH}}b_{\text{inc},C}^{2}}{c_{\text{C}}^{\text{EtOH}}b_{\text{inc},C}^{2} + \cdots}.

The total incoherent intermediate scattering functions
is a weighted sum of all molecular terms

.. math::
    :label:

    F_{\text{inc}}(\mathbf{q},t ) = c_{\text{EtOH}}F_{\text{inc},\text{EtOH}}(\mathbf{q},t ) + c_{\text{H2O}}F_{\text{inc},\text{H2O}}(\mathbf{q},t )

where here the weight are the atom contrations of the atoms in ethanol
and water.

Scaling Factors
~~~~~~~~~~~~~~~
Similarly to the atom grouping, the total results with molecule grouping
are a weighted sum of atomic or molecular terms. In MDANSE either, scaled or
unscaled results can be plotted and may be more useful for the specific
results that has been calculated. For example




Other Groupings
^^^^^^^^^^^^^^^

