Weighting Scheme
================

In MDANSE, most properties are split by atom-type
and the total results is a sum of these atom-type based
properties. For example, the coherent and incoherent intermediate
scattering functions are

.. math::
   :label: ws1

   F_{\text{coh},\alpha\beta}{(\mathbf{q},t) =  \frac{W_{\alpha\beta}}{N c_{\alpha}c_{\beta}}}{\sum\limits_{i}^{N_{\alpha}}{\sum\limits_{j}^{N_{\beta}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle}},

.. math::
   :label: ws2

   F_{\text{inc},\alpha}{(\mathbf{q},t ) = \frac{W_{\alpha}}{Nc_{\alpha}}}{\sum\limits_{i}^{N_{\alpha}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{i}\left( t \right)} \right\rbrack} \right\rangle}

where :math:`\alpha` and :math:`\beta` are the atom-types.
:math:`W_{\alpha\beta}` and :math:`W_{\alpha}` are the weights of the
atom-type pairs :math:`\alpha\beta` and the atom type :math:`\alpha`.
:math:`c_{\alpha} = N_{\alpha} / N` and :math:`c_{\beta} = N_{\beta} / N` are the concentrations of atoms of
atom-types :math:`\alpha` and :math:`\beta` and :math: N_{\alpha}`,
:math:`N_{\beta}`, and :math:`N` are the :math:`\alpha`, :math:`\beta`,
and the total number of atoms. The total is now a sum of these atomic terms

.. math::
   :label: ws3

    F_{\text{coh}}(\mathbf{q},t) = \sum_{\alpha}\sum_{\beta \geq \alpha} \text{F}_{\text{coh},\alpha\beta}(\mathbf{q},t),

.. math::
   :label: ws4

    F_{\text{inc}}(\mathbf{q},t) = \sum_{\alpha} \text{F}_{\text{inc},\alpha}(\mathbf{q},t).

Note that for summation involving two atom-types only the unique pairs
are summed up. This is because in MDANSE the off-diagonal weight
terms are symmetrised so that we assumed that
:math:`F_{\text{coh},\alpha\beta} = F_{\text{coh},\beta\alpha}`.

.. _water-dos-weighted:

.. figure:: ./Pictures/water_dos_weighted.png
   :align: center
   :width: 11.748cm
   :height: 9.393cm

   The total and atomic DOS of water, atomic DOS are **weighted** so that the
   sum of atomic DOS equals to the total.

The atom properties can also be scaled without the weights

.. math::
   :label: ws5

   \mathcal{F}_{\text{coh},\alpha\beta}{(\mathbf{q},t) = \frac{1}{N c_{\alpha} c_{\beta}}}{\sum\limits_{i}^{N_{\alpha}}{\sum\limits_{j}^{N_{\beta}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle}},

.. math::
   :label: ws6

   \mathcal{F}_{\text{inc},\alpha}{(\mathbf{q},t ) = \frac{1}{N c_{\alpha}}}{\sum\limits_{i}^{N_{\alpha}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{i}\left( t \right)} \right\rbrack} \right\rangle}

so the total will now be a weighted sum of these atomic terms

.. math::
   :label: ws7

    F_{\text{coh}}(\mathbf{q},t) = \sum_{\alpha}\sum_{\beta \geq \alpha} W_{\alpha\beta} \mathcal{F}_{\text{coh},\alpha\beta}(\mathbf{q},t),

.. math::
   :label: ws8

    F_{\text{inc}}(\mathbf{q},t) = \sum_{\alpha} W_{\alpha} \mathcal{F}_{\text{inc},\alpha}(\mathbf{q},t).

In the MDANSE_GUI you have the option to plot either weighted (e.g. :math:`F_{\text{coh},\alpha\beta}`
and :math:`F_{\text{inc},\alpha}`) or unweighted (e.g. :math:`\mathcal{F}_{\text{coh},\alpha\beta}`
and :math:`\mathcal{F}_{\text{inc},\alpha}`) atomic based properties.

.. _water-pdf-unweighted:

.. figure:: ./Pictures/water_pdf_unweighted.png
   :align: center
   :width: 11.748cm
   :height: 9.393cm

   The total and atomic intermolecular PDF of water, atomic DOS are
   **unweighted** so that the weighted sum of atomic DOS equals to the total.

The weighted and unweighted options are more useful for different cases, for example,
it might be more useful to use the weighted terms for the DOS calculations (:numref:`water-dos-weighted`)
while the unweighted terms might be more useful of the PDF calculations (:numref:`water-pdf-unweighted`).

Rescaled Weights
^^^^^^^^^^^^^^^^

MDANSE weights are rescaled so that weights for DISF calculation using the ``b_incoherent2`` will be

.. math::
   :label: ws9

   W_{\alpha} = \frac{c_{\alpha} b_{\mathrm{inc},\alpha}^2}{\sum_{\beta} c_{\beta} b_{\mathrm{inc},\beta}^2}

while the DCSF calculation using the ``b_coherent`` will be

.. math::
   :label: ws10

   W_{\alpha\beta} = \frac{c_{\alpha}c_{\beta} b_{\mathrm{coh},\alpha}b_{\mathrm{coh},\beta}}{\sum_{\gamma\delta} c_{\gamma}c_{\delta}  b_{\mathrm{coh},\gamma}b_{\mathrm{coh},\delta}}.

where :math:`b_{\mathrm{inc},\alpha}^2` is the square of the incoherent
scattering length of the atom type :math:`\alpha`.
:math:`b_{\mathrm{coh},\alpha}` and :math:`b_{\mathrm{coh},\beta}` are
the coherent scattering lengths of the atoms of types :math:`\alpha` and :math:`\beta`.
By using these rescaled weights the total incoherent and coherent intermediate
scattering functions becomes

.. math::
   :label: ws11

   \mathcal{F}_{\text{inc}}{(\mathbf{q},t ) = \frac{1}{\sum_{\alpha} c_{\alpha}   b_{\mathrm{inc},\alpha}^2 } \frac{1}{N}}{\sum\limits_{i} b_{\mathrm{inc},\alpha}^2 \left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{i}\left( t \right)} \right\rbrack} \right\rangle}


.. math::
   :label: ws12

   \mathcal{F}_{\text{coh}}{(\mathbf{q},t) = \frac{1}{\sum_{\alpha\beta} c_{\alpha} c_{\beta}  b_{\mathrm{coh},\alpha}b_{\mathrm{coh},\beta}} \frac{1}{N}}{{\sum\limits_{ij} b_{\mathrm{coh},i}b_{\mathrm{coh},j} \left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle}},


where :math:`b_{\mathrm{inc},i}^2` is the square of the incoherent
scattering length of atom :math:`i`. :math:`b_{\mathrm{coh},i}` and
:math:`b_{\mathrm{coh},j}` are the coherent scattering lengths of
atoms :math:`i` and :math:`j`. Notice that by using these weights the proper
total intermediate scattering functions will get rescaled by a number. Notice that by using this
weight scheme the total DISF has the property that

.. math::
   :label: ws13

   \mathcal{F}_{\text{inc}}(\mathbf{q},t=0) = 1.


However the total intermediate scattering function
(sum of the incoherent and coherent parts) will not be equal (or equal to the
sum by some scaling factor) to the to the sum of intermediate scattering function
from the DISF and DCSF calculations using this weight scheme since they
are not scaled in the same way.
