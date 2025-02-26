Weighting Scheme
================

In MDANSE, all properties that are calculated are split by atom-type
and the total results is a weighted sum of these atom-type based
properties. For example the coherent and incoherent are

.. math::
   :label: ws1

   F_{\alpha\beta}^{\text{coh}}{(\mathbf{q},t) = \frac{W_{\alpha\beta}}{\sqrt{N_{\alpha}N_{\beta}}}}{\sum\limits_{i}^{N_{\alpha}}{\sum\limits_{j}^{N_{\beta}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{\alpha}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle}},

.. math::
   :label: ws2

   F_{\alpha}^{\text{inc}}{(\mathbf{q},t ) = \frac{W_{\alpha}}{N_{\alpha}}}{\sum\limits_{i}^{N_{\alpha}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{i}\left( t \right)} \right\rbrack} \right\rangle}

where :math:`\alpha` and :math:`\beta` are the atom-types.
:math:`W_{\alpha\beta}` and :math:`W_{\alpha}` are the weights of the
atom-type pairs :math:`\alpha\beta` and the atom type :math:`\alpha`.
:math:`N_{\alpha}` and :math:`N_{\beta}` are the number of atoms of
atom-types :math:`\alpha` and :math:`\beta`. Notice that MDANSE 2 scales
the results using the prefactors :math:`W_{\alpha\beta} / \sqrt{N_{\alpha}N_{\beta}}`
and :math:`W_{\alpha} / N_{\alpha}`. The total is now a sum of
these atomic terms

.. math::
   :label: ws3

    F^{\text{coh}}(\mathbf{q},t) = \sum_{\alpha}\sum_{\beta \geq \alpha} \text{F}_{\alpha\beta}^{\text{coh}}(\mathbf{q},t),

.. math::
   :label: ws4

    F^{\text{inc}}(\mathbf{q},t) = \sum_{\alpha} \text{F}_{\alpha}^{\text{inc}}(\mathbf{q},t).

Note that for summation involving two atom-types only the unique pairs
are summed up. This is because in MDANSE the off-diagonal weight
terms are symmetrised so that we assumed that
:math:`F_{\alpha\beta}^{\text{coh}} = F_{\beta\alpha}^{\text{coh}}`.

.. _water-dos-weighted:

.. figure:: ./Pictures/water_dos_weighted.png
   :align: center
   :width: 11.748cm
   :height: 9.393cm

   The total and atomic DOS of water, atomic DOS are **weighted** so that the
   sum of atomic DOS equals to the total.

In the old MDANSE these atom properties were scaled slightly differently

.. math::
   :label: ws5

   \mathcal{F}_{\alpha\beta}^{\text{coh}}{(\mathbf{q},t) = \frac{1}{\sqrt{N_{\alpha}N_{\beta}}}}{\sum\limits_{i}^{N_{\alpha}}{\sum\limits_{j}^{N_{\beta}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{\alpha}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{j}\left( t \right)} \right\rbrack} \right\rangle}},

.. math::
   :label: ws6

   \mathcal{F}_{\alpha}^{\text{inc}}{(\mathbf{q},t ) = \frac{1}{N_{\alpha}}}{\sum\limits_{i}^{N_{\alpha}}\left\langle {\exp\left\lbrack {{- i}\mathbf{q}\cdot\mathbf{r}_{i}\left( 0 \right)} \right\rbrack\exp\left\lbrack {i\mathbf{q}\cdot\mathbf{r}_{i}\left( t \right)} \right\rbrack} \right\rangle}

so the total will now be a weighted sum of these atomic terms

.. math::
   :label: ws7

    F^{\text{coh}}(\mathbf{q},t) = \sum_{\alpha}\sum_{\beta \geq \alpha} W_{\alpha\beta} \mathcal{F}_{\alpha\beta}^{\text{coh}}(\mathbf{q},t),

.. math::
   :label: ws8

    F^{\text{inc}}(\mathbf{q},t) = \sum_{\alpha} W_{\alpha} \mathcal{F}_{\alpha}^{\text{inc}}(\mathbf{q},t).

In the MDANSE_GUI you have the option to plot either weighted (e.g. :math:`F_{\alpha\beta}^{\text{coh}}`
and :math:`F_{\alpha}^{\text{inc}}`) or unweighted (e.g. :math:`\mathcal{F}_{\alpha\beta}^{\text{coh}}`
and :math:`\mathcal{F}_{\alpha}^{\text{inc}}`) atomic based properties.

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

MDANSE weights are rescaled so that weights for the :math:`F_{\alpha}^{\text{inc}}`
using the ``b_incoherent2`` will be

.. math::
   :label: ws9

   W_{\alpha} = \frac{N_{\alpha} b_{\mathrm{inc},\alpha}^2}{\sum_{\beta} N_{\beta} b_{\mathrm{inc},\beta}^2}

while :math:`F_{\alpha}^{\text{coh}}` using the ``b_coherent`` will be

.. math::
   :label: ws10

   W_{\alpha\beta} = \frac{N_{\alpha}N_{\beta} b_{\mathrm{coh},\alpha}b_{\mathrm{coh},\beta}}{\sum_{\gamma\delta} N_{\gamma}N_{\delta}  b_{\mathrm{coh},\gamma}b_{\mathrm{coh},\delta}}.
