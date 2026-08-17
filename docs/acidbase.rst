Acid/base conjugation
=====================


Setting ph_low and ph_high
__________________________

We have mostly used molscrub with pH=7.4 (the default). However. it possible to specify a range with
CLI options ``--ph_low 8 --ph_high 10``, or directly in Python:

.. code-block:: python

    scrub = Scrub(ph_low=8, ph_high=10)

Here, 8 and 10 are just examples.

MolScrub predicts the pKa with either the SMARTS based rules or the Extra Trees model. The default is ``rules``
and can be changed with ``--pka_model``. Then, the selected protonation states depend on how the predicted pKa
compares with the pH limits:

.. list-table:: Selection of acid/base states
   :widths: 25 15
   :header-rows: 1

   * - Condition
     - Returned states
   * - pKa < ph_low
     - base
   * - pKa > ph_high
     - acid
   * - ph_low <= pKa <= ph_high
     - both acid and base
 
Each acid/base site is evaluated independently, and the final number of molecular states is the combination of selected states for each site.
If three sites have their pKa within ``ph_low`` and ``ph_high``, eight protonation states are returned (:math:`2^3`).


When ``ph_low == ph_high`` it is rare to select two states for a single site. That happens only if
the predicted pKa matches exactly the pH limits. For most molecules, only one state is returned
by the acid/base conjugation module when the pH range is zero.
