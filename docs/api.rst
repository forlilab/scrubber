Using the Scrub API
=====================

The python API revolves around the ``Scrub`` object which instantiates a 
function that is used to scrub molecules. The ``Scrub`` class takes more 
or less the same as the :doc:`./cli`. See the signature of  
:class:`molscrub.core.Scrub` class in the :doc:`API Docs <molscrub>`

To use the class:

.. code:: python

   from rdkit import Chem
   from molscrub import Scrub

   scrub = Scrub(
       ph_low=7.4,
       ph_high=7.4,
       # ... other options... 
   )

   mol = Chem.MolFromSmiles("Clc1c(OCCC3)c3ccc1C(=O)Nc2nc[nH]c2")

   scrubbed_mols = scrub(mol)

The instatiated ``scrub`` function takes an RDKit mol as input and 
produces an list of RDKit mols as output (with embedded 3D coords). 
Then you can handle these or perform any operations that you would 
on RDKit mols. 