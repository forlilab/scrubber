.. molscrub documentation master file, created by
   sphinx-quickstart on Mon May 11 10:46:13 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MolScrub documentation
======================

Process large numbers of small molecules for docking with AutoDock. May
be useful for structure-based modeling in general.

What happens:

-  generate 3D coordinates using RDKit's ETKDGv3 and UFF minimization
-  enumerate tautomers (aiming at low energy states only)
-  enumerate pH corrections
-  convert boats to chairs (6-member rings) and enumerate both chair
   states
-  enumerate chiral centers (not implemented right now)

Installation
============

Direct installation:

.. code:: sh

   conda activate <desired-environment>    # if you are using conda environments

   pip install git+https://github.com/forlilab/molscrub.git

For developers:

.. code:: sh

   conda activate <desired-environment>    # if you are using conda environments

   git clone git@github.com:forlilab/molscrub.git
   cd molscrub
   pip install -e .

Depends on the RDKit, which can be installed from conda-forge in the
desired environment:

.. code:: sh

   conda activate <desired-environment>
   conda install rdkit -c conda-forge

Python scripting
----------------

.. code:: python

   from rdkit import Chem
   from molscrub import Scrub

   scrub = Scrub(
       ph_low=7.4,
       ph_high=7.4,
   )

   mol = Chem.MolFromSmiles("Clc1c(OCCC3)c3ccc1C(=O)Nc2nc[nH]c2")

   # each state (e.g. tautomer) an rdkit mol and may have multiple conformers
   for mol_state in scrub(mol):
       print(Chem.MolToSmiles(mol_state), "nr conformers: %d" % mol_state.GetNumConformers())

Command line tool examples
--------------------------

.. code:: sh

   scrub.py "c1cc[nH]c(=O)c1" -o scrubbed.sdf --ph 5 --skip_gen3d
   scrub.py input_mols.sdf -o scrubbed.sdf
   scrub.py input_mols.smi -o scrubbed.sdf

Other options described in the help message:

.. code:: sh

   scrub.py -h

Where "input_mols.smi" can look like this:

::

   CC(=O)O aceticacid
   CN(C)C trimethylamine 
   Clc1cc(O)ccc1C(=O)Nc2nc[nH]c2 hello_mol
   c1cccc1 rdkit_will_cry
   CCC good4bbq
   CCO alsogood4bbq
   c1cccnc1CC(=O)C a_ketone

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Contents

   Command Line Options <cli>
   Using the API <api>
   API Modules <modules>


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: How it works

    Acid/base conjugation <acidbase>
