
import sys
import os
import scrub as s

"""
Simple script to test/tune the performance of the 
Scrub pipeline and find optimal parameters
"""

# read input file (SMILES)
fname = sys.argv[1]
name, ext = os.path.splitext(fname)
ext = ext[1:]
# read the first line of the SMILES
with open(fname, 'r') as fp:
    data = fp.readlines()[1]
data = "".join(data)
# set output file name
outname = sys.argv[2]

# define options to test
opts = {'inFormat':ext.lower(),
        'sdsteps':1000,
        'cgsteps':1000,
        'sdconv':'1e-1',
        'cgconv':'1e-2',
        'rotamer_conf':0,
        'sdsteps_extra':500,
        'cgsteps_extra':500,
        'sdconv_extra':'1e-2',
        'cgconv_extra':'1e-2',
        }

# create an instance of Scrub
scrub = s.Scrub(**opts)
# run the full pipeline
out = scrub.process(data, in_type = ext, out_type='mol2')
# write the output
with open(outname, 'w') as fp:
    fp.write(out)
