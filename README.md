# Raccoon Scrubber
Raccoon Scrubber is a ligand library processing tool to generate 3D molecules ready to dock starting from any supported file format (OB3)

## Dependencies
The code depends on OpenBabel v.3.x or newer.

## Usage
```
usage: scrubber.py [-h] --infile INPUT_FILE[.EXT] [--outfile OUTPUT_FILE[.EXT]] [--informat [mol2|sdf|smi|pdb|...]] [--outformat [mol2|sdf|pdb|...]] [--usemolname]
                   [--usemolnamesafe] [--usefieldname FIELD] [--single STRUCTURE_NUMBER] [--begin FIRST_STRUCTURE_NUMBER] [--end LAST_STRUCTURE_NUMBER] [--split]
                   [--byname MOL_NAME] [--pH [7.4]] [--nopH] [--noflipamide] [--nostripsalts] [--nocheckhydro] [--noprocess] [--enumchiral [undefined|protomer|all|off]]
                   [--maxenumchiral MAX_ENANTIOMERS] [--exclude SMARTS] [--excludefromfile FILENAME] [--sdsteps [ 300]] [--sdconv [1e-05]] [--cgsteps [300]]
                   [--cgconv [1e-06]] [--forcefield [mmff94s]] [--sdsteps_extra [1000]] [--sdconv_extra [1e-05]] [--cgsteps_extra [1000]] [--cgconv_extra [1e-06]]
                   [--forcefield_extra [uff]] [--rotamer_conf [10]] [--rotamer_steps [5]] [--nomini] [--noextra] [--norotamer] [--chargemodel [gasteiger]] [--strict]
                   [--multicore [8]] [--nice NICE] [--log LOGFILENAME] [--verbose] [--help_advanced]
```
## TODO
-[ ] disable rotameric search by default
-[ ] find optimal SD/CG default parameters
-[ ] test and activate automatic heuristic for minimization parameters
