import multiprocessing
import argparse

# provide method to convert strings into boolean
from .transform.isomer import MoleculeIsomers
from .storage import MoleculeProvider, MoleculeStorage, MoleculeIssueStorage
from .geom.geometry import ParallelGeometryGenerator
from .core import ScrubberCore

def strtobool(val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))


# TODO source:
# https://www.golinuxcloud.com/python-argparse/

__doc__ = """
This file contains all the information related to the CLI.
Command-line options are defined here, together with their description. All options are initialized to the default values of ther corresponding classes for each function.
"""


"""These are the tags (prefixes) that get automatically prepended to all options
pertaining a given function class.
"""
tags = {
    "input": "in_",
    "output": "out_",
    "isomers": "iso_",
    "geometry": "geom_",
    "errors": "err_",
    "general": "",
}
tags_reverse = {v[:-1]: k for k, v in tags.items() if v}


description = """
DESCRIPTION
This is the Raccoon Scrubber. This is the new version based on RDKit
"""

usage = """R|

In order to get something done with the Scrubber this is the command:

    $ %s [options] filename

       This is indented
           and this indented more

""" % (
    "scrubber.py"
)

epilog = """ (C) 2022 ForliLab, Scripps Research """


advanced_help = """

This is the advanced help

                      _(\-/)_
                     {(#b^d#)}
                     `-.(Y).-`

This is where the full pipeline is described (scrubber.core.process())
"""

from pprint import pprint as pp

isomer_default = MoleculeIsomers.get_defaults()
molprovider_default = MoleculeProvider.get_defaults()
molstorage_default = MoleculeStorage.get_defaults()
molerror_default = MoleculeIssueStorage.get_defaults()
geom_default = ParallelGeometryGenerator.get_defaults()
general_default = ScrubberCore.get_defaults()["general"]["values"]

cli_options = {
    "input": {
        "description": "options for input definition and parsing.",
        "values": {
            "--in_fname": {
                "help": """ [ REQUIRED ] input file to process (SMI, SDF). The file type is
                guessed from the extension, unless the --in_ftype option is used. """,
                "action": "store",
                "metavar": "INPUT_FILE[.EXT]",
                # "required": True,
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--in_ftype": {
                "help": """specify the format of the input file, overriding
                the extension (i.e., parsing SMILES files with '.txt' extension);
                the format has to be one of the supported types (SMI or SDF); NOTE:
                if a molecule does not have a name defined, the \"MOL\" name will
                be set by default [ default: guessed from --in_fname ]""",
                "action": "store",
                "metavar": "EXT",
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--in_sanitize": {
                "help": """perform automatic sanitization of the input
                using RDKit [ default: %s ]"""
                % (str(molprovider_default["sanitize"])),
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(molprovider_default["sanitize"]),
                "default": argparse.SUPPRESS,
            },
            "--in_removeHs": {
                "help": """remove hydrogens from input structure, if present
                [ default: %s ]"""
                % str(molprovider_default["removeHs"]),
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(molprovider_default["removeHs"]),
                "default": argparse.SUPPRESS,
            },
            "--in_strictParsing": {
                "help": """perform strict parsing with RDKit; molecules will be skipped
                if any inaccuracies are found [CLARIFY THIS!]  [ default: %s ]"""
                % (
                    str(
                        molprovider_default["strictParsing"],
                    )
                ),
                "action": "store",
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(molprovider_default["strictParsing"]),
                # "default": molprovider_default["strictParsing"],
                "default": argparse.SUPPRESS,
            },
            "--in_name_property": {
                "help": """[SDF ONLY] use the specified property field to name the
                molecule. If the --out_split option is used, this will be the filename.""",
                "action": "store",
                "metavar": '"property name"',
                "required": False,
                "type": str,
                # "default": None,
                "default": argparse.SUPPRESS,
            },
            "--in_safeparsing": {
                "help": """parse input files in safe mode, instead of using RDKit
                native parsers; enable saving problematic raw input text into log files;
                it could be slower than using unsafe parsing [ default: %s ]"""
                % molprovider_default["safeparsing"],
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(molprovider_default["safeparsing"]),
                # "default": molprovider_default["safeparsing"],
                "default": argparse.SUPPRESS,
            },
            "--in_discarded_datafile": {
                "help": """specify a file where problematic data from the input file
                will be saved if parsing errors occurr. The file will be in text format,
                to preserve the original input data [ default: %s ]"""
                % molprovider_default["discarded_datafile"],
                "action": "store",
                "metavar": "PROBLEMATIC_INPUT.EXT",
                "required": False,
                "type": type(molprovider_default["discarded_datafile"]),
                # "default": None,
                "default": argparse.SUPPRESS,
            },
            "--in_start_count": {
                "help": """start processing the input file from the specified n-th
                molecule (included); molecules before the specified count will be
                ignored [ default: %s ]"""
                % molprovider_default["start_count"],
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": type(molprovider_default["start_count"]),
                # "default": None,
                "default": argparse.SUPPRESS,
            },
            "--in_end_count": {
                "help": """stop processing the file at the specified n-th molecule
                (included); molecules after the specified count will be ignored
                [ default: %s ]"""
                % molprovider_default["end_count"],
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": type(molprovider_default["end_count"]),
                # "default": None,
                "default": argparse.SUPPRESS,
            },
        },
    },
    "output": {
        "description": "options for output definition and saving.",
        "values": {
            "--out_fname": {
                "help": """[ REQUIRED ] specify the output filename; by default the file
                extension will be used to set the file format, unless the option
                \'--out_ftype\' is used; if a fullpath with explicit directories
                is specified, the directories will be created automatically; in
                \"split\" mode, the file name will be prepended to whatever value
                is specified in the \"--out_naming\" option""",
                "action": "store",
                "metavar": "OUTPUT_FNAME[.EXT]",
                # "required": True,
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--out_ftype": {
                "help": """specify the format of the output file, overriding the
                extension (i.e., writing SMILES files with '.txt' extension); the
                format has to be one of the supported types (SMI or SDF); by default
                the format type is guessed from the file specified in --out_fname""",
                "action": "store",
                "metavar": "EXT",
                "choices": ("smi", "sdf"),
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--out_format_opts": {
                "help": """format-specific options: SMI and SDF (some stuff from
                RDKIT?) -> SMI:titleLine option""",
                "action": "store",
                "metavar": "SOMETHING HERE",
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--out_mode": {
                "help": """define how the output data is saved; by default the
                \'single\' mode is used, writing all the output in a single file;
                if the \'split\' mode is used, then each molecule will be saved in
                individual files; the output filenames can be controlled with the
                \"--out_naming\" option [ default: %s ] """
                % molstorage_default["mode"],
                "action": "store",
                # "metavar": "single|split",
                "choices": ("single", "split"),
                "type": type(molstorage_default["mode"]),
                # "default": molstorage_default["mode"],
                "default": argparse.SUPPRESS,
            },
            "--out_naming": {
                "help": """[SPLIT MODE ONLY] when saving a file for each molecule,
                define which naming scheme to use, by default 'auto' mode is used;
                alternatively the 'name' mode can be useed XXXX EXPAND """,
                "action": "store",
                # "metavar": "auto|name",
                "choices": ("auto", "name"),
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--out_max_lig_per_dir": {
                "help": """[SPLIT MODE ONLY] create progressively numbered
                sub-directories containing no more than the specified number
                of ligands [ default: disabled ]""",
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": int,
                "default": argparse.SUPPRESS,
            },
            "--out_disable_name_sanitize": {
                "help": """[SPLIT MODE ONLY] by default the molecule name used for
                the output file name is sanitized by removing spaces, parentheses
                ( "{}", "[]" ) and other unsafe characters; this option disables it
                [ default: %s ]"""
                % str(molstorage_default["disable_name_sanitize"]),
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": bool,
                # "default": False,
                "default": argparse.SUPPRESS,
            },
            "--out_disable_preserve_properties": {
                "help": 'by default, if the input molecules have properties associated (i.e. SDF extra fields, "PUBCHEM_COMPOUND_CID"), they will be preserved in the output file; this option disables it so no properties will be saved in the output [SDF ONLY]',
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": bool,
                # "default": False,
                "default": argparse.SUPPRESS,
            },
        },
    },
    "isomers": {
        "description": """Options for all isomeric transformations (tautomers,
        protomers, stereoisomers)""",
        "values": {
            #
            # stereoisomers
            #
            "--iso_stereo_enum": {
                "help": """enable stereoisomer enumeration; it is possible to process
                only stereo centers with unspecified chirality (\"undefined\"), or all
                stereocenters (\"all\") (including those in for which explicit chirality
                is defined) [ default: disabled ]""",
                "action": "store",
                # "metavar": "undefined|all",
                "choices": ("undefined", "all"),
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--iso_stereo_max_results": {
                "help": """maximum number of stereoisomers to enumerate; if this number is
                exceeded, the enumeration is stopped [ default : %d ]"""
                % isomer_default["stereo_max_results"],
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(isomer_default["stereo_max_results"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_stereo_gen3d": {
                "help": """enable the generation of 3D coordinates for the newly
                enumerated steroisomers (TODO: MAYBE REMOVE? WHEN NEEDED?)[ default: %s ]"""
                % str(isomer_default["stereo_gen3d"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["stereo_gen3d"]),
                # "default": isomer_default["stereo_gen3d"],
                "default": argparse.SUPPRESS,
            },
            #
            # protomers
            #
            "--iso_proto_enum": {
                "help": """enable/disable protomer enumeration; the pH value/ranges
                can be specified with the \"--iso_proto_pH\" option [ default: %s ]"""
                % str(isomer_default["proto_enum"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["proto_enum"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_proto_pH": {
                "help": """[TODO FIX BY SPECIFYING A RANGE] specify the pH for the protomer enumeration; if a
                single value is specified, then all reactions above that pH
                will be performed (e.g. '--iso_proto_pH 7.4'); if two
                comma-separated values are specified then all protomer
                transformation within that range will be peformed (e.g.:
                '--iso_proto_pH 6.4,8.4) [ default: %2.1f ]"""
                % isomer_default["proto_pH"],
                "action": "store",
                "metavar": "pH_value|pH_min,pH_max",
                "required": False,
                # "type": type(isomer_default["proto_pH"]),
                # this should be kept as a str type because it can be
                # either a float or a tuple of floats
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--iso_proto_max_results": {
                "help": """maximum number of protomers to enumerate; if this
                number if exceeded, the generation is stopped and only the
                [max_results] results enumerated will be processed further [ default: %d ]"""
                % isomer_default["proto_max_net_charge"],
                "action": "store",
                "metavar": "INT",  # % isomer_default["proto_max_results"],
                "required": False,
                "type": type(isomer_default["proto_max_results"]),
                # "default": isomer_default["proto_max_results"],
                "default": argparse.SUPPRESS,
            },
            "--iso_proto_keep_all": {
                "help": """by default, exaustive protomer generation is performed
                and only final products are reported; if this flag is used, all
                intermediate protomers are kept, too [ default: %s ]"""
                % isomer_default["proto_keep_all"],
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["proto_keep_all"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_proto_max_net_charge": {
                "help": """molecules which absolute formal charges exceed this
                value are discarded [ default: %d ]"""
                % isomer_default["proto_max_net_charge"],
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(isomer_default["proto_max_net_charge"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_proto_neutralize_only": {
                "help": """generate neutral form of input molecules. WARNING:
                when used, all other protomer transformations will be disabled
                [default: %s]"""
                % isomer_default["proto_neutralize_only"],
                "action": "store_true",
                "required": False,
                # "default": isomer_default["proto_neutralize"],
                "default": argparse.SUPPRESS,
            },
            #
            # tautomers
            #
            "--iso_tauto_enum": {
                "help": "enable tautomer enumeration [ default: %s ]"
                % isomer_default["tauto_enum"],
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["tauto_enum"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_tauto_max_results": {
                "help": """max number of tautomers to enumerate
                    [ default: %d ]"""
                % isomer_default["tauto_max_results"],
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(isomer_default["tauto_max_results"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_tauto_protect_aromatic": {
                "help": """when generating tautomers, discard those reducing
                the number of aromatic atoms in the system; by setting this
                option to 'False', less stable tautomers with fewer aromatic
                atoms will be kept [ default: %s ]"""
                % str(isomer_default["tauto_protect_aromatic"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["tauto_protect_aromatic"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_tauto_protect_amide": {
                "help": """when generating tautomers, discard those reducing
                the number of amide groups in the system; by setting this
                option to 'False', less stable tautomers with fewer amide
                groups atoms will be kept [ default: %s ]"""
                % str(isomer_default["tauto_protect_amide"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["tauto_protect_amide"]),
                "default": argparse.SUPPRESS,
            },
            #
            # generic
            #
            "--iso_add_hydrogens": {
                "help": """when enumerating isomers, add hydrogens to the
                    final products [ default: %s ]"""
                % str(isomer_default["add_hydrogens"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["add_hydrogens"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_max_cycles": {
                "help": """when generating both protomers and tautmers, this
                option defines how many times the process (protomer+tautomers)
                is repeated to make sure all products are generated [ default: %d ]"""
                % isomer_default["max_cycles"],
                "action": "store",
                "required": False,
                "metavar": "INT",
                "type": type(isomer_default["max_cycles"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_verbose": {
                "help": """enable verbose mode for all isomeric transformations
                to track the procress -- WARNING: this can generate a lot of
                data when processing many molecules [ default: %s ]"""
                % str(isomer_default["verbose"]),
                "action": "store_true",
                # "metavar": "TRUE|FALSE",
                # "type": type(isomer_default["verbose"]),
                "default": argparse.SUPPRESS,
            },
            "--iso_ph_datafile": {
                "help": """specify a custom file for the pH model transformations
                (see advanced help)""",
                "action": "store",
                "metavar": "FILENAME",
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--iso_tauto_datafile": {
                "help": """specify a custom file for the tautomer model transformations
                (see advanced help)""",
                "action": "store",
                "metavar": "FILENAME",
                "type": str,
                "default": argparse.SUPPRESS,
            },
            "--iso_suppress_rdkit_warnings": {
                "help": """disable warning messages from RDKit when sanitizing molecules [ default: %s ]"""
                % str(isomer_default["suppress_rdkit_warnings"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(isomer_default["add_hydrogens"]),
                "default": argparse.SUPPRESS,
            },
        },
    },
    "geometry": {
        "description": "Options the generation of accurate 3D coordinates",
        "values": {
            "--geom_add_h": {
                "help": """add hydrogens to the newly generated 3D structures [ default: %s ]"""
                % str(geom_default["add_h"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(geom_default["add_h"]),
                # "default": geom_default["add_h"],
                "default": argparse.SUPPRESS,
            },
            "--geom_force_trans_amide": {
                "help": """guarantee that acyclic secondary amide are in a
                trans conformation [ default: %s ]"""
                % str(geom_default["force_trans_amide"]),
                "action": "store",
                "required": False,
                # "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(geom_default["force_trans_amide"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_force_field": {
                "help": """set the force field to be used for the 3D coordinates
                optimization [ default: %s ]"""
                % geom_default["force_field"],
                "action": "store",
                # "metavar": "uff|mmff94|mmff94s",
                "choices": ["uff", "mmff94", "mmff94s"],
                "required": False,
                "type": type(geom_default["force_field"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_max_iterations": {
                "help": """maximum number of iterations for to perform at each
                minimization cycle [ default: %d ]"""
                % geom_default["max_iterations"],
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(geom_default["max_iterations"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_auto_iter_cycles": {
                "help": """set the maximum number of minimization cycles (or attempts)
                performed on molecules to achieve convergence; the number of iterations
                at each cycle is defined by the option \"--gen_max_iterations\"
                [ defaults: %d ]"""
                % geom_default["auto_iter_cycles"],
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(geom_default["auto_iter_cycles"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_gen3d": {
                "help": """generate new 3D coordinates prior to optimizing them using
                the ETKDG implemented in RDKit; if set to False, this option can be
                used to further minimize molecules with valid 3D coordinates [ default: %s ]"""
                % geom_default["gen3d"],
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(geom_default["gen3d"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_gen3d_max_attempts": {
                "help": """define how many times initial generation of 3D coordinates is
                performed, in case of failure [ default: %d ]"""
                % geom_default["gen3d_max_attempts"],
                "action": "store",
                "metavar": "INT",  # geom_default["gen3d_max_attempts"],
                "required": False,
                "type": type(geom_default["gen3d_max_attempts"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_fix_ring_corners": {
                "help": """process 6-membered rings to prevent boat conformatons and
                identify the optimal chair conformations [ default: %s ]"""
                % str(geom_default["fix_ring_corners"]),
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(geom_default["fix_ring_corners"]),
                "default": argparse.SUPPRESS,
            },
            "--geom_preserve_mol_properties": {
                "help": """[SDF ONLY] preserve properties in the input
                    molecule [ default: %s ]"""
                % str(geom_default["preserve_mol_properties"]),
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                "default": argparse.SUPPRESS,
            },
            "--geom_strict": {
                "help": """reject molecules for which geometry optimization did not
                converge (to rescue them, increase --geom_auto_iter_cycles and/or
                --geom_max_iterations) [ default: %s ]"""
                % str(geom_default["strict"]),
                "action": "store",
                "required": False,
                "metavar": "TRUE|FALSE",
                "choices": [True, False],
                "type": lambda x: bool(strtobool(x)),
                # "type": type(geom_default["strict"]),
                "default": argparse.SUPPRESS,
            },
        },
    },
    "general": {
        "description": "General options for the calculation setup",
        "values": {
            "--max_proc": {
                "help": """maximum number of processors/cores to use [ default: %d on this machine ]"""
                % multiprocessing.cpu_count(),
                "action": "store",
                "metavar": "PROC_NUM",
                "required": False,
                "type": type(general_default["max_proc"]),
                "default": argparse.SUPPRESS,
            },
            "--nice_level": {
                "help": """set the nice level (priority) of the processing; values range
                from 0 (standard priority) to 40 (lowest priority) [ default: 0 ]""",
                "action": "store",
                "metavar": "NICE_LEVEL",
                "required": False,
                "type": int,
                # "type": type(general_default["nice_level"]),
                # "default": general_default["nice_level"],
                "default": argparse.SUPPRESS,
            },
        },
    },
    "errors": {
        "description": "Options to save problematic data",
        "values": {
            "--err_log_basename": {
                "help": """the file basename that will be used to save
                problematic data as following: BASENAME_input.txt (raw data
                text from problematic input); BASENAME_processing.smi
                (problematic molecules processed without coordinates);
                BASENAME_processing.sdf (problematic molecules processed with
                coordinates)  where to save problematic input structures; the
                file will be the same type as the input file [ default:
                \"scrub_errors\" ]""",
                # % str(molerror_default["log_basename"]),
                "action": "store",
                "metavar": "FILENAME",
                "required": False,
                "type": str,
                "default": argparse.SUPPRESS,
            },
            # "--err_log_from_input": {
            #     "help": """filename where to save problematic input structures; the file
            #     will be the same type as the input file [ default: %s ]"""
            #     % str(molerror_default["log_from_input"]),
            #     "action": "store",
            #     "metavar": "FILENAME",
            #     "required": False,
            #     "type": str,
            #     "default": argparse.SUPPRESS,
            # },
            # "--err_log_from_process": {
            #     "help": """filename where to save problematic structures encountered during
            #     the processing (e.g.: geometry minimization); the file format (SDF, SMI)
            #     will be automatically selected from the dimensionality of the first molecule
            #     (TO BE CLAFIRIED [ default: %s ]"""
            #     % str(molerror_default["log_from_process"]),
            #     "action": "store",
            #     "metavar": "FILENAME",
            #     "required": False,
            #     "type": str,
            #     "default": argparse.SUPPRESS,
            # },
            # "--process_err_ftype": {
            #     "help": "file format of the processing errors file; allowed formats are SMI and SDF. ",
            #     "action": "store",
            #     "metavar": "smi|sdf",
            #     "required": False,
            #     "type": type(molerror_default["process_err_ftype"]),
            #     "default": molerror_default["process_err_ftype"],
            # },
        },
    },
}


"""these are extra options providing convenient CLI functionalities that are
absent in the ScrubberCore"""
extra_options = {
    "--config_file": {
        "help": (
            "specify a JSON config file to read to configure the processing;"
            " if used, additional command-line options will superseed the config file values; "
            "a template can be generated using the '--generate_template_config' option"
        ),
        "action": "store",
        "metavar": "CONFIG_FILE.JSON",
        "required": False,
        "type": str,
        "default": argparse.SUPPRESS,
    },
    "--save_config_template": {
        "help": (
            "save a template file containing all configuration options with their "
            "default values, to be customized and used with the '--config_file' option. "
            "If the file exists, the program will cowardly refuse to overwrite it"
        ),
        "action": "store",
        "metavar": "CONFIG_FILE.JSON",
        "required": False,
        "type": str,
        "default": argparse.SUPPRESS,
    },
    # check if we want to keep this
    "--help_advanced": {
        "help": "print the advanced help",
        "action": "store_const",
        "required": None,
        "const": True,
        "default": argparse.SUPPRESS,
    },
    "--quiet": {
        "help": "suppress printing of the progress and final summary",
        "action": "store_const",
        "required": None,
        "const": True,
        # "default": argparse.SUPPRESS,
        "default": False,
    },
    "--overwrite": {
        "help": "Overwrite output file, if exists",
        "action": "store_const",
        "required": None,
        "const": True,
        # "default": argparse.SUPPRESS,
        "default": False,
    },
}
