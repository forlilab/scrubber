import multiprocessing

from .transform.isomer import MoleculeIsomers
from .storage import MoleculeProvider, MoleculeStorage
from .geom.geometry import ParallelGeometryGenerator
from .core import ScrubberCore

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

epilog = """ (C) 2022 ForliLab Industries, Scripps Research """


advanced_help = """

This is the advanced help

                      _(\-/)_
                     {(#b^d#)}
                     `-.(Y).-`

This is where the full pipeline is described (scrubber.core.process())
"""

# TODO make sure they all have reasonable or accurate defaults

isomer_default = MoleculeIsomers.get_defaults()
molprovider_default = MoleculeProvider.get_defaults()
molstorage_default = MoleculeStorage.get_defaults()
geom_default = ParallelGeometryGenerator.get_defaults()
general_default = ScrubberCore.get_defaults()["general"]["values"]

cli_options = {
    "input": {
        "description": "options for input definition and parsing.",
        "values": {
            "--in_fname": {
                "help": "input file to process (SMI, SDF). The file type is guessed from the extension, unless the --in_ftype option is used.",
                "action": "store",
                "metavar": "INPUT_FILE[.EXT]",
                "required": True,
                "type": str,
            },
            "--in_ftype": {
                "help": "specify the format of the input file, overriding the extension (i.e., parsing SMILES files with '.txt' extension); the format has to be one of the supported types (SMI or SDF); NOTE: if a molecule does not have a name defined, the \"MOL\" name will be set by default",
                "action": "store",
                "metavar": "EXT",
                "required": False,
                "type": str,
            },
            "--in_sanitize": {
                "help": "perform automatic sanitization of the input using RDKit",
                "action": "store",
                "required": False,
                "metavar": "[%s]"
                % (
                    str(
                        molprovider_default["sanitize"],
                    )
                ),
                "type": type(molprovider_default["sanitize"]),
                "default": molprovider_default["sanitize"],
            },
            "--in_removeHs": {
                "help": "remove hydrogens from input structure, if present",
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "metavar": "[%s]"
                % (
                    str(
                        molprovider_default["removeHs"],
                    )
                ),
                "type": type(molprovider_default["removeHs"]),
                "default": molprovider_default["removeHs"],
            },
            "--in_strictParsing": {
                "help": "perform strict parsing with RDKit; molecules will be skipped if any inaccuracies are found [CLARIFY THIS!]",
                "action": "store",
                "metavar": "[%s]"
                % (
                    str(
                        molprovider_default["strictParsing"],
                    )
                ),
                "type": type(molprovider_default["strictParsing"]),
                "default": molprovider_default["strictParsing"],
            },
            "--in_name_property": {
                "help": "[SDF only] use the specified property field to name the molecule. If the --out_split option is used, this will be the filename.",
                "action": "store",
                "metavar": '"property name"',
                "required": False,
                "type": str,
                "default": None,
            },
            "--in_safeparsing": {
                "help": "parse input files in safe mode, instead of using RDKit native parsers; enable saving problematic raw input text into log files; it can be slower than using unsafe parsing",
                "action": "store",
                "metavar": "%s" % (str()),
                "required": False,
                "metavar": "[%s]"
                % (
                    str(
                        molprovider_default["safeparsing"],
                    )
                ),
                "type": type(molprovider_default["safeparsing"]),
                "default": molprovider_default["safeparsing"],
            },
            "--in_discarded_datafile": {
                "help": "save raw input text of problematic molecules in the specified datafile for further inspection",
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": str,
                "default": None,
            },
            "--in_start_count": {
                "help": "start processing the file from the specified n-th molecule (included)",
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": int,
                "default": None,
            },
            "--in_end_count": {
                "help": "stop processing the file at the specified n-th molecule (included)",
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": int,
                "default": None,
            },
        },
    },
    "output": {
        "description": "options for output definition and saving.",
        "values": {
            "--out_mode": {
                "help": """define how the output data is saved; by default the
            \'single\' mode is used, writing all the output in a single file;
            if the \'split\' mode is used, then each molecule will be saved in
            individual files; the output filenames can be controlled with the XXXX option """,
                "action": "store",
                "metavar": "[single] | split",
                "type": type(molstorage_default["mode"]),
                "default": molstorage_default["mode"],
            },
            "--out_fname": {
                "help": """specify the output filename; by default the file
            extension will be used to set the file format, unless the option
            \'--out_ftype\' is used; if a fullpath with explicit directories
            will be specified, the directories will be created automatically,
            if necessary """,
                "action": "store",
                "metavar": "OUTPUT_FNAME[.EXT]",
                "required": True,
                "type": str,
            },
            "--out_ftype": {
                "help": "specify the format of the output file, overriding the extension (i.e., writing SMILES files with '.txt' extension); the format has to be one of the supported types (SMI or SDF).",
                "action": "store",
                "metavar": "EXT",
                "required": False,
                "type": str,
            },
            "--out_format_opts": {
                "help": "format-specific options: SMI and SDF (some stuff from RDKIT?) -> SMI:titleLine option",
                "action": "store",
                "metavar": "SOMETHING HERE",
                "required": False,
                "type": str,
            },
            "--out_naming": {
                "help": "[split mode only] when saving a file for each molecule, define which naming scheme to use, by default 'auto' mode is used, ; alternatively the 'name' mode can be useed",
                "action": "store",
                "metavar": "'auto'|'name'",
                "required": False,
                "type": str,
            },
            "--out_max_lig_per_dir": {
                "help": "[SPLIT MODE ONLY] create sub-directories containing no more than the specified number of ligands",
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": int,
            },
            "--out_disable_name_sanitize": {
                "help": '[SPLIT MODE ONLY] by default the molecule name used for the output file name is sanitized by removing spaces, parentheses ( "{}", "[]" ) and other unsafe characters; this option disables it',
                "action": "store",
                "metavar": "NUMBER",
                "required": False,
                "type": bool,
                "default": False,
            },
            "--out_disable_preserve_properties": {
                "help": 'by default, if the input molecules have properties associated (i.e. SDF extra fields, "PUBCHEM_COMPOUND_CID"), they will be preserved in the output file; this option disables it so no properties will be saved in the output [SDF only]',
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": bool,
                "default": False,
            },
        },
    },
    "filter_pre": {
        # TODO this needs to be updated once the classes are completed
        "description": "Options for pre-filter settings (i.e. before any processing occurs)",
        "values": {
            "--filtpre_mw_max": {
                "help": "max mw to accept",
                "action": "store",
                "metavar": "MAX_MW",
                "required": False,
                "type": float,
            },
            "--filtpre_mw_min": {
                "help": "min mw to accept",
                "action": "store",
                "metavar": "MIN_MW",
                "required": False,
                "type": float,
            },
            "--filtpre_num_at_min": {
                "help": "min num atoms to accept",
                "action": "store",
                "metavar": "MIN_NUM_ATOMS",
                "required": False,
                "type": float,
            },
            "--filtpre_num_at_max": {
                "help": "max num atoms to accept",
                "action": "store",
                "metavar": "MAX_NUM_ATOMS",
                "required": False,
                "type": float,
            },
            "--filtpre_smarts_wanted": {
                "help": "SMARTS pattern of wanted molecules (DEFINE A FILE HERE)",
                "action": "store",
                "metavar": "SMARTS_STRING",
                "required": False,
                "type": str,
            },
            "--filtpre_smarts_not_wanted": {
                "help": "SMARTS pattern of not wanted molecules (DEFINE A FILE HERE)",
                "action": "store",
                "metavar": "SMARTS_STRING",
                "required": False,
                "type": str,
            },
            "--filtpre_pains_family": {
                "help": "PAINS family of patterns to reject (allowed: a, b, all)",
                "action": "store",
                "metavar": "[a|b|c|all]",
                "required": False,
                "type": str,
            },
        },
    },
    "filter_post": {
        "description": "Options for post-filter settings (i.e. before any processing occurs)",
        "values": {
            "--filtpost_mw_max": {
                "help": "max mw to accept",
                "action": "store",
                "metavar": "MAX_MW",
                "required": False,
                "type": float,
            },
            "--filtpost_mw_min": {
                "help": "min mw to accept",
                "action": "store",
                "metavar": "MIN_MW",
                "required": False,
                "type": float,
            },
            "--filtpost_num_at_min": {
                "help": "min num atoms to accept",
                "action": "store",
                "metavar": "MIN_NUM_ATOMS",
                "required": False,
                "type": float,
            },
            "--filtpost_num_at_max": {
                "help": "max num atoms to accept",
                "action": "store",
                "metavar": "MAX_NUM_ATOMS",
                "required": False,
                "type": float,
            },
            "--filtpost_smarts_wanted": {
                "help": "SMARTS pattern of wanted molecules (DEFINE A FILE HERE)",
                "action": "store",
                "metavar": "SMARTS_STRING",
                "required": False,
                "type": str,
            },
            "--filtpost_smarts_not_wanted": {
                "help": "SMARTS pattern of not wanted molecules (DEFINE A FILE HERE)",
                "action": "store",
                "metavar": "SMARTS_STRING",
                "required": False,
                "type": str,
            },
            "--filtpost_pains_family": {
                "help": "PAINS family of patterns to reject (allowed: a, b, all)",
                "action": "store",
                "metavar": "[a|b|c|all]",
                "required": False,
                "type": str,
            },
        },
    },
    "reaction": {
        "description": "Options for controlling chemical reactions",
        "values": {
            "--react_mode": {
                "help": 'specify the mode of reaction to be exhaustive ("all") or individual enumeration of each reaction site ("single"))',
                "action": "store",
                "metavar": "all|single",
                "required": False,
                "type": str,
            },
            "--react_reactions_file": {
                "help": "specify a file containing SMIRKS reactions (one per line) to be used for the transformations",
                "action": "store",
                "metavar": "FILENAME",
                "required": False,
                "type": str,
            },
            "--react_verbose": {
                "help": "enable verbose mode for the reaction to track the procress -- WARNING: this can generate a lot of data when processing many molecules",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": bool,
            },
            "--react_keep_only_reacted": {
                "help": "keep only molecules that reacted; molecules that do not react are discarded",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": bool,
            },
        },
    },
    "isomers": {
        "description": "Options for all isomeric transformations (tautomers, protomers, stereoisomers)",
        "values": {
            #
            # stereoisomers
            #
            "--iso_stereo_enum": {
                "help": "enable stereoisomer enumeration; it is possible to process only stereo centers with unspecified chirality, or all stereocenters (including those in for which explicit chirality is defined)",
                "action": "store",
                "metavar": "False|undefined|all",
                "required": False,
                "default": False,
                "type": str,
            },
            "--iso_stereo_max_results": {
                "help": "maximum number of stereoisomers to enumerate",
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(isomer_default["stereo_max_results"]),
                "default": isomer_default["stereo_max_results"],
            },
            "--iso_stereo_gen3d": {
                "help": "attempt to generate 3D coordinates for the enumerated molecules",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": type(isomer_default["stereo_gen3d"]),
                "default": isomer_default["stereo_gen3d"],
            },
            #
            # protomers
            #
            "--iso_proto_enum": {
                "help": "enable protomer enumeration",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": type(isomer_default["proto_enum"]),
                "default": isomer_default["proto_enum"],
            },
            "--iso_proto_pH": {
                "help": "specify the pH for the protomer enumeration; if a single value is specified, then all reactions above that pH will be performed (e.g. '--iso_proto_pH 7.4'); if two comma-separated values are specified then all protomer transformation within that range will be peformed (e.g.: '--iso_proto_pH 6.4,8.4)",
                "action": "store",
                "metavar": "pH_value|pH_min,pH_max",
                "required": False,
                "type": type(isomer_default["proto_pH"]),
                "default": isomer_default["proto_pH"],
            },
            "--iso_proto_max_results": {
                "help": "maximum number of protomers to enumerate",
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(isomer_default["proto_max_results"]),
                "default": isomer_default["proto_max_results"],
            },
            "--iso_proto_keep_all": {
                "help": "by default, exaustive protomer generation is performed and only final products are reported; if this flag is used, all intermediate protomers are kept, too",
                "action": "store_true",
                # "metavar": "[ %s ]" % str(isomer_default["proto_keep_all"]),
                "required": False,
                # "type": type(isomer_default["proto_keep_all"]),
                "default": isomer_default["proto_keep_all"],
            },
            "--iso_proto_max_net_charge": {
                "help": "molecules which absolute formal charges exceed this value are discarded",
                "action": "store",
                "metavar": "CHARGE",
                "required": False,
                "type": type(isomer_default["proto_max_net_charge"]),
                "default": isomer_default["proto_max_net_charge"],
            },
            "--iso_proto_neutralize": {
                "help": "generate neutral form of input molecules. WARNING: when used, protomer transformations will be disabled",
                "action": "store_true",
                # "metavar": "",
                "required": False,
                # "type": type(isomer_default["proto_neutralize"]),
                "default": isomer_default["proto_neutralize"],
            },
            #
            # tautomers
            #
            "--iso_tauto_enum": {
                "help": "enable tautomer enumeration",
                "action": "store",
                "metavar": "INT",
                "required": False,
                "type": type(isomer_default["tauto_enum"]),
                "default": isomer_default["tauto_enum"],
            },
            "--iso_tauto_max_results": {
                "help": "max number of tautomers to enumerate",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": type(isomer_default["tauto_max_results"]),
                "default": isomer_default["tauto_max_results"],
            },
            "--iso_tauto_protect_aromatic": {
                "help": "when generating tautomers, discard those reducing the number of aromatic atoms in the system; by setting this option to 'False', less stable tautomers with fewer aromatic atoms will be kept",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": type(isomer_default["tauto_protect_aromatic"]),
                "default": isomer_default["tauto_protect_aromatic"],
            },
            "--iso_tauto_protect_amide": {
                "help": "when generating tautomers, discard those reducing the number of amide groups in the system; by setting this option to 'False', less stable tautomers with fewer amide groups atoms will be kept",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": type(isomer_default["tauto_protect_amide"]),
                "default": isomer_default["tauto_protect_amide"],
            },
            #
            # generic
            #
            "--iso_add_hydrogens": {
                "help": "when enumerating isomers, add hydrogens to the final products",
                "action": "store",
                "metavar": "",
                "required": False,
                "type": type(isomer_default["add_hydrogens"]),
                "default": isomer_default["add_hydrogens"],
            },
            "--iso_max_cycles": {
                "help": "when generating both protomers and tautmers, this option defines how many times the process (protomer+tautomers) is repeated to make sure all products are generated ",
                "action": "store",
                "required": False,
                # "metavar": isomer_default["max_cycles"],
                # "type": type(isomer_default["max_cycles"]),
                "default": isomer_default["max_cycles"],
            },
            "--iso_verbose": {
                "help": "enable verbose mode for all isomeric transformations to track the procress -- WARNING: this can generate a lot of data when processing many molecules",
                "action": "store_true",
                # "metavar": "[ %s ]" % str(isomer_default["verbose"]),
                # "type": type(isomer_default["verbose"]),
                "default": isomer_default["verbose"],
            },
            "--iso_ph_datafile": {
                "help": "specify a custom file for the pH model transformations",
                "action": "store",
                "metavar": "FILENAME",
                "type": str,
                "default": isomer_default["ph_datafile"],
            },
            "--iso_tauto_datafile": {
                "help": "specify a custom file for the tautomer model transformations",
                "action": "store",
                "metavar": "FILENAME",
                "type": str,
                "default": isomer_default["tauto_datafile"],
            },
            "--iso_suppress_rdkit_warnings": {
                "help": "disable warning messages from RDKit when sanitizing molecules",
                "action": "store",
                "metavar": "True|False",
                "required": False,
                "type": type(isomer_default["add_hydrogens"]),
                "default": isomer_default["add_hydrogens"],
            },
        },
    },
    "geometry": {
        "description": "Options the generation of accurate 3D coordinates",
        "values": {
            "--geom_add_h": {
                "help": "add hydrogens to the newly generated 3D structures",
                "action": "store",
                "metavar": "True|False",
                "required": False,
                "type": type(geom_default["add_h"]),
                "default": geom_default["add_h"],
            },
            "--geom_force_trans_amide": {
                "help": "guarantee that secondary amide are in a trans conformation",
                "action": "store",
                "metavar": "True|False",
                "required": False,
                "type": type(geom_default["force_trans_amide"]),
                "default": geom_default["force_trans_amide"],
            },
            "--geom_force_field": {
                "help": "set the force field to be used for the 3D coordinates optimization",
                "action": "store",
                "metavar": "uff|mmff94|mmff94s    [ default: %s ]"
                % geom_default["force_field"],
                "required": False,
                "type": type(geom_default["force_field"]),
                "default": geom_default["force_field"],
            },
            "--geom_max_iterations": {
                "help": "maximum number of iterations for to perform at each minimization cycle",
                "action": "store",
                "metavar": geom_default["max_iterations"],
                "required": False,
                "type": type(geom_default["max_iterations"]),
                "default": geom_default["max_iterations"],
            },
            "--geom_auto_iter_cycles": {
                "help": "set the maximum number of minimization cycles (or attempts) performed on molecules to achieve convergence; the number of iterations at each cycle is defined by the option '--gen_max_iterations'",
                "action": "store",
                "metavar": geom_default["auto_iter_cycles"],
                "required": False,
                "type": type(geom_default["auto_iter_cycles"]),
                "default": geom_default["auto_iter_cycles"],
            },
            "--geom_gen3d": {
                "help": "generate new 3D coordinates prior to optimizing them using the ETKDG implemented in RDKit; if set to False, this option can be used to further minimize molecules with valid 3D coordinates",
                "action": "store",
                "metavar": "True|False    [ default: %s ]" % geom_default["gen3d"],
                "required": False,
                "type": type(geom_default["gen3d"]),
                "default": geom_default["gen3d"],
            },
            "--geom_gen3d_max_attempts": {
                "help": "define how many times initial generation of 3D coordinates is performed, in case of failure",
                "action": "store",
                "metavar": geom_default["gen3d_max_attempts"],
                "required": False,
                "type": type(geom_default["gen3d_max_attempts"]),
                "default": geom_default["gen3d_max_attempts"],
            },
            "--geom_fix_ring_corners": {
                "help": "process 6-membered rings to prevent boat conformatons and identify the optimal chair conformations",
                "action": "store",
                "metavar": "True|False  [default: %s]"
                % str(geom_default["fix_ring_corners"]),
                "required": False,
                "type": type(geom_default["fix_ring_corners"]),
                "default": geom_default["fix_ring_corners"],
            },
            "--geom_preserve_mol_properties": {
                "help": "[SDF ONLY] preserve properties in the input molecule",
                "action": "store",
                "metavar": "True|False  [default: %s]"
                % str(geom_default["preserve_mol_properties"]),
                "required": False,
                "type": type(geom_default["preserve_mol_properties"]),
                "default": geom_default["preserve_mol_properties"],
            },
            "--geom_strict": {
                "help": "reject molecules for which geometry optimization did not converge (to rescue them, increase --geom_auto_iter_cycles and/or --geom_max_iterations)",
                "action": "store",
                "metavar": "True|False  [default: %s]" % str(geom_default["strict"]),
                "required": False,
                "type": type(geom_default["strict"]),
                "default": geom_default["strict"],
            },
        },
    },
    "general": {
        "description": "General options for the calculation setup",
        "values": {
            "--max_proc": {
                "help": "maximum number of processors/cores to use (default: all %d available)"
                % multiprocessing.cpu_count(),
                "action": "store",
                "metavar": "PROC_NUM",
                "required": False,
                "type": type(general_default["max_proc"]),
                "default": general_default["max_proc"],
            },
            "--nice_level": {
                "help": "set the nice level (priority) of the processing; values range from 0 (standard priority; default) to 40 (lowest priority)",
                "action": "store",
                "metavar": "NICE_LEVEL",
                "required": False,
                "type": type(general_default["nice_level"]),
                "default": general_default["nice_level"],
            },
            "--config_file": {
                "help": "specify a JSON config file to be used for the processing; if used, additional command-line options will superseed the config file values; a template can be generated using the '--generate_template_config' option",
                "action": "store",
                "metavar": "config_file",
                "required": False,
                "type": str,
                "default": None,
            },
            "--generate_template_config": {
                "help": "save the full configuration options set into the specified file to be used to create custom config files for the '--config_file' option. If the file exists, the program will cowardly refuse to overwrite it",
                "action": "store",
                "metavar": "config_file",
                "required": False,
                "type": str,
                "default": None,
            },
            "--help_advanced": {
                "help": "print the advanced help",
                "action": "store_true",
                "required": False,
            },
        },
    },
}
