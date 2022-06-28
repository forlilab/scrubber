description = """
DESCRIPTION
This is the Raccoon Scrubber. This is the new version based on RDKit
"""

usage = """

In order to get something done with the Scrubber this is the command:

    $ %s [options] filename

""" % (
    "scrubber.py"
)

epilog = """

(C) 2022 ForliLab Industries, Scripps Research

"""


advanced_help = """

This is the advanced help

                      _(\-/)_
                     {(#b^d#)}
                     `-.(Y).-`

This is where the full pipeline is described (scrubber.core.process())
"""


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
                "metavar": "",
                "required": False,
                "type": bool,
                "default": True,
            },
            "--in_removeHs": {
                "help": "remove hydrogens from input structure, if present",
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": bool,
                "default": False,
            },
            "--in_strictParsing": {
                "help": "perform strict parsing with RDKit; molecules will be skipped if any inaccuracies are found [CLARIFY THIS!]",
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": bool,
                "default": True,
            },
            "--in_name_property": {
                "help": "[SDF only] use the specified property field to name the molecule. If the --out_split option is used, this will be the filename.",
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": str,
                "default": None,
            },
            "--in_safeparsing": {
                "help": "parse input files in safe mode, instead of using RDKit native parsers; enable saving problematic raw input text into log files; it can be slower than using unsafe parsing",
                "action": "store",
                # "metavar": "[.EXT]",
                "required": False,
                "type": bool,
                "default": True,
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
                "metavar": "'single'|'split'",
                "required": False,
                "type": str,
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
            # "--out_ftype": {
            #     "help": """specify the output file format; by default the file
            # extension will be used to set the file format, unless the option
            # \'--out_ftype\' is used; if a fullpath with explicit directories
            # will be specified, the directories will be created automatically,
            # if necessary """,
            #     "action": "store",
            #     "metavar": "OUTPUT_FNAME[.EXT]",
            #     "required": True,
            #     "type": str,
            # },
            "--out_ftype": {
                "help": "specify the format of the output file, overriding the extension (i.e., writing SMILES files with '.txt' extension); the format has to be one of the supported types (SMI or SDF).",
                "action": "store",
                "metavar": "EXT",
                "required": False,
                "type": str,
            },
            "--out_format_opts": {
                "help": "format-specific options: SMI and SDF (some stuff from RDKIT?)",
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
        "description": "Options for pre-filter settings (i.e. before any processing occurs)",
        "values": {
            "--filtpre_mw_max": {
                "help": "specify the format of the input file, overriding the extension (i.e., parsing SMILES files with '.txt' extension); the format has to be one of the supported types (SMI or SDF); NOTE: if a molecule does not have a name defined, the \"MOL\" name will be set by default",
                "action": "store",
                "metavar": "EXT",
                "required": False,
                "type": str,
            },
        },
    },
}
