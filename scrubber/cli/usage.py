

usage = """USAGE GOES HERE
line 1
line 2
line 3
...
"""

epilig = """ EPILOG GOES HERE
line 1
line 2
line 3
...
"""


description = """ DESCRIPTION GOES HERE
line 1
line 2
line 3
...
"""


options = {
    {
        "INPUT/OUTPUT": {
            "desc": "Control input/output data files",
            "opts": {
                "--infile": {
                    "help": "input file; format is guessed from file extension; supported formats are %s ",
                    #'this is the only required argument',
                    "action": "store",
                    "metavar": "INPUT_FILE[.EXT]",
                    "required": True,
                    "type": str,
                },
                '--import' : { 'help': "import from external file formats using OpenBabel; supported formats are %s" }
            },
        }
        # TODO fix missing types

    }
}
