import sys
import argparse
import json
import os


# DEBUG
from scrubber.core import ScrubberCore, MoleculeIssueStorage

from scrubber.cli import (
    description,
    usage,
    epilog,
    advanced_help,
    cli_options,
    tags,
    tags_reverse,
    extra_options,
)


class SmartFormatter(argparse.HelpFormatter):
    """custom argparse class for smart formatting (including raw text)
    # http://stackoverflow.com/questions/3853722/python-argparse-
    # how-to-insert-newline-the-help-text
    """

    # def __init__(self, prog, indent_increment=1, max_help_position=10, width=80):
    #     argparse.HelpFormatter.__init__(
    #         self,
    #         prog=prog,
    #         indent_increment=indent_increment,
    #         max_help_position=max_help_position,
    #         width=width
    #     )

    def _split_lines(self, text, width=60):
        # this is the RawTextHelpFormatter._split_lines
        if text.startswith("R|"):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)


class ScrubberCLI:
    """This class creates the working instance of the CLI interface of the Scrubber.
    The functions performed are:
        - generate CLI options for all options 'published' by the ScrubberCore
        - add options for CLI-specific functions:
            - read a JSON config file
            - write a JSON config file with all the defaults of the ScrubberCore
            - add a "--help_advanced" CLI option to show advanced documentation
        - parse options passed at CLI and make sure there are consistent
    """

    def __init__(self):
        """something"""
        self.options = ScrubberCore.get_defaults()
        self.init_opt_parser()
        self.quiet = False

    def init_opt_parser(self):
        """initialize the parser"""
        self.parser = argparse.ArgumentParser(
            description=description,
            # usage=usage,
            epilog=epilog,
            # formatter_class=argparse.RawDescriptionHelpFormatter,
            # formatter_class=SmartFormatter,
        )
        for group_name, group_info in self.options.items():
            values = group_info["values"]
            ignore = group_info["ignore"]
            desc = cli_options[group_name]["description"]
            group = self.parser.add_argument_group(group_name.upper(), desc)
            options = cli_options[group_name]["values"]
            # TODO add a line to check if active or not?
            for name, _ in values.items():
                if name in ignore:
                    continue
                tag = tags[group_name]
                tag_name = "--%s%s" % (tag, name)
                arg_settings = options[tag_name]
                group.add_argument(tag_name, **arg_settings)
            # add custom options for general group
            if group_name == "general":
                for k, v in extra_options.items():
                    group.add_argument(k, **v)
                # group.add_argument("--config_file", **options["--config_file"])
                # group.add_argument(
                #     "--save_config_template", **options["--save_config_template"]
                # )
                # group.add_argument("--help_advanced", **options["--help_advanced"])
        # help advanced checking
        # if "--help_advanced" in sys.argv:
        #     self.show_advanced_help()
        # self.parser.print_usage()
        self.process_args()

    def process_args(self):
        """process the options passed to the command line, with higher priority
        to the extra options and generate the ScrubberCore options dictionary
        from the argparse namespace.
        Extra options have higher priority and need to be processed as first because they can
        either be overridden by CLI options (--config_file), or will trigger an
        early termination of the processing (--help_advanced and
        --save_config_template)"""
        self.args = self.parser.parse_args()
        cli_opts = {}
        file_opts = {}
        extra_actions = {
            "config_file": None,
            "save_config_template": None,
            "help_advanced": None,
        }
        # print("XXXX INTERESTING TERST")
        # print(vars(self.args))
        for (
            opt,
            value,
        ) in self.args.__dict__.items():
            # process extra options
            if opt in extra_actions:
                extra_actions[opt] = value
            # process standard options
            else:
                group = "general"
                if "_" in opt:
                    tag, kw = opt.split("_", 1)
                    if tag in tags_reverse:
                        group = tags_reverse[tag]
                        opt = kw
                if not group in cli_opts:
                    cli_opts[group] = {}
                cli_opts[group][opt] = value
                if opt == "quiet":
                    self.quiet = value

        # check that no more than one special option is used at the same time
        # print("CLIOPTS", cli_options)
        check = sum([int(v is not None) for v in extra_actions.values()])
        # print("CHECK", check)
        if check > 1:
            print(
                "*** ERROR *** conflicting options used. "
                "Use one of either --config_file, "
                "--save_config_template, --help_advanced."
            )
            sys.exit(1)
        # help advanced
        if not extra_actions["help_advanced"] is None:
            print("ADVANCED")
            self.show_advanced_help()
            sys.exit(0)
        # save config template
        if not extra_actions["save_config_template"] is None:
            fname = extra_actions["save_config_template"]
            if os.path.exists(fname):
                print(
                    "*** ERROR *** cannrt save config template, "
                    "file [%s] already exists. Cowardly refusing "
                    "to overwrite it." % fname
                )
                sys.exit(1)
            try:
                with open(fname, "w") as fp:
                    # remove the "values"  field
                    data = ScrubberCore.get_defaults(terse=True)
                    for k, v in data.items():
                        data[k] = v["values"]
                    json.dump(ScrubberCore.get_defaults(terse=True), fp, indent=4)
            except Exception as e:
                print(
                    "*** ERROR *** impossible to write the "
                    "file [%s]. Error message: %s" % (fname, str(e))
                )
                sys.exit(1)
            sys.exit(0)
        # read config file
        if not extra_actions["config_file"] is None:
            fname = extra_actions["config_file"]
            print("FNAME IS", fname)
            try:
                with open(fname, "r") as fp:
                    file_opts = json.load(fp)
                # restore the "values" field
                for k, v in file_opts.items():
                    file_opts[k] = {"values": v}
            except Exception as e:
                print(
                    "*** ERROR *** impossible to read the "
                    "file [%s]. Error message: %s" % (fname, str(e))
                )
                sys.exit(1)
        # set file options first
        for group, opts in file_opts.items():
            if not group in self.options:
                print("*** ERROR *** unknown option type in the file: [%s]" % group)
                sys.exit(1)
            for key, value in opts["values"].items():
                if not key in self.options[group]["values"]:
                    print("*** ERROR *** unknown option in the file: [%s]" % key)
                    sys.exit(1)
                    print("setting option,", group, key, value)
                self.options[group]["values"][key] = value
        # set CLI options after files
        # strip the CLI tag to extract the kw, e.g.:
        # --in_fname -> "in_fname" -> ("in_" in tags_reverse): True -> kw = "fname"
        # --max_proc -> "max_proc" -> ["max_" in tags_reverse): False - > kw = "max_proc"]
        for group, values in cli_opts.items():
            for opt, v in values.items():
                # print("RAW:", opt, v)
                if "_" in opt:
                    tag, kw = opt.split("_", 1)
                    if tag in tags_reverse:
                        opt = kw
                self.options[group]["values"][opt] = v
        # after parsing  both CLI and config file options, check that required
        # options are defined:
        # -fname (input)
        # -fname (output)
        missing_options = []
        if self.options["input"]["values"]["fname"] is None:
            missing_options.append("--%s%s" % (tags["input"], "fname"))
        if self.options["output"]["values"]["fname"] is None:
            missing_options.append("--%s%s" % (tags["output"], "fname"))
        if len(missing_options) > 0:
            if extra_actions["config_file"] is None:
                print(
                    "ERROR: required options are missing: %s"
                    % (", ".join(missing_options))
                )
                sys.exit(1)
            else:
                print(
                    "ERROR: input and/or output files are not specified; "
                    "either add them to the config file or set them "
                    "using the following options: %s" % (", ".join(missing_options))
                )
                sys.exit(1)
        # check that files do not exist already
        if os.path.exists(self.options["output"]["values"]["fname"]):
            if not cli_opts["general"]["overwrite"]:
                print(
                    "ERROR: the output file [ %s ] exists. Specify a different "
                    "filename or use the '--overwrite' option"
                    % (self.options["output"]["values"]["fname"])
                )
                sys.exit(1)
        found = []
        if not cli_opts["general"]["overwrite"]:
            for f in MoleculeIssueStorage.generate_filenames(
                    self.options["errors"]["values"]["log_basename"]
                    ):
                if os.path.exists(f):
                    found.append(f)
            if len(found):
                print(
                    "ERROR: the following log files already exist. Specify a different "
                    "basename (--err_log_basename) or use the '--overwrite' option:\n%s"
                    % ("\n".join(found)))
                sys.exit(1)

        self.core = ScrubberCore(self.options)

    def start(self):
        """start processing the files specified in the CLI options"""
        self.core.process_file(self.quiet)
        # if self.options['quiet']:
        #     print("QUIET")
        #     return
        # print(self.core.queue_err.close())

    def show_advanced_help(self):
        """show advanced help and exit"""
        print(advanced_help)

        sys.exit(0)


if __name__ == "__main__":
    scrub = ScrubberCLI()
    scrub.start()
