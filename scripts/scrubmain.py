import sys
import argparse
import json
import os


# DEBUG
from pprint import pprint as pp

sys.path.append("../")
try:
    from scrubber.core import ScrubberCore
except Exception as e:
    print("The script must be executed from the script/ directory")
    raise e
# import cli
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

# from scrubber.core import ScrubberCore, ScrubberClass


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


class ScrubbCLI(object):
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

    def init_opt_parser(self):
        """initialize the parser"""
        self.parser = argparse.ArgumentParser(
            description=description,
            usage=usage,
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
            for name, opt in values.items():
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
        cli_options = {}
        file_options = {}
        extra_actions = {
            "config_file": None,
            "save_config_template": None,
            "help_advanced": None,
        }
        print("XXXX INTERESTING TERST")
        pp(vars(self.args))
        for (
            opt,
            value,
        ) in self.args.__dict__.items():
            # process extra options
            if opt in extra_actions:
                extra_actions[opt] = value
            # process standard options
            else:
                tag, kw = opt.split("_", 1)
                if tag in tags_reverse:
                    group = tags_reverse[tag]
                    # opt is still original
                else:
                    group = "general"
                    opt = kw
                if not group in cli_options:
                    cli_options[group] = {}
                cli_options[group][opt] = value
        # check that no more than one special option is used at the same time
        print("EXTRA", extra_actions)
        check = sum( [ int(v is not None) for _, v in extra_actions.items()] )
        print("CHECK", check)
        if  check  > 1:
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
            # if True:
                with open(fname, "w") as fp:
                    json.dump(ScrubberCore.get_defaults(terse=True), fp, indent=4)
            except Exception as e:
                print(
                    "*** ERROR *** impossible to write the "
                    "file [%s]. Error message: %s" % (fname, e.__str__())
                )
                sys.exit(1)
            sys.exit(0)
        # read config file
        if not extra_actions["config_file"] is None:
            fname = extra_actions["config_file"]
            print("FNAME IS", fname)
            try:
                with open(fname, "r") as fp:
                    options = json.load(fp)
            except Exception as e:
                print(
                    "*** ERROR *** impossible to read the "
                    "file [%s]. Error message: %s" % (fname, e.__str__())
                )
                sys.exit(1)
        # set file options first
        for group, values in file_options.items():
            if not group in self.options:
                print(
                        "*** ERROR *** unknown option type in the file: [%s]" % group
                )
                sys.exit(1)
            for opt, v in values.items():
                if not opt in self.options[group]['values']:
                    print(
                            "*** ERROR *** unknown option in the file: [%s]" % opt
                    )
                    sys.exit(1)
                self.options[group]["values"][opt] = v
        # no other CLI options are allowed if --config_file is used
        if len(sys.argv)>3:
            print(
                    "*** ERROR *** \"--config_file\" option cannot be used with other command-line options"
            )
            sys.exit(1)

        # # set CLI options after files
        # for group, values in cli_options.items():
        #     for opt, v in values.items():
        #         if group in file_options:
        #             if
        #         self.options[group]["values"][opt] = v
        # print("XXXXX")# return cli_options, file_options

    def show_advanced_help(self):
        """show advanced help and exit"""
        print(advanced_help)
        sys.exit(0)


if __name__ == "__main__":
    from pprint import pprint as pp
    ss = ScrubbCLI()
    ss.init_opt_parser()
    ss.process_args()
    print("=======")
    # pp(ss.options)
