import sys
import argparse

sys.path.append("../")
try:
    from scrubber.core import ScrubberCore
except Exception as e:
    print("The script must be executed from the script/ directory")
    raise e
# import cli
from scrubber.cli import description, usage, epilog, advanced_help, cli_options

# from scrubber.core import ScrubberCore, ScrubberClass

config_file_opts ={"help": "specify a JSON config file "}

class ScrubbCLI(object):
    """testing"""

    def __init__(self):
        """something"""
        self.defaults = ScrubberCore.get_defaults()
        self._tags = {
            "input": "in_",
            "output": "out_",
            "isomers": "iso_",
            "geometry": "geom_",
            "log": "log_",
            "filter_pre": "filtpre_",
            "filter_post": "filtpost_",
            "reaction": "react_",
            "general": "",
        }



        self._tag_reverse = {v: k for k, v in self._tags.items()}
        self._init_opt_parser()

    def _init_opt_parser(self):
        """initialize the parser"""
        # http://stackoverflow.com/questions/3853722/python-argparse-
        # how-to-insert-newline-the-help-text
        class SmartFormatter(argparse.HelpFormatter):
            def __init__(
                self, prog, indent_increment=1, max_help_position=10, width=80
            ):
                argparse.HelpFormatter.__init__(
                    self,
                    prog=prog,
                    indent_increment=indent_increment,
                    max_help_position=max_help_position,
                    width=width,
                )

            def _split_lines(self, text, width=60):
                # this is the RawTextHelpFormatter._split_lines
                if text.startswith("R|"):
                    return text[2:].splitlines()
                return argparse.HelpFormatter._split_lines(self, text, width)

        self.parser = argparse.ArgumentParser(
            description=description,
            usage=usage,
            epilog=epilog,
            # formatter_class=argparse.RawDescriptionHelpFormatter,
            formatter_class=SmartFormatter,
        )
        # if True:
        try:
            print("CLI OPTIONS", cli_options.keys())
            for group_name, group_info in self.defaults.items():
                # desc = group_info["desc"]
                values = group_info["values"]
                ignore = group_info["ignore"]
                desc = cli_options[group_name]["description"]
                group = self.parser.add_argument_group(group_name.upper(), desc)
                options = cli_options[group_name]['values']
                # TODO add a line to check if active or not?
                for name, opt in values.items():
                    if name in ignore:
                        # print("[ignoring: %s]" % name)
                        continue
                    # print("NAME_VALUES", group_name, name, opt)
                    tag = self._tags[group_name]
                    tag_name = "--%s%s" % (tag, name)
                    # print("TAGNAME", tag_name)
                    arg_settings = options[tag_name]
                    # TODO default values are updated here
                    # opts = self._args_dict[name]
                    # print("OOPTS", opts)
                    group.add_argument(tag_name, **arg_settings)
                # add custom options for general group
                if group_name =='general':
                    group.add_argument("--config_file", **options['--config_file'])
                    group.add_argument("--help_advanced", **options['--help_advanced'])
        except Exception as e:
            print("EXCEPTION IS", e)
            pass
        # help advanced checking
        if "--help_advanced" in sys.argv:
            self.show_advanced_help()
            sys.exit(0)
        self.args = self.parser.parse_args()
        # activate verbose as soon as possible, if requested
        # self.verbose = self.args.verbose
        self._data = []
        print("CKECKING", self.args)
        # self.parser.print_help()
        self.parser.print_usage()
        # sys.exit(0)

    def show_advanced_help(self):
        """show advanced help and exit"""
        print(advanced_help)
        sys.exit(0)


if __name__ == "__main__":
    ScrubbCLI()
