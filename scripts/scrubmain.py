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
        }
        self._tag_reverse = {v: k for k, v in self._tags.items()}
        self._init_opt_parser()

    def _init_opt_parser(self):
        """initialize the parser"""
        # http://stackoverflow.com/questions/3853722/python-argparse-
        # how-to-insert-newline-the-help-text
        class SmartFormatter(argparse.HelpFormatter):
            def _split_lines(self, text, width):
                # this is the RawTextHelpFormatter._split_lines
                if text.startswith("R|"):
                    return text[2:].splitlines()
                return argparse.HelpFormatter._split_lines(self, text, width)

        self.parser = argparse.ArgumentParser(
            description=description,
            usage=usage,
            epilog=epilog,
            formatter_class=argparse.RawDescriptionHelpFormatter,
        )

        try:
            for group_name, group_info in self.defaults.items():
                # desc = group_info["desc"]
                values = group_info["values"]
                ignore = group_info["ignore"]
                desc = cli_options[group_name]["description"]
                group = self.parser.add_argument_group(group_name.upper(), desc)
                # TODO add a line to check if active or not?
                for name, opt in values.items():
                        if name in ignore:
                            print("[ignoring: %s]" % name)
                            continue
                        print("NAME_VALUES", group_name, name, opt)
                        tag = self._tags[group_name]
                        tag_name = "--%s%s" % (tag, name)
                        print("TAGNAME", tag_name)
                        opts = cli_options[group_name]["values"][tag_name]
                        # TODO default values are updated here
                        # opts = self._args_dict[name]
                        group.add_argument(tag_name, **opts)
        except:
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
