from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from .version import __version__
import subprocess
import sys


def get_args():
    parser = ArgumentParser(description="A wrapper around SURVIVOR, with convenience functions",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version",
                        action="version",
                        version='surpyvor: {}, SURVIVOR {}'.format(
                            __version__, get_survivor_version()),
                        help="Print version and quit.")
    subparsers = parser.add_subparsers(dest='command',
                                       title='[sub-commands]')
    merge = subparsers.add_parser("merge", help="merging vcf files of SVs")
    merge_req = merge.add_argument_group('required arguments')
    merge_req.add_argument("-f", "--files",
                           nargs='+',
                           required=True,
                           help="vcf files to merge")
    merge_req.add_argument("-o", "--output",
                           help="output file",
                           required=True)
    merge_opt = merge.add_argument_group('optional arguments')
    merge_opt.add_argument("-d", "--distance",
                           type=int,
                           default=500,
                           help="distance between variants to merge")
    merge_opt.add_argument("-l", "--minlength",
                           type=int,
                           default=50,
                           help="Minimum length of variants to consider")
    merge_opt.add_argument("-c", "--callers",
                           type=int,
                           default=1,
                           help="Minimum number of callers to support a variant")
    merge_opt.add_argument("-t", "--type",
                           action="store_true",
                           default=False,
                           help="Take type into account")
    merge_opt.add_argument("-s", "--strand",
                           action="store_true",
                           default=False,
                           help="Take strand into account")
    merge_opt.add_argument("-e", "--estimate_distance",
                           action="store_true",
                           default=False,
                           help="Estimate distance between calls")
    highsens = subparsers.add_parser("highsens", help="get union of SV vcfs")
    highsens_req = highsens.add_argument_group('required arguments')
    highsens_req.add_argument("-f", "--files",
                              nargs='+',
                              required=True,
                              help="vcf files to merge")
    highsens_req.add_argument("-o", "--output",
                              help="output file",
                              required=True)
    highconf = subparsers.add_parser("highconf", help="get intersection of SV vcfs")
    highconf_req = highconf.add_argument_group('required arguments')
    highconf_req.add_argument("-f", "--files",
                              nargs='+',
                              required=True,
                              help="vcf files to merge")
    highconf_req.add_argument("-o", "--output",
                              help="output file",
                              required=True)
    prf = subparsers.add_parser('prf', help="calculate precision, recall and F-measure")
    prf_req = prf.add_argument_group('required arguments')
    prf_req.add_argument("--truth",
                         help="vcf containing truth set",
                         required=True)
    prf_req.add_argument("--test",
                         help="vcf containing test set",
                         required=True)
    prf_opt = prf.add_argument_group('optional arguments')
    prf_opt.add_argument("-d", "--distance",
                         help="maximum distance between test and truth call",
                         default=500)
    prf_opt.add_argument("--minlength",
                         help="Minimum length of SVs to be taken into account",
                         default=50)
    prf_opt.add_argument("-i", "--ignore_type",
                         help="Ignore the type of the structural variant",
                         action="store_true")
    prf_opt.add_argument("--bar",
                         help="Make stacked bar chart of SV lengths coloured by validation status",
                         action="store_true")
    prf_opt.add_argument("--ignore_chroms",
                         help="Chromosomes to ignore for prf calculation.",
                         nargs='*',
                         default=['chrEBV'])

    args = parser.parse_args()
    if not args.command:
        sys.stderr.write("INPUT ERROR: sub-command required\n\n")
        parser.print_help()
        sys.exit()
    return args


def get_survivor_version():
    for line in subprocess.check_output(args="SURVIVOR",
                                        stderr=subprocess.STDOUT,
                                        universal_newlines=True).split('\n'):
        if line.startswith("Version:"):
            return line.strip().split(' ')[1]
    else:
        return "version not found"
