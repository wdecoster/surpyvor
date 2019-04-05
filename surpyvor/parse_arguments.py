from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from .version import __version__
import subprocess
import sys
from os import path


def get_args():
    parser = ArgumentParser(description="A wrapper around SURVIVOR, with convenience functions.",
                            formatter_class=ArgumentDefaultsHelpFormatter,)
    parser.add_argument("-v", "--version",
                        action="version",
                        version='surpyvor: {}, SURVIVOR {}'.format(
                            __version__, get_survivor_version()),
                        help="Print version and quit.")
    subparsers = parser.add_subparsers(dest='command',
                                       title='[sub-commands]')
    merge = subparsers.add_parser("merge",
                                  help="merging vcf files of SVs",
                                  formatter_class=ArgumentDefaultsHelpFormatter)
    merge_req = merge.add_argument_group('required arguments')
    merge_req.add_argument("--variants",
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
    merge_opt.add_argument("-i", "--ignore_type",
                           help="Ignore the type of the structural variant",
                           action="store_true",
                           default=False)
    merge_opt.add_argument("-s", "--strand",
                           action="store_true",
                           default=False,
                           help="Take strand into account")
    merge_opt.add_argument("-e", "--estimate_distance",
                           action="store_true",
                           default=False,
                           help="Estimate distance between calls")
    highsens = subparsers.add_parser("highsens",
                                     help="get union of SV vcfs",
                                     formatter_class=ArgumentDefaultsHelpFormatter)
    highsens_req = highsens.add_argument_group('required arguments')
    highsens_req.add_argument("--variants",
                              nargs='+',
                              required=True,
                              help="vcf files to merge")
    highsens_req.add_argument("-o", "--output",
                              help="output file",
                              required=True)
    highconf = subparsers.add_parser("highconf",
                                     help="get intersection of SV vcfs",
                                     formatter_class=ArgumentDefaultsHelpFormatter)
    highconf_req = highconf.add_argument_group('required arguments')
    highconf_req.add_argument("--variants",
                              nargs='+',
                              required=True,
                              help="vcf files to merge")
    highconf_req.add_argument("-o", "--output",
                              help="output file",
                              required=True)
    prf = subparsers.add_parser('prf',
                                help="calculate precision, recall and F-measure",
                                formatter_class=ArgumentDefaultsHelpFormatter)
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
                         action="store_true",
                         default=False)
    prf_opt.add_argument("--ignore_chroms",
                         help="Chromosomes to ignore for prf calculation.",
                         nargs='*',
                         default=['chrEBV'])
    prf_opt.add_argument("--keepmerged",
                         help="Save merged vcf file.",
                         default=False)
    prf_opt.add_argument("--bar",
                         help="Make stacked bar chart of SV lengths coloured by validation status",
                         action="store_true")
    prf_opt.add_argument("--matrix",
                         help="Make a confusion matrix.",
                         action="store_true")
    venn = subparsers.add_parser('venn',
                                 help="Make venn diagram for 2 or 3 SV vcf files",
                                 formatter_class=ArgumentDefaultsHelpFormatter)
    venn_req = venn.add_argument_group('required arguments')
    venn_req.add_argument("--variants",
                          help="vcfs containing structural variants",
                          required=True,
                          nargs="*")
    venn_opt = venn.add_argument_group('optional arguments')
    venn_opt.add_argument("--names",
                          help="Names of datasets in --variants",
                          nargs="*")
    venn_opt.add_argument("-d", "--distance",
                          help="maximum distance between test and truth call",
                          default=500)
    venn_opt.add_argument("--minlength",
                          help="Minimum length of SVs to be taken into account",
                          default=50)
    venn_opt.add_argument("-i", "--ignore_type",
                          help="Ignore the type of the structural variant",
                          action="store_true",
                          default=False)
    venn_opt.add_argument("--keepmerged",
                          help="Save merged vcf file")
    venn_opt.add_argument("--plotout",
                          help="Name of output plot",
                          default="venn.png")
    upset = subparsers.add_parser('upset',
                                  help="Make upset plot for multiple SV vcf files",
                                  formatter_class=ArgumentDefaultsHelpFormatter)
    upset_req = upset.add_argument_group('required arguments')
    upset_req.add_argument("--variants",
                           help="vcfs containing structural variants",
                           required=True,
                           nargs="*")
    upset_opt = upset.add_argument_group('optional arguments')
    upset_opt.add_argument("--names",
                           help="Names of datasets in --variants",
                           nargs="*")
    upset_opt.add_argument("-d", "--distance",
                           help="maximum distance between test and truth call",
                           default=500)
    upset_opt.add_argument("--minlength",
                           help="Minimum length of SVs to be taken into account",
                           default=50)
    upset_opt.add_argument("-i", "--ignore_type",
                           help="Ignore the type of the structural variant",
                           action="store_true",
                           default=False)
    upset_opt.add_argument("--keepmerged",
                           help="Save merged vcf file")
    upset_opt.add_argument("--plotout",
                           help="Name of output plot",
                           default="UpSetPlot.png")
    args = parser.parse_args()
    validate_args(parser, args)
    return args


def validate_args(parser, args):
    if not args.command:
        sys.stderr.write("INPUT ERROR: sub-command required\n\n")
        parser.print_help()
        sys.exit()
    if args.command in ['upset', 'venn']:
        if args.names:
            if not len(args.variants) == len(args.names):
                sys.exit("INPUT ERROR: "
                         "Need to have same number of values in --names as --variants!")
    if args.command == 'venn':
        if len(args.variants) > 3:
            sys.exit("INPUT ERROR: "
                     "Venn diagrams are only created for 2 or 3 vcf files!")
    if hasattr(args, 'variants'):
        for f in args.variants:
            if not path.isfile(f):
                sys.exit("File not found: {}".format(f))
    if hasattr(args, 'truth'):
        if not path.isfile(args.truth):
            sys.exit("File not found: {}".format(args.truth))
    if hasattr(args, 'test'):
        if not path.isfile(args.test):
            sys.exit("File not found: {}".format(args.test))


def get_survivor_version():
    for line in subprocess.check_output(args="SURVIVOR",
                                        stderr=subprocess.STDOUT,
                                        universal_newlines=True).split('\n'):
        if line.startswith("Version:"):
            return line.strip().split(' ')[1]
    else:
        return "version not found"
