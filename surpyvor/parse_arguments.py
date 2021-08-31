from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from .version import __version__
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
    parent_parser = ArgumentParser(add_help=False)
    parent_parser.add_argument("--verbose",
                               help="Print out more information while running.",
                               action='store_true')
    subparsers = parser.add_subparsers(dest='command',
                                       title='[sub-commands]')
    merge = subparsers.add_parser("merge",
                                  help="merging vcf files of SVs",
                                  parents=[parent_parser])
    merge_req = merge.add_argument_group('required arguments')
    merge_req.add_argument("--variants",
                           nargs='+',
                           required=True,
                           help="vcf files to merge")
    merge_opt = merge.add_argument_group('optional arguments')
    merge_opt.add_argument("-o", "--output",
                           help="output file",
                           default="stdout")
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
                                     parents=[parent_parser])
    highsens_req = highsens.add_argument_group('required arguments')
    highsens_req.add_argument("--variants",
                              nargs='+',
                              required=True,
                              help="vcf files to merge")
    highsens_opt = highsens.add_argument_group('optional arguments')
    highsens_opt.add_argument("-o", "--output",
                              help="output file",
                              default="stdout")
    highsens_opt.add_argument("-d", "--distance",
                              type=int,
                              default=100,
                              help="distance between variants to merge")
    highsens_opt.add_argument("-l", "--minlength",
                              type=int,
                              default=50,
                              help="Minimum length of variants to consider")
    highconf = subparsers.add_parser("highconf",
                                     help="get intersection of SV vcfs",
                                     parents=[parent_parser])
    highconf_req = highconf.add_argument_group('required arguments')
    highconf_req.add_argument("--variants",
                              nargs='+',
                              required=True,
                              help="vcf files to merge")
    highconf_opt = highconf.add_argument_group('optional arguments')
    highconf_opt.add_argument("-o", "--output",
                              help="output file",
                              default="stdout")
    highconf_opt.add_argument("-d", "--distance",
                              type=int,
                              default=500,
                              help="distance between variants to merge")
    highconf_opt.add_argument("-l", "--minlength",
                              type=int,
                              default=50,
                              help="Minimum length of variants to consider")
    highconf_opt.add_argument("-s", "--strand",
                              action="store_true",
                              default=False,
                              help="Take strand into account")
    prf = subparsers.add_parser('prf',
                                help="calculate precision, recall and F-measure",
                                parents=[parent_parser])
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
    prf_opt.add_argument("--venn",
                         help="Make a venn diagram.",
                         action="store_true")

    venn = subparsers.add_parser('venn',
                                 help="Make venn diagram for 2 or 3 SV vcf files",
                                 parents=[parent_parser])
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
                                  parents=[parent_parser])
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

    haplomerge = subparsers.add_parser("haplomerge",
                                       help="merging vcf files of SVs from two haplotypes",
                                       formatter_class=ArgumentDefaultsHelpFormatter,
                                       parents=[parent_parser])
    haplomerge_req = haplomerge.add_argument_group('required arguments')
    haplomerge_req.add_argument("--variants",
                                required=True,
                                nargs="*",
                                help="vcf files to merge")
    haplomerge_opt = haplomerge.add_argument_group('optional arguments')
    haplomerge_opt.add_argument("-o", "--output",
                                help="output file",
                                default="stdout")
    haplomerge_opt.add_argument("-n", "--name",
                                help="name of sample in output VCF",
                                default="stdout")
    haplomerge_opt.add_argument("-d", "--distance",
                                type=int,
                                default=200,
                                help="distance between variants to merge")
    haplomerge_opt.add_argument("-l", "--minlength",
                                type=int,
                                default=50,
                                help="Minimum length of variants to consider")
    haplomerge_opt.add_argument("-c", "--callers",
                                type=int,
                                default=1,
                                help="Minimum number of callers to support a variant")
    haplomerge_opt.add_argument("-i", "--ignore_type",
                                help="Ignore the type of the structural variant",
                                action="store_true",
                                default=False)
    haplomerge_opt.add_argument("-s", "--strand",
                                action="store_true",
                                default=False,
                                help="Take strand into account")
    haplomerge_opt.add_argument("-e", "--estimate_distance",
                                action="store_true",
                                default=False,
                                help="Estimate distance between calls")

    lengthplot = subparsers.add_parser('lengthplot',
                                       help="create stacked bar plot of SV lengths split by type",
                                       parents=[parent_parser])
    lengthplot_req = lengthplot.add_argument_group('required arguments')
    lengthplot_req.add_argument("vcf", help="vcf file to parse")
    lengthplot_opt = lengthplot.add_argument_group('optional arguments')
    lengthplot_opt.add_argument("--plotout",
                                help="output file to write figure to",
                                default="SV-length.png")
    lengthplot_opt.add_argument("-c", "--counts",
                                help="output file to write counts to",
                                default="SV-length.txt")

    minlength = subparsers.add_parser('minlen',
                                      help="filter a SV vcf file by minimal variant length",
                                      parents=[parent_parser])
    minlength_req = minlength.add_argument_group('required arguments')
    minlength_req.add_argument("vcf", help="vcf file to parse")

    minlength_opt = minlength.add_argument_group('optional arguments')
    minlength_opt.add_argument("-l", "--length",
                               help="minimal SV length",
                               type=int,
                               default=50)
    minlength_opt.add_argument("-o", "--output", help="vcf file to write to", default=None)

    truncate_svlen = subparsers.add_parser('svlentruncate',
                                           help="limit the SVLEN to a certain (positive) length",
                                           parents=[parent_parser])
    truncate_svlen_req = truncate_svlen.add_argument_group('required arguments')
    truncate_svlen_req.add_argument("vcf", help="vcf file to parse")

    truncate_svlen_opt = truncate_svlen.add_argument_group('optional arguments')
    truncate_svlen_opt.add_argument("-l", "--length",
                                    help="maximal SVLEN, replace SVLEN by this value if larger",
                                    type=int,
                                    default=1e5)
    truncate_svlen_opt.add_argument("-o", "--output", help="vcf file to write to", default=None)

    fixvcf = subparsers.add_parser('fixvcf',
                                   help="Some fixes to make compatible with e.g. vcfanno",
                                   parents=[parent_parser])
    fixvcf_req = fixvcf.add_argument_group('required arguments')
    fixvcf_req.add_argument("vcf", help="vcf file to parse")
    fixvcf_req.add_argument("--fai", help="index of corresponding fasta file", required=True)
    fixvcf_opt = fixvcf.add_argument_group('optional arguments')
    fixvcf_opt.add_argument("-o", "--output", help="vcf file to write to", default=None)

    purge2d = subparsers.add_parser('purge2d',
                                    help="Remove accidental 2D reads from a bam file",
                                    parents=[parent_parser])
    purge2d_req = purge2d.add_argument_group('required arguments')
    purge2d_req.add_argument("bam", help="bam file to filter")

    purge2d_opt = purge2d.add_argument_group('optional arguments')
    purge2d_opt.add_argument("-o", "--output",
                             help="sam/bam file to write filtered alignments to [stdout]",
                             default="-")

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
    if args.command == 'haplomerge':
        if not len(args.variants) in [2, 3]:
            sys.exit("INPUT ERROR: "
                     "haplomerge can only be used on 2 or 3 vcf files!")
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
    import subprocess

    for line in subprocess.check_output(args="SURVIVOR",
                                        stderr=subprocess.STDOUT,
                                        universal_newlines=True).split('\n'):
        if line.startswith("Version:"):
            return line.strip().split(' ')[1]
    else:
        return "version not found"
