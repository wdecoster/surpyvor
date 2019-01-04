from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from .version import __version__
import subprocess
import tempfile
import sys
import shlex
import os


def main():
    args = get_args()

    if args.command == "merge":
        sv_merge(args.files, args.distance, args.callers, args.type, args.strand,
                 args.estimate_distance, args.minlength, args.output)
    elif args.command == "highsens":
        sv_merge(samples=args.files,
                 distance=100,
                 callers=1,
                 type_arg=True,
                 strand_arg=False,
                 estimate_distance_arg=False,
                 minlength=50,
                 output=args.output)
    elif args.command == "highconf":
        sv_merge(samples=args.files,
                 distance=500,
                 callers=2,
                 type_arg=True,
                 strand_arg=False,
                 estimate_distance_arg=False,
                 minlength=50,
                 output=args.output)
    else:
        print("Unrecognized command")


def sv_merge(samples, distance, callers, type_arg, strand_arg,
             estimate_distance_arg, minlength, output):
    """
    Executes SURVIVOR merge, with parameters:
    -samples.fofn (samples, list)
    -distance between calls (distance, int)
    -number of callers to support call (callers, int)
    -require variants to have sampe type (type, boolean)
    -require variants to be on same strand (strand, boolean)
    -estimate distance between calls (estimate_distance, boolean)
    -specify minimal size of SV event (minlength, int)
    """
    fhf, fofn_f = tempfile.mkstemp()
    with open(fofn_f, 'w') as fofn:
        for s in samples:
            fofn.write(s + "\n")
    type, strand, estimate_distance = (0, 0, 0)
    if type_arg:
        type = 1
    if strand_arg:
        strand = 1
    if estimate_distance_arg:
        estimate_distance = 1
    survivor_cmd = "SURVIVOR merge {} {} {} {} {} {} {} {}".format(
        fofn_f,
        distance,
        callers,
        type,
        strand,
        estimate_distance,
        minlength,
        output)
    sys.stderr.write("Executing SURVIVOR...\n")
    subprocess.call(shlex.split(survivor_cmd), stdout=subprocess.DEVNULL)
    os.close(fhf)


def get_args():
    survivor_version = get_survivor_version()
    parser = ArgumentParser(description="A wrapper around SURVIVOR",
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version",
                        action="version",
                        version='surpyvor: {}, SURVIVOR {}'.format(__version__, survivor_version),
                        help="Print version and quit.")
    subparsers = parser.add_subparsers(help='Available subcommands',
                                       dest='command',
                                       title='[sub-commands]')
    merge = subparsers.add_parser("merge")
    merge.add_argument("-f", "--files",
                       nargs='+',
                       required=True,
                       help="vcf files to merge")
    merge.add_argument("-o", "--output",
                       help="output file",
                       required=True)
    merge.add_argument("-d", "--distance",
                       type=int,
                       default=500,
                       help="distance between variants to merge")
    merge.add_argument("-l", "--minlength",
                       type=int,
                       default=50,
                       help="Minimum length of variants to consider")
    merge.add_argument("-c", "--callers",
                       type=int,
                       default=1,
                       help="Minimum number of callers to support a variant")
    merge.add_argument("-t", "--type",
                       action="store_true",
                       default=False,
                       help="Take type into account")
    merge.add_argument("-s", "--strand",
                       action="store_true",
                       default=False,
                       help="Take strand into account")
    merge.add_argument("-e", "--estimate_distance",
                       action="store_true",
                       default=False,
                       help="Estimate distance between calls")
    highsens = subparsers.add_parser("highsens")
    highsens.add_argument("-f", "--files",
                          nargs='+',
                          required=True,
                          help="vcf files to merge")
    highsens.add_argument("-o", "--output",
                          help="output file",
                          required=True)
    highsens = subparsers.add_parser("highconf")
    highsens.add_argument("-f", "--files",
                          nargs='+',
                          required=True,
                          help="vcf files to merge")
    highsens.add_argument("-o", "--output",
                          help="output file",
                          required=True)
    return parser.parse_args()


def get_survivor_version():
    for line in subprocess.check_output(args="SURVIVOR",
                                        stderr=subprocess.STDOUT,
                                        universal_newlines=True).split('\n'):
        if line.startswith("Version:"):
            return line.strip().split(' ')[1]
    else:
        return "version not found"


if __name__ == '__main__':
    main()
