from surpyvor import plots, utils, parse_arguments
import subprocess
import tempfile
import sys
import shlex
import os


def main():
    args = parse_arguments.get_args()

    if args.command == "merge":
        sv_merge(samples=args.files,
                 distance=args.distance,
                 callers=args.callers,
                 type_arg=args.type,
                 strand_arg=args.strand,
                 estimate_distance_arg=args.estimate_distance,
                 minlength=args.minlength,
                 output=args.output)
    elif args.command == "highsens":
        sv_merge(samples=utils.vcf_concat(samples=args.files),
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
                 callers=len(args.files),
                 type_arg=True,
                 strand_arg=False,
                 estimate_distance_arg=False,
                 minlength=50,
                 output=args.output)
    elif args.command == 'prf':
        precision_recall_fmeasure(args)


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
    survivor_cmd = "SURVIVOR merge {fof} {dist} {call} {typ} {str} {estm} {ml} {out}".format(
        fof=fofn_f,
        dist=distance,
        call=callers,
        typ=type,
        str=strand,
        estm=estimate_distance,
        ml=minlength,
        out=output)
    sys.stderr.write("Executing SURVIVOR...\n")
    subprocess.call(shlex.split(survivor_cmd), stdout=subprocess.DEVNULL)
    os.close(fhf)


def precision_recall_fmeasure(args):
    if args.keepmerged:
        vcf_out = args.keepmerged
    else:
        fhv, vcf_out = tempfile.mkstemp()
    sv_merge(samples=[utils.normalize_vcf(s) for s in [args.truth, args.test]],
             distance=args.distance,
             callers=1,
             type_arg=-1 if args.ignore_type else 1,
             strand_arg=-1,
             estimate_distance_arg=-1,
             minlength=args.minlength,
             output=vcf_out)
    truth_set, test_set = utils.get_variant_identifiers(vcf=vcf_out,
                                                        ignore_chroms=args.ignore_chroms)
    plots.venn((truth_set, test_set))
    tp = len(truth_set & test_set)
    precision = tp / len(test_set)
    print(f"Precision: {round(precision, ndigits=4)}")
    recall = tp / len(truth_set)
    print(f"Recall: {round(recall, ndigits=4)}")
    fmeasure = 2*(precision*recall)/(precision + recall)
    print(f"F-measure: {round(fmeasure, ndigits=4)}")
    if args.bar:
        plots.bar_chart(vcf_out)


if __name__ == '__main__':
    main()
