from surpyvor import plots, utils, parse_arguments
import subprocess
import tempfile
import sys
import shlex
import os


def main():
    args = parse_arguments.get_args()
    utils.test_dependencies()
    if args.command == "merge":
        sv_merge(samples=args.variants,
                 distance=args.distance,
                 callers=args.callers,
                 require_type=not args.ignore_type,
                 require_strand=args.strand,
                 estimate_distance=args.estimate_distance,
                 minlength=args.minlength,
                 output=args.output)
    elif args.command == "highsens":
        sv_merge(samples=[utils.vcf_concat(args.variants)],
                 distance=100,
                 callers=1,
                 require_type=True,
                 require_strand=False,
                 estimate_distance=False,
                 minlength=50,
                 output=args.output)
    elif args.command == "highconf":
        sv_merge(samples=args.variants,
                 distance=500,
                 callers=len(args.variants),
                 require_type=True,
                 require_strand=False,
                 estimate_distance=False,
                 minlength=50,
                 output=args.output)
    elif args.command == 'prf':
        precision_recall_fmeasure(args)
    elif args.command == 'upset':
        upset(args)
    elif args.command == 'venn':
        venn(args)
    elif args.command == 'lengthplot':
        lengthplot(args)


def sv_merge(samples, distance, callers, require_type, require_strand,
             estimate_distance, minlength, output):
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
    fhs, interm_out = tempfile.mkstemp(suffix=".vcf")
    with open(fofn_f, 'w') as fofn:
        for s in [utils.decompress(s) for s in samples]:
            fofn.write(s + "\n")
    survivor_cmd = "SURVIVOR merge {fof} {dist} {call} {typ} {str} {estm} {ml} {out}".format(
        fof=fofn_f,
        dist=distance,
        call=callers,
        typ=1 if require_type else -1,
        str=1 if require_strand else -1,
        estm=1 if estimate_distance else -1,
        ml=minlength,
        out=interm_out)
    print("Executing SURVIVOR...", end="", flush=True, file=sys.stderr)
    subprocess.call(shlex.split(survivor_cmd), stdout=subprocess.DEVNULL)
    print("DONE", file=sys.stderr)
    utils.vcf_sort(interm_out, output)
    os.close(fhf)
    os.close(fhs)


def default_merge(args, variants):
    if args.keepmerged:
        vcf_out = args.keepmerged
    else:
        _, vcf_out = tempfile.mkstemp()
    sv_merge(samples=[utils.normalize_vcf(s) for s in variants],
             distance=args.distance,
             callers=1,
             require_type=not args.ignore_type,
             require_strand=False,
             estimate_distance=False,
             minlength=args.minlength,
             output=vcf_out)
    return vcf_out


def precision_recall_fmeasure(args):
    vcf_out = default_merge(args, variants=[args.truth, args.test])
    truth_set, test_set = utils.get_variant_identifiers(vcf=vcf_out,
                                                        ignore_chroms=args.ignore_chroms)

    tp = len(truth_set & test_set)
    precision = tp / len(test_set)
    print(f"Precision: {round(precision, ndigits=4)}")
    recall = tp / len(truth_set)
    print(f"Recall: {round(recall, ndigits=4)}")
    fmeasure = 2 * (precision * recall) / (precision + recall)
    print(f"F-measure: {round(fmeasure, ndigits=4)}")

    if args.venn:
        plots.venn_diagram((truth_set, test_set), labels=('Truth', 'Test'))
    if args.bar:
        plots.bar_chart(vcf_out)
    if args.matrix:
        utils.confusion_matrix(vcf_out, names=['truth', 'test'])


def upset(args):
    vcf_out = default_merge(args, args.variants)
    upsets = utils.make_sets(vcf=vcf_out,
                             names=args.names or args.variants)
    plots.upset_plot(upsets, outname=args.plotout)


def venn(args):
    vcf_out = default_merge(args, args.variants)
    sets = utils.get_variant_identifiers(vcf=vcf_out,
                                         ignore_chroms=[],
                                         num_samples=len(args.variants))
    plots.venn_diagram(sets,
                       labels=args.names or args.variants,
                       num_samples=len(args.variants),
                       outname=args.plotout)


def lengthplot(args):
    len_dict = utils.get_svlengths(args.vcf)
    with open(args.counts, 'w') as counts:
        counts.write("Number of nucleotides affected by SV:\n")
        for svtype, lengths in len_dict.items():
            counts.write("{}:\t{} variants\t{}bp\n".format(
                svtype, len(lengths), sum(lengths)))
    plots.length_plot(dict_of_lengths=len_dict,
                      output=args.plotout)


if __name__ == '__main__':
    main()
