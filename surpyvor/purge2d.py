import pysam
import sys
from cigar import Cigar
import pandas as pd


def process(bamfile, output="-", write_candidates=False):
    twod = []
    bam = pysam.AlignmentFile(bamfile, "rb")
    mode = 'wb' if output.endswith('.bam') else 'w'
    filtered_alignments = pysam.AlignmentFile(output, mode, template=bam)
    for read in bam.fetch():
        if is_accidental_2d(read):
            twod.append(read)
        else:
            filtered_alignments.write(read)
    sys.stderr.write(f"Detected {len(twod)} potential artefacts "
                     f"out of {bam.mapped} alignments ({100*(len(twod))/bam.mapped}%)\n")
    if write_candidates:
        twod_bam = pysam.AlignmentFile("2D-candidates.bam", "wb", template=bam)
        for r in twod:
            twod_bam.write(r)
    keep_these_candidates = follow_up_2d_candidates(twod, distance=500)
    if keep_these_candidates:
        sys.stderr.write("WARNING: Some potential artefacts are close to another.\n")
        sys.stderr.write("WARNING: As this could be an SV, these reads are kept.\n")
        sys.stderr.write("WARNING: Surpyvor is producing an UNSORTED sam/bam.\n")
        for read in keep_these_candidates:
            filtered_alignments.write(read)


def get_sa_attributes(sa_tag):
    '''
    An entry in the SA tag consist of rname, POS, strand, CIGAR, mapQ, NM
    Returned are POS, strand and aligned length (based on cigar)
    '''
    sasplit = sa_tag.split(',')
    start = int(sasplit[1])
    end = start + Cigar(sasplit[3]).reference_length()
    strand = sasplit[2]
    return start, end, strand


def get_strand(read):
    return '-' if read.is_reverse else '+'


def reverse_overlaps(sa_start, sa_end, sa_strand, read):
    """
    Return True for reads overlapping their supplementary alignment on the reverse strand
    """
    return sa_strand != get_strand(read) \
        and max(sa_start, read.reference_start) <= min(sa_end, read.reference_end)


def is_accidental_2d(read):
    """
    Return True for supplementary alignments for which a single other (the primary) alignment exists
    for which the strand is reverse and the reference span is overlapping
    """
    if read.has_tag('SA') and read.is_supplementary:
        supps = [get_sa_attributes(s) for s in read.get_tag('SA').split(';') if s]
        if len(supps) > 1:
            return False
        start_sa, end_sa, strand_sa = supps[0]
        if reverse_overlaps(*supps[0], read):
            return True
        else:
            return False
    else:
        return False


def follow_up_2d_candidates(candidates, distance=500):
    df = pd.DataFrame.from_records(
        data=[(r.query_name, r.reference_name, r.reference_start, r.reference_end)
              for r in candidates],
        columns=["RNAME", "Chromosome", "Start", "End"])
    keep_these_reads = []
    for chrom in df["Chromosome"].unique():
        keep_these_reads.extend(
            follow_up_2d_candidates_per_chrom(
                df[df["Chromosome"] == chrom].sort_values(by=["Start", "End"]),
                distance=distance)
        )
    return [r for r in candidates if r.query_name in keep_these_reads]


def follow_up_2d_candidates_per_chrom(df, distance=500):
    """
    Per chromosome, calculate the minimum of the distance with the previous and next read
    Distances are calculated based on read start position.
    Return only those read identifiers that are <= <distance> from each other.
    """
    df["distance"] = pd.concat([df["Start"].diff().abs(),
                                df["Start"].diff(periods=-1).abs()],
                               axis="columns").min(axis='columns')
    return df.loc[df["distance"] < distance, "RNAME"].to_list()
