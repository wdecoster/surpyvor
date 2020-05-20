import subprocess
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna
import pysam
from cyvcf2 import VCF, Writer
import shlex
from cigar import Cigar


def process(vcff, bamf, output):
    vcf = VCF(vcff)
    output = Writer(output, vcf)
    bam = pysam.AlignmentFile(bamf, "rb")
    for v in vcf:
        if v.INFO.get('SVTYPE') == 'INS':
            seq = improve_insertion_sequence(bam, v.CHROM, v.start, v.end, v.INFO.get('SVLEN'))
            if seq:
                v.ALT = seq
        output.write_record(v)


def improve_insertion_sequence(bam, chrom, start, end, svlen, window=250):
    sequences = [get_inserted_sequence(read, minpos=start-window, maxpos=end+window)
                 for read in bam.fetch(contig=chrom, start=start-window, end=end+window)
                 if not read.is_secondary]
    return assemble_and_consensus([s for s in sequences if svlen*0.75 < len(s) < svlen*1.25])


cigar_dict = {i: o for i, o in zip(range(0, 10), 'MIDNSHP=XB')}


def get_inserted_sequence(read, minpos, maxpos, min_ins_length=25):
    '''
    Get the inserted sequences in a read, requiring that the CIGAR operation is at least
    of length min_ins_length and the CIGAR operations reference coordinates
    are between minpos and maxpos
    '''
    position = 0
    seq = []
    cigar_list = []
    for operation, length in read.cigartuples:
        cigar_list.append(f"{length}{cigar_dict[operation]}")
        if operation == 1 and length > min_ins_length \
                and minpos <= get_ref_pos(read.reference_start, cigar_list) <= maxpos:
            seq.append(read.query_sequence[position:position + length])
        if operation not in [2, 3, 5, 6]:  # these operations do not consume the query
            position += length
    return reverse_complement_table(''.join(seq)) if read.is_reverse else ''.join(seq)


def get_ref_pos(refstart, cigar_list):
    '''
    Get the reference position of a position based
    on the start of the read and a (partial) CIGAR string
    '''
    return refstart + Cigar(''.join(cigar_list)).reference_length()


tab = str.maketrans("ACTG", "TGAC")


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


def reverse_complement(seq):
    return


def assemble_and_consensus(sequences):
    child = subprocess.Popen(shlex.split("muscle -clwstrict -maxiters 2"),
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    SeqIO.write(
        sequences=[SeqRecord(Seq(s, generic_dna), id=f"ins{i}") for i, s in enumerate(sequences)],
        handle=child.stdin,
        format="fasta")
    child.stdin.close()
    align = AlignIO.read(child.stdout, "clustal")
    summary_align = AlignInfo.SummaryInfo(align)
    consensus = summary_align.dumb_consensus(ambiguous='N', require_multiple=True)
