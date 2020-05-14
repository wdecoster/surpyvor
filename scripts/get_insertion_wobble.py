import pysam
from cyvcf2 import VCF
import sys
from cigar import Cigar
import matplotlib.pyplot as plt


def main():
    process(sys.argv[1], sys.argv[2])


def process(vcff, bamf):
    vcf = VCF(vcff)
    bam = pysam.AlignmentFile(bamf, "rb")
    distances = []
    for v in vcf:
        if v.INFO.get('SVTYPE') == 'INS':
            distances.extend(get_distances(v.CHROM, (v.start+v.end)/2, bam))
    plt.hist([i for i in distances if not abs(i) > 5000])
    plt.xlabel('Distance from call')
    plt.xlim(-5000, 5000)
    plt.title("Distance between call and individual insertions")
    plt.savefig('insertion_distances.png')


cigar_dict = {i: o for i, o in zip(range(0, 10), 'MIDNSHP=XB')}


def get_distances(chrom, position, bam, window=1000):
    distances = []
    for read in bam.fetch(contig=chrom, start=position-window, end=position+window):
        cigar_list = []
        for operation, length in read.cigartuples:
            if (operation == 1 and length > 25):
                distances.append(position - (read.reference_start +
                                             Cigar(''.join(cigar_list)).reference_length()))
                break
            else:
                cigar_list.append(f"{length}{cigar_dict[operation]}")
    return distances


if __name__ == '__main__':
    main()
