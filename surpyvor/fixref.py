from cyvcf2 import VCF, Writer
from pyfaidx import Fasta


def fixref(vcf, fasta):
    """
    Fix reference alleles in VCF file.

    Parameters
    ----------
    vcf : str
        Path to VCF file.
    fasta : str
        Path to FASTA file.
    """
    fas = Fasta(fasta)
    vcf = VCF(vcf)
    w = Writer("-", vcf)
    for v in vcf:
        v.REF = fas[v.CHROM][v.start : v.end].seq
        w.write_record(v)
