import os
import tempfile
from cyvcf2 import VCF
import subprocess
import shlex
import gzip


def is_variant(call):
    """Check if a variant position qualifies as a variant

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT"""
    if call == 1 or call == 3:
        return True
    else:
        return False


def normalize_vcf(vcff):
    handle, name = tempfile.mkstemp()
    out = open(name, 'w')
    if vcff.endswith('.gz'):
        vcf = gzip.open(vcff, 'rt')
    else:
        vcf = open(vcff)
    for line in vcf:
        out.write(line.replace('DUP', 'INS'))
    os.close(handle)
    return name


def get_variant_identifiers(vcf, ignore_chroms):
    positions = [[], []]
    for v in VCF(vcf):
        if v.CHROM not in ignore_chroms:
            for index, call in enumerate(v.gt_types):
                if is_variant(call):
                    positions[index].append(
                        "{}:{}-{}".format(v.CHROM, v.start, v.INFO.get('SVTYPE')))
    identifier_sets = [set(i) for i in positions]

    return identifier_sets


def vcf_concat(vcffiles):
    fhf, concatenated = tempfile.mkstemp()
    bcftools_concat_cmd = "bcftools concat {inp} | bcftools sort -o {out}".format(
        out=concatenated,
        inp=' '.join(vcffiles))
    subprocess.call(shlex.split(bcftools_concat_cmd))
    return concatenated
