import os
import tempfile
from cyvcf2 import VCF
import subprocess
import shlex
import pandas as pd
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


def confusion_matrix(vcff, names):
    """
    First level of the dict is the "first" call, second level is the "second" sample
    0: hom_ref
    1: heterozygous
    2: unknown/nocall
    3: hom_alt
    """
    zygosities = {0: {0: 0, 1: 0, 2: 0, 3: 0},
                  1: {0: 0, 1: 0, 2: 0, 3: 0},
                  2: {0: 0, 1: 0, 2: 0, 3: 0},
                  3: {0: 0, 1: 0, 2: 0, 3: 0},
                  }
    for v in VCF(vcff):
        zygosities[v.gt_types[0]][v.gt_types[1]] += 1
    zygs = [2, 0, 1, 3]
    df = pd.DataFrame(index=zygs, columns=zygs)
    for tr in zygs:
        for te in zygs:
            df.loc[tr, te] = zygosities[tr][te]
    df.columns = ['nocall', 'hom_ref', 'het', 'hom_alt']
    df.columns.name = names[1]
    df.index = ['nocall', 'hom_ref', 'het', 'hom_alt']
    df.index.name = names[0]
    print(df)
