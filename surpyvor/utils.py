import os
import sys
from shutil import which
import tempfile
from cyvcf2 import VCF
import subprocess
import shlex
import pandas as pd
import gzip
from collections import defaultdict


def is_variant(call):
    """Check if a variant position qualifies as a variant

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT"""
    if call == 1 or call == 3:
        return True
    else:
        return False


def normalize_vcf(vcff):
    """Normalize a vcf by changing DUP to INS"""
    handle, name = tempfile.mkstemp(suffix='.vcf')
    out = open(name, 'w')
    if vcff.endswith('.gz'):
        vcf = gzip.open(vcff, 'rt')
    else:
        vcf = open(vcff)
    for line in vcf:
        out.write(line.replace('DUP', 'INS'))
    os.close(handle)
    return name


def get_variant_identifiers(vcf, ignore_chroms, num_samples=2):
    """Get sets of variants for each sample in a merged vcf.

    Loop over the vcf file, adding a unique identifier to the
    respective list if the sample has a variant for that position
    return as set
    """
    positions = [[] for _ in range(num_samples)]
    for v in VCF(vcf):
        if v.CHROM not in ignore_chroms:
            for index, call in enumerate(v.gt_types):
                if is_variant(call):
                    positions[index].append(
                        "{}:{}-{}".format(v.CHROM, v.start, v.INFO.get('SVTYPE')))
    identifier_sets = [set(i) for i in positions]
    return identifier_sets


def gt_types_to_binary_comparison(calls):
    """From an array of calls, check if a variant position qualifies as a variant.

    0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT
    Return string of 1s and 0s to represent position"""
    binary_calls = []
    for call in calls:
        if call == 1 or call == 3:
            binary_calls.append(1)
        else:
            binary_calls.append(0)
    return ''.join([str(i) for i in binary_calls])


def make_sets(vcf, names):
    """From the merged SV file, return pd.Series of overlapping sets.

    Intended for making an upset plot"""
    calls = defaultdict(int)
    for v in VCF(vcf):
        calls[gt_types_to_binary_comparison(v.gt_types)] += 1
    tf_array = [[True, False]] * len(list(calls.keys())[0])
    index = pd.MultiIndex.from_product(tf_array, names=names)
    values = [calls[''.join([str(int(j)) for j in i])] for i in index]
    return pd.Series(values, index=index)


def vcf_concat(vcffiles):
    _, concatenated = tempfile.mkstemp(suffix=".vcf")
    sample = [get_sample(f) for f in vcffiles][0]
    vcffiles = [reheader(f, sample=sample) for f in vcffiles]
    vcffiles = [compress_and_tabix(f) for f in vcffiles]
    c = subprocess.Popen(shlex.split("bcftools concat -a {}".format(' '.join(vcffiles))),
                         stdout=subprocess.PIPE)
    subprocess.call(shlex.split("bcftools sort -o {}".format(concatenated)), stdin=c.stdout)
    return concatenated


def get_sample(vcffile):
    vcf = VCF(vcffile)
    return vcf.samples[0]


def reheader(vcf, sample):
    _, output = tempfile.mkstemp(suffix=".vcf")
    handle, samplef = tempfile.mkstemp()
    open(samplef, 'w').write(sample)
    os.close(handle)
    subprocess.call(shlex.split("bcftools reheader -s {} {} -o {}".format(samplef, vcf, output)))
    return output


def compress_and_tabix(vcf):
    if vcf.endswith('.vcf'):
        handle, output = tempfile.mkstemp(suffix=".vcf.gz")
        subprocess.call(shlex.split("bgzip -c {}".format(vcf)), stdout=handle)
        subprocess.call(shlex.split("tabix -p vcf {}".format(output)))
        return output
    else:
        return vcf


def decompress(vcf):
    """
    Decompress output to temporary file if filename endswith .gz or .bgz
    """
    if vcf.endswith(('.gz', '.bgz')):
        handle, output = tempfile.mkstemp(suffix=".vcf")
        subprocess.call(shlex.split("bgzip -cd {}".format(vcf)), stdout=handle)
        return output
    else:
        return vcf


def test_dependencies():
    for dependency in ['bcftools', 'bgzip', 'tabix', 'SURVIVOR']:
        if not which(dependency):
            sys.exit("ERROR: Could not find required executable '{}'.\n"
                     "Make sure it is installed and in $PATH".format(dependency))


def vcf_sort(input, output):
    if output in ["stdout", "-"]:
        subprocess.call(shlex.split('bcftools sort {}'.format(input)))
    else:
        subprocess.call(shlex.split('bcftools sort {} -o {}'.format(input, output)))


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


def merge_split_called_haplotypes(merged):
    vcf = VCF(merged)
    write_header(vcf)

    for v in vcf:
        info = {'SVLEN': v.INFO.get('AVGLEN'), 'END': v.end, 'SVTYPE': v.INFO.get('SVTYPE')}
        print("{chrom}\t{pos}\t{idf}\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t{form}\t{sam}".format(
            chrom=v.CHROM,
            pos=v.start,
            idf=v.ID,
            ref=v.REF,
            alt=','.join(v.ALT),
            qual=v.QUAL or '.',
            filt='.',
            info=';'.join(['{}={}'.format(k, v) for k, v in info.items()]),
            form='GT',
            sam=get_genotype(v.gt_types)
        ))


def get_genotype(alleles):
    alternative_alleles = sum([is_variant(c) for c in alleles])
    if alternative_alleles == 1:
        return '0/1'
    elif alternative_alleles == 2:
        return '1/1'
    else:
        return '0/0'


def write_header(vcf):
    print('##fileformat=VCFv4.1')
    print('##source=surpyvor haplomerge')
    for line in vcf.header_iter():
        if line["HeaderType"] == 'CONTIG':
            print('##contig=<ID={}>)'.format(line['ID']))
    print('##ALT=<ID=DEL,Description="Deletion">')
    print('##ALT=<ID=DUP,Description="Duplication">')
    print('##ALT=<ID=INV,Description="Inversion">')
    print('##ALT=<ID=BND,Description="Translocation">')
    print('##ALT=<ID=INS,Description="Insertion">')
    print('##INFO=<ID=END,Number=1,Type=Integer,Description="End of the structural variant">')
    print('##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">')
    print('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">')
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE')
