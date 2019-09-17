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


def get_svlengths(vcf):
    len_dict = defaultdict(list)
    vcf = VCF(vcf)
    if len(vcf.samples) > 1:
        sys.stderr.write("\n\nWarning: this script does not support multiple samples in a vcf.\n")
        sys.stderr.write("Plotting and counting only for {}.".format(vcf.samples[0]))
    for v in vcf:
        if is_variant(v.gt_types[0]) and not v.INFO.get('SVTYPE') == 'TRA':
            try:
                if get_svlen(v) >= 50:
                    len_dict[get_svtype(v)].append(get_svlen(v))
            except TypeError:
                if v.INFO.get('SVTYPE') == 'INV':
                    if (v.end - v.start) >= 50:
                        len_dict[get_svtype(v)].append(v.end - v.start)
                elif v.INFO.get('SVTYPE') == 'BND':
                    try:
                        len_dict['BND'].append(abs(v.INFO.get('MATEDIST')))
                    except TypeError:
                        len_dict['BND'].append(0)
                else:
                    sys.stderr.write("Exception when parsing variant:\n{}\n\n".format(v))
                    len_dict["parse_error"].append(0)
    return len_dict


def get_svtype(v):
    if v.INFO.get('SVTYPE') == "INVDUP":
        return "INV"
    else:
        return v.INFO.get('SVTYPE').split(':')[0].split('/')[0]


def get_svlen(v):
    return abs(v.INFO.get('SVLEN'))
