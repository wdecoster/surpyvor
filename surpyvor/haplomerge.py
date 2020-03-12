import tempfile
from cyvcf2 import VCF
from surpyvor import utils
import sys


def merge_split_called_haplotypes(merged, output, name=None):
    _, name = tempfile.mkstemp(suffix='.vcf')
    vcf = VCF(merged)
    if len(vcf.samples) == 3:
        get_genotype = get_genotype_from_three
    elif len(vcf.samples) == 2:
        get_genotype = get_genotype_from_two
    else:
        sys.exit("ERROR: Unexpected number of samples in haplomerge intermediate VCF!")
    with open(name, 'a') as tmpoutput:
        tmpoutput.write("{}\n".format('\n'.join(make_header(vcf, name=name))))
        allele_dict = {0: 'HOM_REF', 1: 'HET', 2: 'UNKNOWN', 3: 'HOM_ALT'}
        for v in vcf:
            info = {'SVLEN': v.INFO.get('SVLEN'),
                    'END': v.end,
                    'SVTYPE': v.INFO.get('SVTYPE'),
                    'HAPSUPPORT': '-'.join([allele_dict[gt] for gt in v.gt_types])}
            print("{chrom}\t{pos}\t{idf}\t{ref}\t{alt}\t{q}\t{filt}\t{info}\t{form}\t{sam}"
                  .format(
                      chrom=v.CHROM,
                      pos=v.start,
                      idf=v.ID,
                      ref=v.REF,
                      alt=','.join(v.ALT),
                      q=v.QUAL or '.',
                      filt='.',
                      info=';'.join(['{}={}'.format(k, v) for k, v in info.items()]),
                      form='GT',
                      sam=get_genotype(v.gt_types)),
                  file=tmpoutput)
    utils.vcf_sort(name, output)


def get_genotype_from_two(alleles):
    '''Combine calls on two haplotypes in one genotype'''
    hap1, hap2 = [int(utils.is_variant(c)) for c in alleles]
    return '{}|{}'.format(hap1, hap2)


def get_genotype_from_three(alleles):
    '''Combine calls on two haplotypes and unphased set in one genotype'''
    calls = [int(utils.is_variant(c)) for c in alleles]
    if sum(calls) >= 2 or calls == [0, 0, 1]:
        return '1|1'
    elif sum(calls) == 1:
        return '{}|{}'.format(calls[0], calls[1])
    else:  # I don't think should happen
        return '0|0'


def make_header(vcf, name=None):
    header = ['##fileformat=VCFv4.1', '##source=surpyvor haplomerge']
    for line in vcf.header_iter():
        if line["HeaderType"] == 'CONTIG':
            header.append('##contig=<ID={}>'.format(line['ID']))
    if name is None:
        name = "SAMPLE"
    header.extend([
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=BND,Description="Translocation">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">',
        '##INFO=<ID=HAPSUPPORT,Number=1,Type=String,Description="Support per haplotype.">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(name)])
    return header
