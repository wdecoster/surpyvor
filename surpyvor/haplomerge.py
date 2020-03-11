import tempfile
from cyvcf2 import VCF
from surpyvor import utils


def merge_split_called_haplotypes(merged, output):
    _, name = tempfile.mkstemp(suffix='.vcf')
    vcf = VCF(merged)

    with open(name, 'a') as tmpoutput:
        tmpoutput.write("{}\n".format('\n'.join(make_header(vcf))))
        for v in vcf:
            info = {'SVLEN': v.INFO.get('SVLEN'), 'END': v.end, 'SVTYPE': v.INFO.get('SVTYPE')}
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


def get_genotype(alleles):
    hap1, hap2 = [int(utils.is_variant(c)) for c in alleles]
    return '{}|{}'.format(hap1, hap2)


def make_header(vcf):
    header = ['##fileformat=VCFv4.1', '##source=surpyvor haplomerge']
    for line in vcf.header_iter():
        if line["HeaderType"] == 'CONTIG':
            header.append('##contig=<ID={}>'.format(line['ID']))
    header.extend([
        '##ALT=<ID=DEL,Description="Deletion">',
        '##ALT=<ID=DUP,Description="Duplication">',
        '##ALT=<ID=INV,Description="Inversion">',
        '##ALT=<ID=BND,Description="Translocation">',
        '##ALT=<ID=INS,Description="Insertion">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the structural variant">',
        '##INFO=<ID=SVLEN,Number=1,Type=Float,Description="Length of the SV">',
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of the SV.">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE'])
    return header
