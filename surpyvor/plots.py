from cyvcf2 import VCF
import matplotlib.pyplot as plt
from surpyvor import utils
import numpy as np
from matplotlib_venn import venn2, venn3
from upsetplot import plot as upsetplot


def bar_chart(vcf, outname="stacked_bar.png"):
    """
    Make a stacked bar chart for length of the SV split by validation status
    This ignores zygosity.
    """
    len_dict = {"True": [], "False": [], "Missed": []}
    for v in VCF(vcf):
        if not v.INFO.get('SVTYPE') == 'TRA' and abs(v.INFO.get('SVLEN')) >= 50:
            calls = [utils.is_variant(call) for call in v.gt_types]
            if calls == [True, True]:
                len_dict['True'].append(v.INFO.get('SVLEN'))
            elif calls == [False, True]:
                len_dict['False'].append(v.INFO.get('SVLEN'))
            elif calls == [True, False]:
                len_dict['Missed'].append(v.INFO.get('SVLEN'))
    plt.subplot(2, 1, 1)
    plt.hist(x=np.array(list(len_dict.values())),
             bins=[i for i in range(0, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=list(len_dict.keys()))
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.subplot(2, 1, 2)
    plt.hist(x=np.array(list(len_dict.values())),
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=list(len_dict.keys()),
             log=True)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(outname)
    plt.close()


def upset_plot(upsets, outname="UpSetPlot.png"):
    upsetplot(upsets, sort_by='cardinality')
    plt.savefig(outname)


def venn_diagram(sets, labels, num_samples=2, outname="venn.png"):
    if num_samples == 2:
        venn = venn2
    else:
        venn = venn3
    venn(sets, set_labels=labels)
    plt.savefig(outname)
    plt.close()
