import matplotlib.pyplot as plt


def bar_chart(vcf, outname="stacked_bar.png"):
    """
    Make a stacked bar chart for length of the SV split by validation status
    This ignores zygosity.
    """
    from cyvcf2 import VCF
    from surpyvor import utils
    import numpy as np

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
    from upsetplot import plot as upsetplot
    upsetplot(upsets, sort_by='cardinality')
    plt.savefig(outname)


def venn_diagram(sets, labels, num_samples=2, outname="venn.png"):
    from matplotlib_venn import venn2, venn3
    if num_samples == 2:
        from matplotlib_venn import venn2 as venn
    else:
        from matplotlib_venn import venn3 as venn
    venn(sets, set_labels=labels)
    plt.savefig(outname)
    plt.close()


def length_plot(dict_of_lengths, output):
    """Makes two stacked bar charts
    Plotting two bar charts of number of SVs by length split by SV type
    Use a consistent colouring scheme for those in "standard_order" to
    make comparison reasonable

    First bar chart is up to 2kb with bins of 10bp
    Second bar chart is up to 20kb, with bins of 100bp
     and uses log scaling on the y-axis
    """
    standard_order = ['DEL', 'INS', 'INV', 'DUP']
    spec_order = sorted([i for i in dict_of_lengths.keys() if i not in standard_order])
    sorter = standard_order + spec_order
    names, lengths = zip(
        *sorted([(svtype, lengths) for svtype, lengths in dict_of_lengths.items()],
                key=lambda x: sorter.index(x[0])))
    plt.subplot(2, 1, 1)
    plt.hist(x=lengths,
             bins=[i for i in range(50, 2000, 10)],
             stacked=True,
             histtype='bar',
             label=names)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")

    plt.subplot(2, 1, 2)
    plt.hist(x=lengths,
             bins=[i for i in range(0, 20000, 100)],
             stacked=True,
             histtype='bar',
             label=names,
             log=True)
    plt.xlabel('Length of structural variant')
    plt.ylabel('Number of variants')
    plt.legend(frameon=False,
               fontsize="small")
    plt.tight_layout()
    plt.savefig(output)
