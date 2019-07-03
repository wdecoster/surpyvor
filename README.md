# surpyvor
A python wrapper around [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), with additional convenience functions.

## Installation and dependencies
surpyvor requires bcftools, bgzip, tabix and SURVIVOR to be installed and in the $PATH.
Required python modules are cyvcf2, matplotlib, numpy, matplotlib-venn and upsetplot



## sub-commands:
    merge               merging vcf files of SVs
    highsens            get union of SV vcfs
    highconf            get intersection of SV vcfs
    prf                 calculate precision, recall and F-measure
    upset               Make upset plot for multiple SV vcf files
    venn                Make venn diagram for multiple SV vcf files

Each sub-command has its own help information, accessible by running `surpyvor <command> -h/--help`

## Citation
If you use this tool, please consider citing our [publication](https://genome.cshlp.org/content/early/2019/06/11/gr.244939.118.abstract) and the [citation for SURVIVOR](https://www.nature.com/articles/ncomms14061).
