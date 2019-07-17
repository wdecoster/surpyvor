# surpyvor
A python wrapper around [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), with additional convenience functions.

## Installation and dependencies
surpyvor requires bcftools, bgzip, tabix and SURVIVOR to be installed and in the $PATH.
Required python modules are cyvcf2, matplotlib, numpy, matplotlib-venn and upsetplot

surpyvor and its dependencies can be installed from [bioconda](https://anaconda.org/bioconda/surpyvor):

`conda install -c bioconda surpyvor`

## USAGE
### sub-commands:
    merge               merging vcf files of SVs
    highsens            get union of SV vcfs
    highconf            get intersection of SV vcfs
    prf                 calculate precision, recall and F-measure
    upset               Make upset plot for multiple SV vcf files
    venn                Make venn diagram for 2 or 3 SV vcf files

Each sub-command has its own help information, accessible by running `surpyvor <command> -h/--help`

### General and common arguments for most sub-commands:
```
-o/--output: output variant file to write. Default: stdout
--plotout: name ouf output plot to write. Default names depending on plot type.
-d/--distance: maximal pairwise distance between coordinates of SVs to be considered concordant. Default: 500
-l/--minlength: minimal SV length to include. Default: 50
--variants: vcf files to combine
```

### Specific arguments

#### surpyvor prf
```
--ignore_chroms: ignore some chromosomes for calculations. Default: chrEBV
--bar: create a stacked bar chart colored by validation status [not created by default]
--matrix: create a confusion matrix [not created by default]
```

## Citation
If you use this tool, please consider citing our [publication](https://genome.cshlp.org/content/early/2019/06/11/gr.244939.118.abstract) and the [citation for SURVIVOR](https://www.nature.com/articles/ncomms14061).
