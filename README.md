# surpyvor
A python wrapper around [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR), with additional convenience functions.


optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version and quit.

## sub-commands:
    merge               merging vcf files of SVs
    highsens            get union of SV vcfs
    highconf            get intersection of SV vcfs
    prf                 calculate precision, recall and F-measure
    upset               Make upset plot for multiple SV vcf files

Each subcommand has its own help information, accessible by running `surpyvor <command> -h/--help`
