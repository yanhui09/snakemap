# Snakemake workflow: `snakemap`

[![Snakemake](https://img.shields.io/badge/snakemake-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/yanhui09/snakemap/actions/workflows/main.yml/badge.svg?branch=master)](https://github.com/yanhui09/snakemap/actions?query=branch%3Amaster+workflow%3ATests)


A Snakemake workflow for `snakemap`

A snakemake workflow for read-mapping against the contig catalogue from metagenomics by illumina and nanopore

# Usage

1. Init
   
```
python init.py -i /path/to/illumina_fqs -n /path/to/nanopore_fqs -w /path/to/workdir
```

2. Run

```
snakemake all -j6 --snakefile /path/to/Snakefile --configfile /path/to/config.yaml --directory /path/to/workdir --use-conda
```