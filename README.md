# clover-Seq
Pipeline for analyses of tRNA-derived small RNAs (tDRs), mature tRNAs, and inference of RNA modification from high-throughput sequencing. 

# Table of Contents
- [Introduction](#introduction)
- [Databases](#databases)
- [Summary](#summary)
- [Directories](#directories)
- [Files](#files)
- [Implementation](#implementation)
- [Contact](#contact)
- [Citation](#citation)

## Introduction
The pipeline was inspired by the [tRAX pipeline](https://github.com/UCSC-LoweLab/tRAX) but adapted to be compatible with Snakemake and allow parallel processing on the [Dartmouth Discovery HPC](https://rc.dartmouth.edu/discoveryhpc/). This pipeline supports the analysis of tDRs, mature-tRNAs and RNA modification for human (hg38), mouse (mm10), and fly(dm6) genomes. 

## Databases
tRNA-genome references encompass the full host genome with additional tRNA-specific gene information obtained from [gtRNAdb](https://gtrnadb.ucsc.edu) through [tRNA-scan](https://lowelab.ucsc.edu/tRNAscan-SE/) experiments. These references are pre-downloaded along with pre-built Bowtie2 indices and hosted on the [Genomics and Molecular Biology Shared Resources](https://geiselmed.dartmouth.edu/gsr/) on Discovery for ease of use and efficiency. However, if you wish to build from scratch, or customize the reference, a workflow for doing so is included. 

References can be accessed at the following path on Discovery: 

`/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases`


