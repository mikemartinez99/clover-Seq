# clover-Seq
Pipeline for analyses of tRNA-derived small RNAs (tDRs), mature tRNAs, and inference of RNA modification from high-throughput sequencing. 

# Table of Contents
- [Introduction](#introduction)
- [Databases](#databases)
- [Development Notes](#development_notes)
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

## Development Notes:
-  Is mapping too stringent (20% of bases are modified in tRNA)
AMino counts not informative
Coverage plots are useful (normalized read count)
Mutations: Linear modeling and bundance (this will always be pairwise) - miRNA-specific data for mining
Estebans data are size selected (tRNA selecting anything under 200 nt in length (snRNAs, snoRNAs, tRNAs, miRNAs, fragments))
CCA, CC, C endpoint tails : shows maturity (pretRNA is usually longer than the mature tRNA)
tRNA modifications heatmap

VOlcano plot for just tRNAs and then a separate one for all other small RNAs
140 - 250 nts returned (120 of those are adapters)
Also might be interested in mapping tRNA fragments

Clean up microbes step....
Might also be working in zebrafish

#----- Jargon...
Isodecoders - isoforms of tRNAs (gtRNAdb) --RSEM????? -> as long as you have a gtf
Isoacceptros - grouping by family (codon rdegeneracy)

Seed length too long


## Helpful commands:

To run getcoverage.py from the command line (within code folder)--- Ensure you have trax_env activated
```shell
python getcoverage.py --bedfile=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.bed --samplefile=MM_Working_Scripts/sample_file.txt --stkfile=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-trnaalign.stk --allcoverage=all-coverage.txt --trnafasta=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-maturetRNAs.fa --sizefactors=test/test-SizeFactors.txt --locistk=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db-trnaloci.stk 
```

To run processsamples.py from the command line
```shell
python 03_processsamples.py --experimentname=test --databasename=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/db --ensemblgtf=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/tRAX_databases/hg38_db/genes.gtf --samplefile=MM_Working_Scripts/sample_file.txt --bamdir=/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/shared-software/workflows/clover-Seq/tRNA_alignment/ --nofrag --exppairs=MM_Working_Scripts/pairfile.txt

```

