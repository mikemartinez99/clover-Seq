#!/bin/bash

#SBATCH --job-name=cloverSeq
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16  
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=cloverSeq_%j.out

#----- START
echo "#------------------------ Initialization ------------------------#"
echo "Starting job: $SLURM_JOB_NAME (Job ID: $SLURM_JOB_ID)"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo -e "#-------------------------------------------------------------#"


#----- Source conda and activate snakemake
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

#----- Run snakemake workflow
snakemake -s Snakefile \
    --use-conda \
    --conda-frontend conda \
    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/GDSC-Clover-Seq \
    --profile cluster_profile \
    --rerun-incomplete \
    --keep-going 
    
#----- END
echo "End time: $(date)"

#----- Run Testing workflow
#snakemake -s database_Snakefile.smk \
#    --configfile database_configs/test_new.yaml \
#    --use-conda \
#    --conda-frontend conda \
#    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline \
#    --profile cluster_profile \
#    --rerun-incomplete \
#    --keep-going 

