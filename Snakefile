#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GDSC-tRAX v2 Pipeline (Claude Version)
# 
# Authors: Mike Martinez
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- TO DO
#- Add columns in read counts to denote isoacceptor (i.e., the AA)
#- Add code and conda environment to make some cool plots at the read counts step
#- Add code to make some QC plots outside of multiqc

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import pandas as pd 

#----- Set config file
configfile: "config.yaml"

#----- Read in the sample data
samples_df = pd.read_table(config["sample_txt"], delimiter = ",").set_index("Sample_ID", drop = False)
sample_list = list(samples_df["Sample_ID"])
genome = config["genome"]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Final Rule
rule all:
    input:
        #----- Rule trimming outputs
        expand("trimming/{sample}.R1.trim.fastq.gz", sample = sample_list),
        expand("trimming/logs/{sample}.cutadapt.report", sample = sample_list),

        #----- Rule alignment outputs
        expand("alignment/{sample}.alignment.log.txt", sample = sample_list),
        expand("alignment/{sample}.srt.bam", sample = sample_list),
        expand("unaligned/{sample}.unalign.fastq", sample = sample_list),

        #----- Rule mark_duplicates outputs
        expand("alignment/{sample}.mkdup.bam", sample = sample_list),
        expand("alignment/{sample}.mkdup.log.txt", sample = sample_list),

        #----- Rule tRNA_map_stats outputs
        expand("tRNA_alignment_stats/{sample}.mkdup.bam.idxstats", sample = sample_list),
        expand("tRNA_alignment_stats/{sample}.mkdup.bam.flagstat", sample = sample_list),

        #----- Rule tRNA_count outputs
        expand("tRNA_counts/tRNA.readcounts.tsv", sample = sample_list),
        expand("tRNA_counts/tRNA.readcounts_tpm.tsv", sample = sample_list),
        expand("tRNA_counts/tRNA.readcounts.ann.tsv", sample = sample_list),
        expand("tRNA_counts/tRNA.readcounts_tpm.ann.tsv", sample = sample_list),
    output:
        "done.txt"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        genome = config["genome"]
    shell:"""
    
        #----- Make dummy file
        touch done.txt
    
    """

#----- Rule to trim
rule trimming:
    output:
        trim_1 = "trimming/{sample}.R1.trim.fastq.gz",
        report = "trimming/logs/{sample}.cutadapt.report"
    conda: "cutadapt"
    resources: cpus="8", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        fastq_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        adapter_1 = config["adapter_1"],
        minlength = config["minlength"]
    shell: """
    
        #----- Run cutadapt
        cutadapt \
            -o {output.trim_1} \
            {params.fastq_1} \
            -m {params.minlength} \
            -a {params.adapter_1} \
            -j {resources.cpus} > {output.report} 2>&1
    """

#----- Rule to align samples to tRNA database (need db-trnatable.txt)
rule align:
    input:
        trim_1 = "trimming/{sample}.R1.trim.fastq.gz",
    output:
        #sam = "alignment/{sample}.aln.sam",
        alignLog = "alignment/{sample}.alignment.log.txt",
        unalign = "unaligned/{sample}.unalign.fastq",
        srtBam = "alignment/{sample}.srt.bam"
    conda: "trax_env"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        database = config["trna_db"],
        bt2_index = config["bt2_index"],
        maxMaps = config["maxMaps"],
        nPenalty = config["nPenalty"],
        minnontrnasize = config["minlength_nontRNA"],
        TMPDIR = "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/shared-software/workflows/tRAX_v2/temp"
    shell: """
    
        #----- Run Bowtie2 alignment (removed the --reorder arg from original code to speed up...don't see the need with SE data)
        bowtie2 \
            -x {params.bt2_index} \
            -U {input.trim_1} \
            -D 20 -R 3 -N 1 -L 12 -i S,1,0.50 \
            -k {params.maxMaps} \
            --very-sensitive \
            --np {params.nPenalty} \
            --reorder \
            --ignore-qual \
            --un {output.unalign} \
            -p {resources.cpus} \
            -S alignment/{params.sample}.aln.sam 2> {output.alignLog}

        #----- subset reads for aligned length > 16 & < 28bp & any reads with gaps (XO/XG tags)
        samtools view -h alignment/{params.sample}.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 15 && length($10) <= 90) && ($0 !~ /XG:i:[^0]/ && $0 !~ /XO:i:[^0]/)) {{print $0}}' | \
            samtools view -Sb - > alignment/{params.sample}.bam

        #----- filter for any reads with MAPQ <=1
        samtools view -h -q 2 alignment/{params.sample}.bam > alignment/{params.sample}.sub.bam

        #----- Sort and filter the bam file
        samtools sort -@ 4 alignment/{params.sample}.sub.bam > {output.srtBam}
        samtools index {output.srtBam}

        #----- Remove temp files
        rm -rf alignment/{params.sample}.aln.sam
        rm -rf alignment/{params.sample}.bam
        rm -rf alignment/{params.sample}.sub.bam
    """

#----- Rule to mark duplicates
rule mark_duplicates:
    input:
        bam = "alignment/{sample}.srt.bam"
    output:
        mkdup = "alignment/{sample}.mkdup.bam",
        mkdupLog = "alignment/{sample}.mkdup.log.txt"
    conda: "rnaseq1"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample
    shell: """
    
        #----- Run Picard mark Duplicates
        picard -Xmx16G -Xms16G  \
            MarkDuplicates \
            I={input.bam} \
            O={output.mkdup} \
            M={output.mkdupLog} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
            CREATE_INDEX=false  \
            MAX_RECORDS_IN_RAM=4000000 \
            ASSUME_SORTED=true \
            MAX_FILE_HANDLES=768

        #----- Index the mkdup bam
        samtools index {output.mkdup}
    """

#----- Rule to collate tRNA mapping statistics
rule tRNA_map_stats:
    input:
        mkdup = "alignment/{sample}.mkdup.bam"
    output:
        idxStats = "tRNA_alignment_stats/{sample}.mkdup.bam.idxstats",
        flagStats = "tRNA_alignment_stats/{sample}.mkdup.bam.flagstat"
    conda: "trax_env"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params: 
        sample = lambda wildcards: wildcards.sample
    shell: """
    
        #----- Collect metrics
        samtools idxstats {input.mkdup} > {output.idxStats}
        samtools flagstat {input.mkdup} > {output.flagStats}
    
    """

#----- Rule to count tRNAs
rule tRNA_count:
    input:
        expand("tRNA_alignment_stats/{sample}.mkdup.bam.idxstats", sample = sample_list)
    output:
        tRNA_counts = "tRNA_counts/tRNA.readcounts.tsv",
        tRNA_tpms = "tRNA_counts/tRNA.readcounts_tpm.tsv",
        isoAcc_anno = "tRNA_counts/tRNA.readcounts.ann.tsv",
        isoAcc_tpms_anno = "tRNA_counts/tRNA.readcounts_tpm.ann.tsv"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    params:
        countScript = config["countScript"],
        tpm_normalization = config["tpm_normalization"]
    shell: """
    
        #----- Count tRNAs
        echo -ne tRNA_ID"\t"Length"\t" > {output.tRNA_counts} 
        echo {input} | tr " " "\t"| sed s/"tRNA_alignment_stats\/"//g| sed s/".mkdup.bam.idxstats"//g >> {output.tRNA_counts}
        paste {input}| awk -f {params.countScript} >> {output.tRNA_counts}

        #----- Run TPM normalization
        python {params.tpm_normalization} {output.tRNA_counts}

        #----- Add isoacceptor information
        awk 'BEGIN {{OFS="\t"}} NR==1 {{print $0, "Isoacceptor"; next}} {{split($1, a, "-"); print $0, a[2]}}' {output.tRNA_counts} > {output.isoAcc_anno}
        awk 'BEGIN {{OFS="\t"}} NR==1 {{print $0, "Isoacceptor"; next}} {{split($1, a, "-"); print $0, a[2]}}' {output.tRNA_tpms} > {output.isoAcc_tpms_anno}

    """
