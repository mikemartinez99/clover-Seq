#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# clover-Seq pipeline for analysis of tRNAs
# 
# Authors: Mike Martinez
# Lab: GDSC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- TO DO

#- Add code and conda environment to make some cool plots at the read counts step
#- Add code to make some QC plots outside of multiqc
#- Add prebuilt configs for different host species
#- Built tRNA-only bowtie2 indices
#- Figure out the mitoDB issue
#- Build tRNA + mitoRNA db if we can figure out the mito issue
#- Break down the coveragePlot source code into new code
#- Conda environment for featureCounts

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import pandas as pd 
import csv

#----- Set config file
configfile: "config.yaml"

#----- Read in the sample data
samples_df = pd.read_table(config["sample_txt"], delimiter = ",").set_index("Sample_ID", drop = False)
sample_list = list(samples_df["Sample_ID"])
genome = config["genome"]

#----- Generate run script for read counting
def generate_runfile(sample_file):
    with open(sample_file, 'r') as infile, open("runfile.txt", 'w') as outfile:
        reader = csv.DictReader(infile)
        for row in reader:
            sample_id = row["Sample_ID"]
            group = row['Group']
            outfile.write(f"{sample_id} {group} tRNA_alignment\n")

#----- Run function
generate_runfile(config["sample_txt"])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Final Rule
rule all:
    input:
        #----- Rule trimming outputs
        expand("trimming/{sample}.R1.trim.fastq.gz", sample = sample_list),
        expand("trimming/logs/{sample}.cutadapt.report", sample = sample_list),

        #----- Rule tRNA_alignment outputs
        expand("tRNA_alignment/{sample}.alignment.log.txt", sample = sample_list),
        expand("tRNA_alignment/{sample}.srt.bam", sample = sample_list),
        expand("tRNA_unaligned/{sample}.unalign.fastq", sample = sample_list),

        #----- Rule tRNA_mark_duplicates outputs
        expand("tRNA_alignment/{sample}.mkdup.bam", sample = sample_list),
        expand("tRNA_alignment/{sample}.mkdup.log.txt", sample = sample_list),

        #----- Rule tRNA_map_stats outputs
        expand("tRNA_alignment_stats/{sample}.mkdup.bam.idxstats", sample = sample_list),
        expand("tRNA_alignment_stats/{sample}.mkdup.bam.flagstat", sample = sample_list),

        #----- Rule tRNA_count outputs
        "tRNA_counts/genetype_counts.txt",
        "tRNA_counts/tRNA_isotype_counts.txt",
        "tRNA_counts/gene_level_counts_detailed.txt",
        "tRNA_counts/gene_level_counts_collapsed.txt",
        "tRNA_counts/tRNA_ends_counts.txt",

        #----- Rule read_length_distribution outputs
        "tRNA_alignment/full_alignment_read_length_distribution.txt",

        #----- Rule count_smRNAs outputs
        "smRNA_counts/raw_amino_counts_by_group.txt",
        "smRNA_counts/read_length_distribution.txt",
        "smRNA_counts/smRNA_raw_counts_by_group.txt",
        "smRNA_counts/smRNA_raw_counts_by_sample.txt",

        #----- Rule exploratory_data_analysis outputs
        #"plots/PCA_Plot.png",
        #"plots/read_length_distribution_by_sample.png",
        
    output:
        "QC/tRNA_multiqc_report.html"
    conda: "r_viz"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        genome = config["genome"],
        vis_script = config["vis_script"]
    shell:"""
    
        #----- Run MultiQC Report
        multiqc \
            trimming/logs \
            tRNA_alignment \
            tRNA_alignment_stats \
            tRNA_counts \
            -n QC/tRNA_multiqc_report.html \
            -c multiqc_config.yaml
    
    """

#----- Rule to trim
rule trimming:
    output:
        trim_1 = "trimming/{sample}.R1.trim.fastq.gz",
        report = "trimming/logs/{sample}.cutadapt.report"
    conda: "clover-seq"
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
rule tRNA_align:
    input:
        trim_1 = "trimming/{sample}.R1.trim.fastq.gz",
    output:
        #sam = "alignment/{sample}.aln.sam",
        alignLog = "tRNA_alignment/{sample}.alignment.log.txt",
        unalign = "tRNA_unaligned/{sample}.unalign.fastq",
        srtBam = "tRNA_alignment/{sample}.srt.bam"
    conda: "clover-seq"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        bt2_index = config["bt2_index"],
        maxMaps = config["maxMaps"],
        nPenalty = config["nPenalty"],
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
            -S tRNA_alignment/{params.sample}.aln.sam 2> {output.alignLog}

        #----- subset reads for aligned length > 15 & < 90bp & any reads with gaps (XO/XG tags)
        samtools view -h tRNA_alignment/{params.sample}.aln.sam | \
            awk 'BEGIN {{OFS="\t"}} $1 ~ /^@/ || ((length($10) > 15 && length($10) <= 90))' | \
            samtools view -Sb - > tRNA_alignment/{params.sample}.bam

        #----- filter for any reads with MAPQ <=1
        samtools view -h -q 2 tRNA_alignment/{params.sample}.bam > tRNA_alignment/{params.sample}.sub.bam

        #----- Sort and filter the bam file
        samtools sort -@ 4 tRNA_alignment/{params.sample}.sub.bam > {output.srtBam}
        samtools index {output.srtBam}

        #----- Remove temp files
        rm -rf tRNA_alignment/{params.sample}.aln.sam
        rm -rf tRNA_alignment/{params.sample}.bam
        rm -rf tRNA_alignment/{params.sample}.sub.bam
    """

#----- Rule to mark duplicates
rule tRNA_mark_duplicates:
    input:
        bam = "tRNA_alignment/{sample}.srt.bam"
    output:
        mkdup = "tRNA_alignment/{sample}.mkdup.bam",
        mkdupLog = "tRNA_alignment/{sample}.mkdup.log.txt"
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
        mkdup = "tRNA_alignment/{sample}.mkdup.bam"
    output:
        idxStats = "tRNA_alignment_stats/{sample}.mkdup.bam.idxstats",
        flagStats = "tRNA_alignment_stats/{sample}.mkdup.bam.flagstat"
    conda: "clover-seq"
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
        genetypeFile = "tRNA_counts/genetype_counts.txt",
        tRNA_isotype_counts = "tRNA_counts/tRNA_isotype_counts.txt",
        trnaCountsDetailed = "tRNA_counts/gene_level_counts_detailed.txt",
        trnaCountsCollapsed = "tRNA_counts/gene_level_counts_collapsed.txt",
        trnaEnds = "tRNA_counts/tRNA_ends_counts.txt"
    conda: "clover-seq"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    params:
        countScript = "code/countreads.py",
        runFile = config["runFile"],
        trna_db = config["trna_db"]
    shell: """
    
        #----- Run the countreads.py script
        python {params.countScript} \
            --samplefile={params.runFile} \
            --ensemblgtf={params.trna_db}/genes.gtf \
            --trnaloci={params.trna_db}/db-trnaloci.bed \
            --maturetrnas={params.trna_db}/db-maturetRNAs.bed \
            --trnatable={params.trna_db}/db-trnatable.txt \
            --genetypefile={output.genetypeFile} \
            --trnaends={output.trnaEnds} \
            --trnacounts={output.tRNA_isotype_counts} > {output.trnaCountsDetailed}

        #----- Collapse the detailed tRNA counts to overall tRNA counts
        awk '
        NR==1 {{ 
            print; 
            next 
        }}
        {{
            split($1, arr, "_");
            base = arr[1];
            if (!(base in seen)) {{
                order[++count] = base;
                seen[base] = 1;
            }}
            for (i=2; i<=NF; i++) {{
                counts[base, i] += $i;
            }}
        }}
        END {{
            for (i=1; i<=count; i++) {{
                tRNA = order[i];
                printf "%s", tRNA;
                for (j=2; j<=NF; j++) {{
                    printf "\t%d", counts[tRNA,j] + 0;
                }}
                print "";
            }}
        }}
        ' {output.trnaCountsDetailed} > {output.trnaCountsCollapsed}   
    """

#----- Rule to plot read-length distributions for ALL reads (not just tRNAs)
rule read_length_distribution:
    input:
        expand("tRNA_alignment/{sample}.mkdup.bam", sample = sample_list)
    output:
        distribution = "tRNA_alignment/full_alignment_read_length_distribution.txt"
    conda: "clover-seq"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
    shell: """
    
        #----- Calculate read length distribution information
        echo -e "Length\tSample\tCount" > {output.distribution}

        # Loop through all BAM files
        for bam in {input}; do
            # Extract the sample name from the filename
            sample=$(basename "$bam" .mkdup.bam)

            # Generate read length distribution using samtools and awk
            samtools view "$bam" | \
                awk -v sample="$sample" '{{ 
                    len = length($10)
                    counts[len]++
                }} 
                END {{ 
                    for (i = 0; i <= 100; i++) {{
                        printf "%d\t%s\t%d\\n", i, sample, counts[i]+0
                    }}
                }}' >> {output.distribution}
        done
    """

#----- Rule to count other smRNAs
rule count_smRNAs:
    input: 
        expand("tRNA_alignment/{sample}.mkdup.bam", sample = sample_list),
    output:
        aminoCounts = "smRNA_counts/raw_amino_counts_by_group.txt",
        readLengths = "smRNA_counts/read_length_distribution.txt",
        groupCounts = "smRNA_counts/smRNA_raw_counts_by_group.txt",
        counts = "smRNA_counts/smRNA_raw_counts_by_sample.txt"
    conda: "clover-seq"
    resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
    params:
        smRNA_count = "code/count_all_smRNA.py",
        runFile = config["runFile"],
        trna_db = config["trna_db"]

    shell: """
    
        #----- Run the code to count all tRNA + smRNA
        python {params.smRNA_count} \
            --samplefile={params.runFile} \
            --trnatable={params.trna_db}/db-trnatable.txt \
            --ensemblgtf={params.trna_db}/genes.gtf \
            --trnaloci={params.trna_db}/db-trnaloci.bed \
            --maturetrnas={params.trna_db}/db-maturetRNAs.bed \
            --trnaaminofile={output.aminoCounts} \
            --readlengthfile={output.readLengths} \
            --realcountfile={output.counts} \
            --countfile={output.groupCounts}
    
    """

#----- Rule to output QC plots
#rule exploratory_analysis:
 #   input:
 #       counts = "tRNA_counts/tRNA.readcounts.ann.tsv",
 #       lengths = "tRNA_alignment/read_length_distribution.txt"
 #   output:
  #      "plots/PCA_Plot.png",
  #      "plots/read_length_distribution_by_sample.png"
  #  conda: "r_viz"
  #  resources: cpus="12", maxtime="6:00:00", mem_mb="60gb"
  #  params:
  #      vis_script = "code/visualizations/CloverSeq_Exploratory_Utils.R",
   #     read_length_script = "code/visualizations/plot_readlength_distribution.R"
   # shell: """
    
        #----- Plot (Arg1 = tRNA annotated counts, Arg2 = output dir)
   #     Rscript {params.vis_script} \
   #         {input.counts} \
   #         plots/

        #----- Plot read length distributions
    #    Rscript {params.read_length_script} \
    #        {input.lengths} \
     #       plots/

    
  #  """

