#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: Full_tRNA_Genome_PCA.R
# Description: Generate PCA for full genome + tRNA counts
#
# Author: Mike Martinez
# Lab: GDSC
# Project: Clover-Seq
# Date created: 05/16/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(cowplot))
suppressMessages(library(DESeq2))
suppressMessages(library(stats))
suppressMessages(library(SummarizedExperiment))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE DATA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Arg 1: sample sheet
# Arg 2: reference level
# Arg 3: gene level counts collapsed
# Arg 4: tRNA isotype counts
# Arg 5: normaliezd data output file path
# Arg 6: PCA output path

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: RScript PCA.R <Sample_list_SE.txt> <reference level> <gene_level_counts_collapsed.txt> <tRNA_isotype_counts.txt> <normalized data output path/> <PCA output path/>")
}

metadata <- args[1]
refLevel <- args[2]
fullData <- args[3]
tRNAData <- args[4]
normalizedDir <- args[5]
pcaDir <- args[6]

if (!dir.exists(normalizedDir)) {
  dir.create(normalizedDir)
}

if (!dir.exists(pcaDir)) {
  dir.create(pcaDir)
}

#----- Debugging
message("--------------------------------------------------")
message(paste0("Metadata: ", metadata))
message(paste0("reference level: ", refLevel))
message(paste0("gene level data: ", fullData))
message(paste0("tRNA isoform data: ", tRNAData))
message(paste0("normalized data output dir: ", normalizedDir))
message(paste0("PCA data output dir: ", pcaDir))
message("--------------------------------------------------")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN DATA FILES FROM INPUT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Read in netadata
meta <- read.csv(metadata, sep = ",")
rownames(meta) <- meta$Sample_ID

#----- Read in the data files
data <- read.csv(fullData, sep = "\t")
trna <- read.csv(tRNAData, sep = "\t")

#----- Check all samples are shared between data and metadata
check1 <- all(rownames(meta) %in% colnames(data))
check2 <- all(rownames(meta) %in% colnames(trna))
message(paste0("All samples shared between gene-level data and metadata: ", check1))
message(paste0("All samples shared between tRNA-isoform data and metadata: ", check2))

#----- Stop if not
if (check1 != TRUE | check2 != TRUE){
    stop("Metadata samples do not match sample information in data!")
} 

#----- Set colors
myColors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
  "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
  "#9C755F", "#BAB0AC")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# CREATE DESEQ2 OBJECT FOR FULL TRNA + SMRNA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

message("Building DESeq2 dataset for gene-level counts...")
ddsFull <- DESeq2::DESeqDataSetFromMatrix(countData = data,
                                  colData = meta,
                                  design = ~Group)
colData(ddsFull)$Group <- stats::relevel(colData(ddsFull)$Group, ref = refLevel)
levels(colData(ddsFull)$Group)

#----- Run DESeq2 on full genome + tRNA counts
ddsFull <- DESeq2::DESeq(ddsFull)

#----- Extract size factors
full_sizeFactors <- DESeq2::sizeFactors(ddsFull)
write.csv(full_sizeFactors, file = paste0(normalizedDir, "gene_level_counts_size_factors.csv"))
message("Size factors saved...")

#----- Extract normalized counts
normData <- DESeq2::counts(ddsFull, normalized = TRUE) %>%
  as.data.frame()
write.csv(normData, file = paste0(normalizedDir, "normalized_gene_level_counts.csv"))
message("Normalized counts saved...")

#----- Perform rlog transformation
ddsFullRlog <- rlog(ddsFull) 
ddsFullRlogMat <- SummarizedExperiment::assay(ddsFullRlog) %>%
  as.matrix()

#----- Explore variance
varianceFull <- apply(ddsFullRlogMat, 1, var)
fullVars <- sort(varianceFull, decreasing = TRUE)

#----- Convert to dataframe
varianceFull <- data.frame(
  Features = seq_along(fullVars),
  Variance = fullVars
)

# Save the plot as a PNG
fullVar <- ggplot2::ggplot(varianceFull, aes(x = Features, y = Variance)) +
  geom_point() +
  theme_classic() +
  labs(
    title = "tRNA + smRNA Gene Variance Plot",
    x = "Number of Features (tRNAs + smRNAs)",
    y = "Variance") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))
ggplot2::ggsave(fullVar, filename = paste0(pcaDir, "gene_level_variance_plot.png"), width = 6, height = 6)

#----- Generate PCs
generatePCs <- function(MAT, VARS, NFEATURES) {
  #----- Set a variable for the number of genes (features) to be used for PCA and clustering
  var_features_n <- NFEATURES
  #----- Order variance and select the rows (genes) with the most variance
  select <- order(VARS, decreasing = TRUE)[1:var_features_n]
  #----- Subset vsd values for genes by top variance ranks
  vsd_sub <- MAT[select,]
  #----- Transpose the matrix
  vsd_sub <- t(vsd_sub)
  #----- Run principal component analysis
  message(paste0("Running PCA on ", NFEATURES, " most variable features..."))
  pca <- prcomp(vsd_sub)
  #----- extract the variance explained by each PC
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  #----- subset for first 5 elements
  percentVar <- percentVar[1:5]
  #----- add names to the percentVar vector
  message("Percent variations:")
  percentVar <- paste(round(percentVar*100,2), "%", sep = " ")
  names(percentVar) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
  message(percentVar)
  #----- Construct a data frame with PC loadings and add sample labels
  pca_res <- list()
  pca_df <- as.data.frame(pca$x)
  pca_eigenvecs <- as.data.frame(pca$rotation[,1:3])
  pca_res[["Loadings"]] <- pca_df
  pca_res[["Eigenvectors"]] <- pca_eigenvecs
  pca_res[["percent_var"]] <- percentVar
  #----- Return PCA values
  return(pca_res)
}

#----- Calculate PCs and extract loadings
PCsFull <- generatePCs(ddsFullRlogMat, fullVars, 500)
loadingsFull <- PCsFull[[1]]
loadingsFull$Sample <- rownames(loadingsFull)

#----- Ensure data is in same order as metadata
loadingsFull <- loadingsFull[match(rownames(meta), rownames(loadingsFull)),]
loadingsFull$Group <- meta$Group

#----- Extract the rotations
rotationsFull <- PCsFull[[2]]
write.csv(rotationsFull, file = paste0(pcaDir, "gene_level_loadings.csv"))
message("PCA loadings saved...")

#----- Plot PCA
PCAplotFull <- ggplot2::ggplot(loadingsFull, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = myColors) +
  labs(title = "tRNA + smRNA PCA", 
       x = paste0("PC1: ", PCsFull[[3]][1]),
       y = paste0("PC2: ", PCsFull[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        title = element_text(face = "bold"))
ggplot2::ggsave(PCAplotFull, file = paste0(pcaDir, "gene_level_PCA.png"), width = 6, height = 6)
message("Done!")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# TRNA ONLY
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

message("Building DESeq2 dataset for tRNA-isoform counts...")
ddstrna <- DESeq2::DESeqDataSetFromMatrix(countData = trna,
                                          colData = meta,
                                          design = ~Group)
colData(ddstrna)$Group <- stats::relevel(colData(ddstrna)$Group, ref = refLevel)
levels(colData(ddstrna)$Group)

#----- Run DESeq2 on trna genome + tRNA counts
ddstrna <- DESeq2::DESeq(ddstrna)

#----- Extract size factors
trna_sizeFactors <- DESeq2::sizeFactors(ddstrna)
write.csv(trna_sizeFactors, file = paste0(normalizedDir, "tRNA_isotype_counts_size_factors.csv"))
message("Size factors saved...")

#----- Extract normalized counts
normData <- DESeq2::counts(ddstrna, normalized = TRUE) %>%
  as.data.frame()
write.csv(normData, file = paste0(normalizedDir, "normalized_tRNA_isotype_counts.csv"))
message("Normalized counts saved...")

#----- Perform rlog transformation
ddstrnaRlog <- rlog(ddstrna) 
ddstrnaRlogMat <- SummarizedExperiment::assay(ddstrnaRlog) %>%
  as.matrix()

#----- Explore variance
variancetrna <- apply(ddstrnaRlogMat, 1, var)
trnaVars <- sort(variancetrna, decreasing = TRUE)

#----- Convert to dataframe
variancetrna <- data.frame(
  Features = seq_along(trnaVars),
  Variance = trnaVars
)

# Save the plot as a PNG
trnaVar <- ggplot2::ggplot(variancetrna, aes(x = Features, y = Variance)) +
  geom_point() +
theme_classic() +
  labs(
    title = "tRNA Isodecoder Variance Plot",
    x = "Number of Features (tRNAs + smRNAs)",
    y = "Variance") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12))
ggplot2::ggsave(trnaVar, filename = paste0(pcaDir, "tRNA_isotype_variance_plot.png"), width = 6, height = 6)

#----- Calculate PCs and extract loadings
PCstrna <- generatePCs(ddstrnaRlogMat, trnaVars, 500)
loadingstrna <- PCstrna[[1]]
loadingstrna$Sample <- rownames(loadingstrna)

#----- Ensure data is in same order as metadata
loadingstrna <- loadingstrna[match(rownames(meta), rownames(loadingstrna)),]
loadingstrna$Group <- meta$Group

#----- Extract the rotations
rotationstrna <- PCstrna[[2]]
write.csv(rotationstrna, file = paste0(pcaDir, "tRNA_isotype_loadings.csv"))
message("PCA loadings saved...")

#----- Plot PCA
PCAplottrna <- ggplot2::ggplot(loadingstrna, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = myColors) +
  labs(title = "tRNA Isodecoder PCA", 
       x = paste0("PC1: ", PCstrna[[3]][1]),
       y = paste0("PC2: ", PCstrna[[3]][2])) +
  theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none",
        title = element_text(face = "bold"))
ggplot2::ggsave(PCAplottrna, file = paste0(pcaDir, "tRNA_isotype_PCA.png"), width = 6, height = 6)
message("Done!")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# COMBINED PCA PLOT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

combinedPCA <- cowplot::plot_grid(PCAplotFull, PCAplottrna)
ggplot2::ggsave(paste0(pcaDir, "PCA_Analysis_Summary.png"), width = 10, height = 10)


