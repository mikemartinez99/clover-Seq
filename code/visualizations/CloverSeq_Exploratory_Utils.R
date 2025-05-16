#~~~~~~~~~~~~~~~~~~~~~~~~ README ~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Title: CloverSeq_Exploratory_Utils.R
# Description: Count abundance barplots for tRNAs and PCA
#
# Author: Mike Martinez
# Lab: Genomic Data Science Core
# Project: CloverSeq-Pipeline
# Date created: 5/11/25
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# LOAD LIBRARIES AND SET PATHS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#----- Libraries
library(ggplot2)
library(reshape2)
library(dplyr)
library(magrittr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)
library(cowplot)

#----- Set command line args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: RScript CloverSeq_Exploratory_Utils.R <input counts> <output_dir>")
}

input_tsv <- args[1]
outputDir <- args[2]

if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# READ IN THE TRNA COUNTS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
trna <- read.csv(input_tsv, sep = "\t")


#----- Remove final row of the data
trna <- trna[trna$tRNA_ID != "*",]
trnas <- subset(trna, grepl("^tRNA|^tRX", trna$tRNA_ID))
trnas$Length <- NULL

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GROUPED BY ISOACCEPTOR
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Calculate relative abundance
trnaRelAbund <- trnas %>%
  tidyr::pivot_longer(-c(tRNA_ID, Isoacceptor), names_to = "Sample", values_to = "count") %>%
  dplyr::group_by(Isoacceptor) %>%
  dplyr::mutate(relative_abundance = count / sum(count)) %>%
  dplyr::ungroup()

#----- Calculate absolute abundance
trnaAbAbund <- trnas %>%
  tidyr::pivot_longer(-c(tRNA_ID, Isoacceptor), names_to = "Sample", values_to = "count") %>%
  dplyr::group_by(Isoacceptor) %>%
  dplyr::mutate(absolute_abundance = sum(count)) %>%
  dplyr::ungroup()

#----- Set color palette
color_palette <- c(
  "#E07B39", "#E8C547", "#72A98F", "#2C6E49", "#BC4B51", "#9C0D38", "#543864", "#805E73",
  "#F56476", "#F7B32B", "#376996", "#50808E", "#D74E09", "#A288E3", "#4E8D7C", "#B89D60",
  "#B13E53", "#5C164E", "#D2D4A9", "#F1AB86", "#ACD8AA", "#FF6F59", "#6A0572", "#FFB8D1"
)


#----- Factor the isoacceptors
trnaRelAbund$Isoacceptor <- factor(trnaRelAbund$Isoacceptor, levels = c(
  "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", "Leu",
  "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Sup", "Thr", "Trp", "Tyr", "Val", "Und"
))

trnaAbAbund$Isoacceptor <- factor(trnaAbAbund$Isoacceptor, levels = c(
  "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", "Leu",
  "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Sup", "Thr", "Trp", "Tyr", "Val", "Und"
))

#----- Plot relative abundance barplot
relAbund <- ggplot2::ggplot(trnaRelAbund, aes(x = Isoacceptor, y = relative_abundance*100, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(y = "Relative Abundance",
       x = "",
       fill = "Sample") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 40, hjust = 1))
ggplot2::ggsave(relAbund, file = paste0(outputDir, "tRNA_relAbund_by_sample.png"), width = 10, height = 8)

#----- Plot absolute abundance barplot
absAbund <- ggplot2::ggplot(trnaAbAbund, aes(x = Isoacceptor, y = absolute_abundance, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(y = "Counts",
       x = "",
       fill = "Sample") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 40, hjust = 1))
ggplot2::ggsave(absAbund, file = paste0(outputDir, "tRNA_absAbund_by_sample.png"), width = 10, height = 8)


#----- Combine the plots
bySample <- cowplot::plot_grid(absAbund, relAbund, ncol = 1)
ggplot2::ggsave(bySample, file = paste0(outputDir, "tRNA_abundance_summary_by_sample.png"), width = 10, height = 8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GROUPED BY SAMPLE
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Calculate relative abundance
sampleRelAbund <- trnas %>%
  tidyr::pivot_longer(-c(tRNA_ID, Isoacceptor), names_to = "Sample", values_to = "count") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(relative_abundance = count / sum(count)) %>%
  dplyr::ungroup()

sampleAbAbund <- trnas %>%
  tidyr::pivot_longer(-c(tRNA_ID, Isoacceptor), names_to = "Sample", values_to = "count") %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(absolute_abundance = sum(count)) %>%
  dplyr::ungroup()

#----- Factor the isoacceptors
sampleRelAbund$Isoacceptor <- factor(sampleRelAbund$Isoacceptor, levels = c(
  "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", "Leu",
  "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Sup", "Thr", "Trp", "Tyr", "Val", "Und"
))

sampleAbAbund$Isoacceptor <- factor(sampleAbAbund$Isoacceptor, levels = c(
  "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile", "iMet", "Leu",
  "Lys", "Met", "Phe", "Pro", "SeC", "Ser", "Sup", "Thr", "Trp", "Tyr", "Val", "Und"
))

#----- Plot relative abundance barplot
relAbundISO <- ggplot2::ggplot(sampleRelAbund, aes(x = Sample, y = relative_abundance*100, fill = Isoacceptor)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(y = "Relative Abundance",
       x = "",
       fill = "Isoacceptor") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 40, hjust = 1))
ggplot2::ggsave(relAbundISO, file = paste0(outputDir, "tRNA_relAbund_by_isoacceptor.png"), width = 10, height = 8)


#----- Plot absolute abundance barplot
absAbundISO <- ggplot2::ggplot(sampleAbAbund, aes(x = Sample, y = absolute_abundance, fill = Isoacceptor)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  labs(y = "Counts",
       x = "",
       fill = "Isoacceptor") +
  theme_classic() +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, angle = 40, hjust = 1))
ggplot2::ggsave(absAbundISO, file = paste0(outputDir, "tRNA_absAbund_by_isoacceptor.png"), width = 10, height = 8)


#----- Combine the plots
byIsoAcc <- cowplot::plot_grid(absAbundISO, relAbundISO, ncol = 1)
ggplot2::ggsave(byIsoAcc, file = paste0(outputDir, "tRNA_abundance_summary_by_isoacceptor.png"), width = 10, height = 8)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# EXPLORATORY DATA ANALYSIS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Prepare the data
trna_eda <- trnas
trna_eda$Isoacceptor <- NULL
rownames(trna_eda) <- trna_eda$tRNA_ID
trna_eda$tRNA_ID <- NULL

#----- Step 1: Calculate geometric mean for each gene
geo_mean <- function(x) {
  exp(mean(log(x[x > 0]), na.rm = TRUE))
}
geomeans <- apply(trna_eda, 1, geo_mean)

#----- Calculate size factors
size_factors <- apply(trna_eda, 2, function(col) {
  ratios <- col / geomeans
  median(ratios, na.rm = TRUE)  # Median of ratios
})

#----- Calculate normalized counts
norm_counts <- sweep(trna_eda, 2, size_factors, "/")

#----- Explore variance
variance <- apply(trna_eda, 1, var)
variance <- sort(variance, decreasing = TRUE)
#plot(variance,
 #    las = 1, 
  #   main = "Sample tRNA Gene Expression Variance",
   #  xlab = "nTNRAs",
   #  cex.lab = 1.4,
    # cex.axis = 1.1,
    # font.lab = 2)

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
  print(percentVar)
  #----- Construct a data frame with PC loadings and add sample labels
  pca_res <- list()
  pca_df <- as.data.frame(pca$x)
  pca_eigenvecs <- as.data.frame(pca$rotation[,1:2])
  pca_res[["Loadings"]] <- pca_df
  pca_res[["Eigenvectors"]] <- pca_eigenvecs
  pca_res[["percent_var"]] <- percentVar
  #----- Return PCA values
  return(pca_res)
}

#----- Calculate PCs and extract loadings
PCs <- generatePCs(norm_counts, variance, 100)
loadings <- PCs[[1]]
loadings$Sample <- rownames(loadings)

#----- Plot PCA
PCAplot <- ggplot2::ggplot(loadings, aes(x = PC1, y = PC2, color = Sample, label = Sample)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = color_palette) +
  labs(x = paste0("PC1: ", PCs[[3]][1]),
       y = paste0("PC2: ", PCs[[3]][2])) +
  theme_classic() +
  theme(axis.text = element_text(size = 14, face = "bold"),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none")
ggplot2::ggsave(PCAplot, file = paste0(outputDir, "PCA_Plot.png"), width = 6, height = 6)

### Add some code to get size factors and normalized counts...

