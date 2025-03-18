# Analysis of Bulk RNA extracted from mouse cerebellum

# Install and load packages for analysis
library(BiocManager)

BiocManager::install("sva")
library(sva)
library(clusterProfiler)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(org.Mm.eg.db)
library(biomaRt)
library(devtools)
library(GeneStructureTools)
library(EnhancedVolcano)

# Differential Expression Analysis using DESeq2
data <- read.table("Raw.txt",sep="\t",header=T,row.names=1)
meta <- read.table("Meta.txt",sep="\t",header=T,row.names=1)

dds <- DESeqDataSetFromMatrix(countData=round(data), 
                                   colData=meta, 
                                   design=~factor(Genotype))
dds
dds <- DESeq(dds)
res <- results(dds)
head(results(dds,tidy=TRUE))
summary <- summary(res)
write.table(res, "DESeq_results.txt", quote = FALSE, sep="\t")
normalised_counts <- counts(dds,normalized = TRUE)
write.table(normalised_counts, file="NormCounts.txt", sep="\t", quote=F, col.names=NA)
sum(res$padj<0.05, na.rm=T) #1552

# Figure S5d - Volcano Plot
pdf(file = "Volcano_Plot.pdf")
EnhancedVolcano(res,
                lab = NA,
                ylab = bquote(~Log[10]~'P-adj'), 
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Control vs ATM-KO; Novogene',
                pCutoff = 0.0276831034463070,
                xlim = c(-4,4),
                ylim = c(0,8),
                legendLabels = c('Not sig.', 'Log(base2)FC', 'p-adj','p-adj & Log(base2)FC'))
dev.off()

# Figure S5d - PCA Plot
vsdata <- vst(Novo_dds, blind=FALSE)
pdf(file="PCA Plot.pdf")
plotPCA(vsdata,intgroup="Genotype")
dev.off()

# Labels on PCA Plot
pdf(file = "2b_PCA Plot.pdf")
z <- plotPCA(vsdata,intgroup="Genotype")
nudge <- position_nudge(y = 1)
z + geom_text(aes(label = name), position = nudge)
dev.off()

# Figure S5e - MA Plot
pdf("MA Plot.pdf")
plotMA(
  Novo_res,
  alpha = 0.05,
  main = "MA Plot",
  xlab = "mean of normalized counts",
  ylim = c(-6,6),
  colNonSig = "gray60",
  colSig = "blue",
  colLine = "red",
  returnData = FALSE,
  MLE = FALSE,
)
dev.off()

# Figure S5f - Histogram of p-values
pdf(file="Histogram of p-values.pdf")
hist(res$pvalue, main = "Histogram of p-values", breaks=20,col = "light grey")
dev.off()
