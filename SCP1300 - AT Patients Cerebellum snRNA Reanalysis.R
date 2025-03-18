# SCP1300 ATM Cerebellum snRNAseq analysis 

# The raw files can be found at this link https://singlecell.broadinstitute.org/single_cell/study/SCP1300/single-cell-atlas-of-the-human-cerebellum-in-ataxia-telangiectasia#study-download

# Install and load necessary packages
install.packages("readr")
library(Seurat)
library(cellranger)
library(DESeq2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gplots)
library(DropletUtils)
library(presto)
library(SingleCellExperiment)
library(readr)
library(shiny)
library(RColorBrewer)
library(cowplot)
library(smplot2)

# SCP1300 Data Loading and Clustering 
setwd("/SCP1300/20241206_SCP1300_Rerun_pct")
counts = readr::read_tsv("20210217_ATCB_expressionmatrix.tsv")
counts = as.data.frame(counts)
rownames(counts) = counts$GENE
counts = counts[,-1]
metadata = read.csv("20210217_ATCB_metadata.csv")
metadata = metadata[-1,]
rownames(metadata) = metadata$NAME

sc_data = CreateSeuratObject(counts = counts, meta.data = metadata, min.features = 200, min.cells = 3)
counts = GetAssayData(object = sc_data, slot = "counts")
metadata = sc_data@meta.data

coord = read.table("20210217_ATCB_clustercoordinates.txt", header = T)
coord = coord[-1,]
rownames(coord) = coord$NAME
coord = coord[,-c(1,4)]
colnames(coord) = c("UMAP_1","UMAP_2")
coord = as.matrix(coord)
class(coord) = "numeric"
sc_data[['umap']] = CreateDimReducObject(embeddings = coord, key = 'UMAP_', assay = 'RNA')

# Figure 4g - UMAP
# UMAP plot
Idents(sc_data) = "celltype"
pdf("UMAP_celltype.pdf")
DimPlot(sc_data, reduction = "umap", raster = FALSE) + theme(aspect.ratio = 1)
dev.off()

# find ATM expression across cell types
DefaultAssay(sc_data) = "RNA"
Idents(sc_data) = "celltype"
Idents(sc_data) = factor(sc_data@active.ident, sort(levels(sc_data@active.ident)))

# Figure 4h - ALDOC expression across all cell types
pdf("ALDOCexpression.pdf")
DotPlot(sc_data, features =c("ALDOC"), scale.by = 'radius',split.by='celltype', dot.scale = 7, cluster.idents = FALSE, col.min = 0, dot.min = 0, cols=c("RdBu")) + theme(axis.text.x=element_text(angle=45,hjust = 0.95, vjust=0.95, size = 7), aspect.ratio = 1/1.5, axis.text.y = element_text(size = 10), legend.title = element_text(size = 10), legend.text = element_text(size = 7))
dev.off()

Purkinje <- subset(sc_data, idents="Purkinje")
Purkinje = NormalizeData(Purkinje)
Purkinje = FindVariableFeatures(Purkinje, selection.method = "vst", nfeatures=2000)
pdf("3_PurkinjeVariableFeaturePlot.pdf", width=8,height=8)
VariableFeaturePlot(Purkinje)
dev.off()

purkinje.genes = rownames(Purkinje)
Purkinje = ScaleData(Purkinje, features=purkinje.genes)
Purkinje = RunPCA(Purkinje, features=(VariableFeatures(object=Purkinje)))
# Selecting the number of PC
pdf("4_PurkinjeElbowPlot.pdf", width=8,height=8)
ElbowPlot(Purkinje)
dev.off()

pct = Purkinje[["pca"]]@stdev / sum(Purkinje[["pca"]]@stdev)*100
cumu = cumsum(pct)
co1 = which(cumu > 90 & pct < 5)[1]
co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)])> 0.1), decreasing = T) [1] +1
min(co1,co2)
# Minimum = 9, so PC9 should be sufficient with the elbow plot too

Purkinje = RunTSNE(Purkinje, reduction = "pca",dims = 1:10,check_duplicates=FALSE)
Purkinje = FindNeighbors(Purkinje, dims=1:10)
Purkinje = FindClusters(Purkinje,resolution = c(0.2,0.3,0.4,0.5,0.6,0.8,1))

# Figure 4i - Cell ratios of ALDOC-positive and ALDOC-negative Purkinje Cells
# Identifying aldoc-positive cells
aldoc_pospc <- WhichCells(Purkinje, expression = ALDOC > 0)
Purkinje$ALDOC <- ifelse(colnames(Purkinje) %in% aldoc_pospc, "ALDOC_pos", "ALDOC_neg")

pdf("ALDOC_expression.pdf",width=4, height=4)
ggplot(Purkinje@meta.data, aes(x=disease__ontology_label,fill=ALDOC)) + geom_bar(position="fill")
dev.off()

# DEG: ALDOC positive PC - Control (NC) vs AT Patients
table(Purkinje@meta.data$ALDOC, Purkinje@meta.data$disease__ontology_label)
#ataxia telangiectasia normal
#ALDOC_neg                   479    367
#ALDOC_pos                   364   1059

Idents(Purkinje) <- "ALDOC"
ALDOC_pos <- subset(Purkinje, idents= "ALDOC_pos")

Idents(ALDOC_pos) = "disease__ontology_label"
aldoc_pos_DEG <- FindMarkers(ALDOC_pos, ident.1="ataxia telangiectasia",ident.2="normal",min.pct=0.25,verbose=TRUE,test.use = "wilcox")
write.table(aldoc_pos_DEG, "aldoc_pos_DEGATvcontrol.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Figure 4j - Volcano Plot of DEGs comparing AT and NC in ALDOC-positive PCs
library(EnhancedVolcano)
pdf(file = "WTvAT_Aldocpos_PC_Volcano_SC.pdf")
EnhancedVolcano(aldoc_pos_DEG,
                lab = NA,
                ylab = bquote(~Log[10]~'P-adj'), 
                pCutoff=0.05,
                FCcutoff = 0.69,
                x = 'avg_log2FC.sc',
                y = 'p_val_adj.sc', 
                xlim = c(-3,3),
                ylim = c(0,20),
                legendLabels = c('Not sig.', 'Log(base2)FC', 'p-adj','p-adj & Log(base2)FC'))
dev.off()

# Supplementary Figure S11d - Figure 4j's counterpart, Volcano Plot of DEGs comparing AT and NC in ALDOC-negative PCs
# DEG: ALDOC negative PC - Control (NC) vs AT Patients
table(Purkinje@meta.data$ALDOC, Purkinje@meta.data$disease__ontology_label)
#ataxia telangiectasia normal
#ALDOC_neg                   479    367
#ALDOC_pos                   364   1059

Idents(Purkinje) <- "ALDOC"
ALDOC_neg <- subset(Purkinje, idents="ALDOC_neg")
Idents(ALDOC_neg) = "disease__ontology_label"
aldoc_neg_DEG <- FindMarkers(ALDOC_neg, ident.1="ataxia telangiectasia",ident.2="normal",min.pct=0.25,verbose=TRUE,test.use = "wilcox")
write.table(aldoc_neg_DEG, "aldoc_neg_DEGATvcontrol.txt", row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

# Figure S11d - Volcano Plot of DEGs comparing AT and NC in ALDOC-negative PCs
library(EnhancedVolcano)
pdf(file = "WTvAT_Aldocneg_PC_Volcano_SC.pdf")
EnhancedVolcano(aldoc_neg_DEG,
                lab = NA,
                ylab = bquote(~Log[10]~'P-adj'), 
                pCutoff=0.05,
                FCcutoff = 0.69,
                x = 'avg_log2FC.sc',
                y = 'p_val_adj.sc', 
                xlim = c(-3,3),
                ylim = c(0,25),
                legendLabels = c('Not sig.', 'Log(base2)FC', 'p-adj','p-adj & Log(base2)FC'))
dev.off()
