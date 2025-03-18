library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(presto)
library(dplyr)
library(tibble)
library(fgsea)

# Figure 4L to 4O - GSE165371 Data Mining Reanalysis 
# The raw data in the .RDS format can be found at https://singlecell.broadinstitute.org/single_cell/study/SCP795/
setwd("Project/Cerebellum")
CB4AdultMouseOldVersion=readRDS("rawData/cb_annotated_object.RDS") #24409 genes across 611034 samples
CB4AdultMouse <- UpdateSeuratObject(object = CB4AdultMouseOldVersion) #update seurat object from v2 to v3
CB4AdultMouseMyUMAP=RunUMAP(CB4AdultMouse,dims=c(1:50))

colorPalette=colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(unique(CB4AdultMouseMyUMAP$cluster))) 
pdf("CB4AdultMouseCellTypeMyUMAP.pdf",width=8)
DimPlot(CB4AdultMouseMyUMAP,group.by="cluster",reduction="umap")&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

# Extract the Purkinje Cell Cluster
PCs=subset(CB4AdultMouseMyUMAP,cluster %in% "Purkinje") 

# Processing the snRNA data and TSNE plot to show Purkinje cell localisation within Vermis vs Hemisphere
PCs$RegionGroup=ifelse(PCs$region%in% c("I","II","III","CUL","VI","VII","VIII","IX","X"),"Vermis","Hemispere")
PCs<- RunTSNE(PCs, reduction = "pca", dims = 1:50)
pdf("RegionGroup4PCs8TSNE.pdf",width=8)
DimPlot(PCs,group.by="RegionGroup",reduction="tsne",cols=c("LightGrey","Orange"),raster=TRUE,pt.size=1)&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

# TSNE plot to Purkinje cells localisation within the anterior vs posterior vermis
Vermis=subset(PCs,RegionGroup %in% "Vermis")
Vermis$regionGroup=ifelse(Vermis$region %in% c("I","II","III","CUL"),"Anterior","Posterior")
Vermis<- RunTSNE(Vermis, reduction = "pca", dims = 1:50)
pdf("RegionGroup4Vermis8TSNE.pdf",width=8)
DimPlot(Vermis,group.by="regionGroup",reduction="tsne",cols=c("Khaki","MediumTurquoise"),raster=TRUE,pt.size=1)&NoAxes()&theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

# Aldoc expression level in anterior vs posterior vermis
pdf("AldocExprVlnPlot.pdf",width=3,height=3)
VlnPlot(object = Vermis, features = c("Aldoc"),group.by="regionGroup")&theme(axis.title.x = element_blank(),axis.title.y = element_blank())
dev.off()

# Aldoc+ve vs Aldoc-ve Purkinje cell ratios within the anterior and posterior vermis
Vermis$AdlocGroup=ifelse(Vermis$subcluster %in% c("Purkinje_Anti_Aldoc_1","Purkinje_Anti_Aldoc_2"),"AldocNeg","AldocPos")
tmp=paste0(Vermis$regionGroup,"_",Vermis$AdlocGroup,sep="")
tmp=data.frame(table(tmp))
colnames(tmp)=c("Group","GroupCount")
TmpInfo=as.matrix(do.call(rbind, strsplit(as.character(tmp$Group),'_')))
resultTmp=data.frame(Region=TmpInfo[,1],Aldoc=TmpInfo[,2],"Number"=tmp$GroupCount)
RegionCount=data.frame(table(Vermis$regionGroup))
colnames(RegionCount)=c("Region","RegionCount")
Result=merge(resultTmp,RegionCount,by="Region")
g=ggplot(Result, aes(Region, Number, fill=Aldoc)) +
  geom_bar(stat="identity",position="fill") +
  scale_y_continuous(expand=c(0,0))+theme_bw()+
  scale_fill_manual(values=c("SlateBlue","Salmon"))+
  theme(axis.text.x = element_text(angle = 0))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())
pdf("AldocDisInRegion.pdf",width=6)
print(g)
dev.off()

# Metabolic pathway analysis of Purkinje cells found in anterior lobe vs posterior lobe using KEGG pathway enrichment analysis
Metabolism=read.table("MetabolismPathway4MouseGene.txt",header=T,sep="\t")
fgsea_sets<- split(Metabolism$Symbol,Metabolism$PathwayName)
VermisRegionDEG <- wilcoxauc(Vermis, 'regionGroup')
table(VermisRegionDEG$group)
regionGroup="Posterior"
clusterCell<- VermisRegionDEG %>% dplyr::filter(group == regionGroup) %>% arrange(desc(logFC)) %>% dplyr::select(feature, logFC)
ranks<- deframe(clusterCell)
fgseaRes<- fgseaMultilevel(fgsea_sets, stats = ranks,eps=0,minSize = 10)
fwrite(fgseaRes, file="VermisRegionKEGG4Graph8GSEA.txt", sep="\t", sep2=c("", " ", ""))

data=read.table("Metabolism_PostvsAnter_PCsVermis.txt",header=T,sep="\t") 
data=data[order(data$pval),]
data=data[1:10,]
PathwayList=factor(data$pathway,levels=rev(unique(data$pathway)))
pdf("VermisRegionKEGG4Graph8GSEATop10Terms.pdf",height=3,width=5)
ggplot(data, aes(x=-log10(pval), y=PathwayList,color=NES,size=size)) +
  geom_point() + scale_color_gradient2(low="navy",mid="white",high = "red")+
  theme_bw()
dev.off()
