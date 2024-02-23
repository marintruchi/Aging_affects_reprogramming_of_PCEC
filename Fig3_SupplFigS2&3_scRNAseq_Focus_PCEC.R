# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to explore the subset of PCEC and produce the plots of figure 3 and supplementary figure S2 & S3
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")

#----------------------------------------------------------------------------------------------------------#
#-------------- Identify PCEC subpopulations signature and compare them with public datasets --------------#
#----------------------------------------------------------------------------------------------------------#

#Loading dataset
#For convenience, load directly the seurat object obtained after the processing and integration steps 
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(aggr) <- "celltype"
group_order <- c("PBS","D14","D28","D60")
celltype_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                    "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                    "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                    "Prolif. EC","gCap","aCap","Venous EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                    "Mesothelial cells","AT1","AT2","Multiciliated cells")

#Colors for plots
colors.group <- c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c")
colors.Bleo<- c("Prolif. AM"="#a6d8f5","AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0","IM"="#225EA8","LCM"="#0e3378",
                "Prolif. DC"="#ffd1ec","cDC1"="#e670b4","cDC2"="#9c119c","Mature DC"="#960c5a","pDC"="#de439a","cMonocytes"="#66119c","ncMonocytes"="#282b5c","Neutrophils"="#4300de","Neutrophil-like Monocytes"="#38015c","Mast cells"="#9315ed","Platelets"="#de1212",
                "Prolif. T cells"="#F7FCB9","CD4 T cells"="#D9F0A3","Reg T cells"="#ADDD8E","T helper cells"="#78C679","Gamma delta T cells"="#41AB5D","ILC2"="#739154","CD8 T cells"="#006837" ,"NKT1"="#004529","NKT2"="#01736b","NK cells"="#46c798","B cells"="#ebd513","Plasmocytes"="#ccbf47",
                "Prolif. EC"="#fee0d2","gCap"="#911414","aCap"="#5e1e33","Venous EC"="#f2ae5a","Arterial EC"="#e36805","Lymphatic EC"="#7d3f15","Fibroblasts"="#d3d5db","Pericytes"="#806c6c",
                "Mesothelial cells"="#120c04","AT1"="#8184b8","AT2"="#3e447a","Multiciliated cells"="#1486f7")



#Subsetting within alveolar capillary cells 
DefaultAssay(aggr) <- "integrated"
Idents(aggr) <- "celltype"
SUB <-subset(x = aggr, idents = c("aCap","gCap"))
DimPlot(object = SUB, reduction = 'umap',label = T)
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.6)
DimPlot(object = SUB, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = SUB, reduction = 'umap',group.by = "group",label=F,cols = colors.group)+NoLegend()

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(4))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Lrg1+ gCap')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(8))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'SV EC')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(10))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Lrg1+ aCap')
DimPlot(object = aggr, reduction = 'umap',label=T)+NoLegend()

aggr$cellstate <- Idents(aggr)
aggr$cellstate <- ifelse(aggr$cellstate="Venous EC", "PV EC",aggr$cellstate)
cellstate_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                     "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                     "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                     "Prolif. EC","SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap" ,"PV EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                     "Mesothelial cells","AT1","AT2","Multiciliated cells")
aggr$cellstate <- factor(aggr$cellstate,levels = cellstate_order)
Idents(aggr) <- "cellstate"




#Capillary endothelial cellstates markers
Idents(aggr) <- "cellstate"
SUB <-subset(x = aggr, idents = c("aCap","gCap","Lrg1+ aCap","Lrg1+ gCap","SV EC","Prolif. EC","PV EC","Arterial EC"))
SUB$cellstate <- factor(SUB$cellstate,levels = c("SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Prolif. EC","PV EC","Arterial EC"))
subcolors <- c("SV EC"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414","Lrg1+ aCap"="#FFD92F","aCap"="#5e1e33","Prolif. EC"="#fee0d2","PV EC"="#f2ae5a","Arterial EC"="#e36805")
SUB <- RunUMAP(SUB,dims=1:80)
DimPlot(object = SUB,group.by = "cellstate", reduction = 'umap',label = F,cols = subcolors,pt.size = 0.2)+NoLegend()+NoAxes()+ggtitle("")
ggsave("Fig3_A.pdf", width = 8, height = 8, units = c("cm"), dpi = 200)


#------------------------------------------------------------------------------------------
#PCEC markers
DefaultAssay(SUB) <- "RNA"
SUB <- SUB %>% NormalizeData()
markers <- FindAllMarkers(SUB,logfc.threshold = 0.5,min.pct = 0.3,only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2*markers$avg_log2FC

#Heatmap of marker genes
SUB$group.cs <- paste(SUB$group,SUB$cellstate,sep = "-")
Idents(SUB) <- "group.cs"
average <- AverageExpression(SUB, return.seurat = T)

data <- GetAssayData(average,slot = "data")
ngenes <- unlist(findnoisygenes.mm(SUB))
markers$cluster <- factor(markers$cluster,levels = c("SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Prolif. EC","PV EC","Arterial EC"))
genes.to.plot <- markers %>% filter(gene%ni%ngenes) %>% group_by(cluster) %>% top_n(8,avg_log2FC) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)

sample.order <- levels(aggr$group)
cellstate_order<- levels(markers$cluster)
data <- data[match(unique(genes.to.plot$gene),rownames(data)),]
annotation <- data.frame(group=colnames(data),sample=gsub("-.*","",colnames(data)),cellstate=gsub(".*-","",colnames(data)))
annotation$sample <- factor(annotation$sample,levels = sample.order)
annotation$cellstate <- factor(annotation$cellstate,levels = cellstate_order)
annotation <-arrange(annotation,cellstate,sample)
data <- data[,match(annotation$group,colnames(data))]
annotation <- annotation %>% column_to_rownames(.,"group")
data <- flexible_normalization(data)
data[data < -2] <- -2
data[data > 2] <- 2

pdf("Fig3_B.pdf", width = 10, height = 3.5)
pheatmap(t(data),annotation_row  = annotation,show_rownames = F,cluster_cols = F,cluster_rows = F,angle_col =90,fontsize_col = 10,
         annotation_colors = list(cellstate= subcolors,sample=colors.group))
dev.off()

genes.to.plot <- markers %>% filter(gene%ni%ngenes) %>% filter(cluster%in%c("SV EC","Lrg1+ gCap","Lrg1+ aCap")) %>% group_by(cluster) %>% top_n(15,ratio) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)
library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- genes.to.plot$gene
df <- listAttributes(mart)
G_list <- getBM(filters= "mgi_symbol", attributes= c("description","external_gene_name"),values=genes,mart= mart) %>% mutate(gene=external_gene_name)

genes.to.plot <- merge(genes.to.plot,G_list,by.x="gene")
genes.to.plot <- genes.to.plot  %>% dplyr::select(cluster,gene,description) %>% arrange(cluster,gene)
genes.to.plot$description <- gsub("\\[Source.*","",genes.to.plot$description)


#------------------------------------------------------------------------------------------
#Vlnplot CEC markers
VlnPlot(SUB,features=c("Lrg1"),group.by = "cellstate",pt.size = 0,cols =  subcolors,sort = F,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("Fig3_C_Lrg1.pdf", width = 7, height = 4, units = c("cm"), dpi = 100)
VlnPlot(SUB,features=c("Col15a1"),group.by = "cellstate",pt.size = 0,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("Fig3_C_Col15a1.pdf", width = 7, height = 4, units = c("cm"), dpi = 100)
VlnPlot(SUB,features=c("Serpine1"),group.by = "cellstate",pt.size = 0,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("Fig3_C_Serpine1.pdf", width = 7, height = 4, units = c("cm"), dpi = 100)
VlnPlot(SUB,features=c("Ednrb"),group.by = "cellstate",pt.size = 0,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("Fig3_C_Ednrb.pdf", width = 7, height = 4, units = c("cm"), dpi = 100)
VlnPlot(SUB,features=c("Aplnr"),group.by = "cellstate",pt.size = 0,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
ggsave("Fig3_C_Aplnr.pdf", width = 7, height = 4, units = c("cm"), dpi = 100)

#------------------------------------------------------------------------------------------
#CEC relative proportions
celltype_order <- c("SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Prolif. EC","PV EC","Arterial EC")
metadata <- SUB@meta.data %>% filter(cellstate %in% celltype_order)
metadata$cellstate <- factor(metadata$cellstate, levels = celltype_order)
metadata = as.data.frame.matrix(table(metadata$group, metadata$cellstate))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$group = rownames(metadata)
metadata$group =factor(metadata$group, levels = group_order)
metadata = gather(metadata, cellstate, percentage, 'SV EC':'Arterial EC')
metadata$cellstate <- factor(metadata$cellstate, levels =celltype_order)

ggplot(metadata, aes(x = group, y = percentage, fill = cellstate)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = subcolors) +
  ggtitle("CEC relative proportions")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+ggtitle("")+NoLegend()
ggsave("Fig3_D.pdf", width = 5, height = 8, units = c("cm"), dpi = 200)


#------------------------------------------------------------------------------------------
#Comparison between SVEC markers in human and mouse fibrosis

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SVEC.list <- c("Col15a1","Hspg2","Nrp2","Gja1","Vwa1","Lamb1","Spry1","Filip1","Ndrg1","Nr5a2","Galnt15","Meox1")

Idents(aggr) <- "cellstate"
aggr <- subset(aggr,idents=c('SV EC',"PV EC","Arterial EC","gCap","aCap","Lrg1+ gCap","Lrg1+ aCap"))
DefaultAssay(aggr) <- "RNA"
aggr <- aggr %>% NormalizeData() %>% ScaleData()
aggr$cellstate <- factor(aggr$cellstate,levels = rev(c('SV EC',"PV EC","Arterial EC","gCap","aCap","Lrg1+ gCap","Lrg1+ aCap")))
DotPlot(aggr, features = SVEC.list,group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)+ylab("")+xlab("")+NoLegend()
ggsave("SupplFig2_A_mice.pdf", width = 9, height = 5, units = c("cm"), dpi = 200)

#Seurat object of the Habermann et al. dataset reintegrated and reannoted
Habermann <- readRDS("~/Habermann_integrated.rds")
Habermann <- subset(Habermann,idents=c('SV EC',"PV EC","Arterial EC","aCap","gCap"))
DefaultAssay(Habermann) <- "RNA" 
Habermann <- Habermann %>% NormalizeData() %>% ScaleData()
Habermann$celltype_MT <- factor(Habermann$celltype_MT,levels = rev(c('SV EC',"PV EC","Arterial EC","gCap","aCap")))
DotPlot(Habermann, features = str_to_upper(SVEC.list),group.by ="celltype_MT",dot.scale = 5,cols = "RdBu",assay = "RNA") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)+ylab("")+xlab("")+NoLegend()
ggsave("SupplFig2_A_human.pdf", width = 9, height = 5, units = c("cm"), dpi = 200)


#------------------------------------------------------------------------------------------
#Lrg1 expression in all celltypes at all time points
colors.Bleo<- c("Prolif. AM"="#a6d8f5","AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0","IM"="#225EA8","LCM"="#0e3378",
                "Prolif. DC"="#ffd1ec","cDC1"="#e670b4","cDC2"="#9c119c","Mature DC"="#960c5a","pDC"="#de439a","cMonocytes"="#66119c","ncMonocytes"="#282b5c","Neutrophils"="#4300de","Neutrophil-like Monocytes"="#38015c","Mast cells"="#9315ed","Platelets"="#de1212",
                "Prolif. T cells"="#F7FCB9","CD4 T cells"="#D9F0A3","Reg T cells"="#ADDD8E","T helper cells"="#78C679","Gamma delta T cells"="#41AB5D","ILC2"="#739154","CD8 T cells"="#006837" ,"NKT1"="#004529","NKT2"="#01736b","NK cells"="#46c798","B cells"="#ebd513","Plasmocytes"="#ccbf47",
                "SV EC"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414","Lrg1+ aCap"="#FFD92F","aCap"="#5e1e33","Prolif. EC"="#fee0d2","PV EC"="#f2ae5a","Arterial EC"="#e36805","Lymphatic EC"="#7d3f15","Fibroblasts"="#d3d5db","Pericytes"="#806c6c",
                "Mesothelial cells"="#120c04","AT1"="#8184b8","AT2"="#3e447a","Multiciliated cells"="#1486f7")
cellstate_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                     "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                     "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                     "Prolif. EC","SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap" ,"PV EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                     "Mesothelial cells","AT1","AT2","Multiciliated cells")


aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")

aggr$celltype_global <- ifelse(aggr$cellstate%in%c("Prolif. DC","cDC1","cDC2","Mature DC","pDC"),"DCs",
                               ifelse(aggr$cellstate%in%c("cMonocytes","ncMonocytes","Neutrophil-like Monocytes"),"Monocytes",
                                      ifelse(aggr$cellstate%in%c("Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes"),"Lymphocytes",as.character(aggr$cellstate))))   
aggr$celltype_global <-factor(aggr$celltype_global,levels = c("aCap","gCap","Prolif. EC","SV EC","Lrg1+ gCap","Lrg1+ aCap","Arterial EC","PV EC","Lymphatic EC","Fibroblasts","Pericytes","Mesothelial cells","AT1","AT2","Multiciliated cells",
                                                              "Prolif. AM","AM1","AM2","AM3","IM","LCM","DCs","Monocytes","Mast cells","Neutrophils","Platelets","Lymphocytes"))

Idents(aggr) <- "group"
SUB <-subset(x = aggr, subset = group=="PBS")
SUB <- SUB %>% NormalizeData()
SUB1 <-subset(x = aggr, idents = "D14")
SUB1 <- SUB1 %>% NormalizeData()
SUB2 <-subset(x = aggr, idents = "D28")
SUB2 <- SUB2 %>% NormalizeData()
SUB3 <-subset(x = aggr, idents = "D60")
SUB3 <- SUB3 %>% NormalizeData()
VlnPlot(SUB,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)
ggsave("SupplFig2_B_PBS.pdf", width = 15, height = 9, units = c("cm"), dpi = 200)


VlnPlot(SUB1,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)
ggsave("SupplFig2_C_D14.pdf", width = 15, height = 9, units = c("cm"), dpi = 200)


VlnPlot(SUB2,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)
ggsave("SupplFig2_C_D28.pdf", width = 15, height = 9, units = c("cm"), dpi = 200)


VlnPlot(SUB3,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)
ggsave("SupplFig2_C_D60.pdf", width = 15, height = 9, units = c("cm"), dpi = 200)


#------------------------------------------------------------------------------------------
#Lrg1 Expression in CEC from Strunz et al.

#Seurat object of the Strunz et al. dataset reintegrated and reannoted
Strunz <- readRDS(file = "~/Strunz_RPCA.rds")
SUB <-subset(x = Strunz, idents = c("9","18"))
DefaultAssay(SUB) <- "integrated" 
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.5)
SUB <-subset(x = SUB, idents = c("3","0","2","5"))
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.5)

#Project annotated dataset on this one
ref <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(ref) <- "cellstate"
celltype_order <- c("SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","PV EC", "Arterial EC", "Lymphatic EC")
ref <- ref %>% subset(idents=celltype_order)
DefaultAssay(ref) <- "integrated"
anchors <- FindTransferAnchors(reference = ref, query = SUB,dims = 1:80, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = ref$cellstate,dims = 1:80)
SUB <- AddMetaData(SUB, metadata = predictions)

p0 <- DimPlot(object = SUB,group.by = "seurat_clusters", label=TRUE)+ ggtitle("Clusters")+NoAxes()+NoLegend()
p1 <- DimPlot(object = SUB, group.by = "predicted.id",cols = colors.Bleo)+ ggtitle("Celltype")+NoAxes()+NoLegend()
plot_grid(p0,p1,ncol = 2)
colors.sample <- c("PBS"="#8dbee0","d3"="#e08d8d","d7"="#e05c5c","d10"="#e62727","d14"="#946666","d21"="#913c3c","d28"="#8f1414")
DimPlot(object = SUB, group.by = "grouping",cols = colors.sample)+ ggtitle("Celltype")+NoAxes()+NoLegend()
DefaultAssay(SUB) <- "RNA"
SUB <- SUB %>% NormalizeData()
FeaturePlot(SUB,features = "Col4a1")
ngenes <- rownames(SUB)[which(rownames(SUB)%ni%unlist(findnoisygenes.mm(SUB)))]
markers <- FindAllMarkers(SUB,features = ngenes,logfc.threshold = 0.4,min.pct = 0.4,only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2*markers$avg_log2FC
genes.to.plot <- markers %>% group_by(cluster) %>% top_n(50,ratio) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)
VlnPlot(SUB,features="nCount_RNA",cols =  colors.Bleo,pt.size = 0)

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(0))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'gCap')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(1))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'PV EC')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(2))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'aCap')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(3))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'Lymphatic EC')
DimPlot(object = SUB, reduction = 'umap',label = T)

SUB$grouping <- factor(SUB$grouping, levels = c("PBS","d3","d7","d10","d14","d21","d28"))
aCap <- subset(SUB,ident=c('aCap'))
gCap <- subset(SUB,ident=c('gCap'))

VlnPlot(aCap,features=c("Lrg1"),group.by = "grouping",pt.size = 0,cols =  colors.sample,sort = F)+FontSize(x.text = 10,y.text = 8)+geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=1.5))+FontSize(x.text = 8,y.text = 8,y.title = 0,x.title=0,main = 0)+NoLegend()
ggsave("SupplFig2_F_aCap.pdf", width = 7, height = 4.5, units = c("cm"), dpi = 200)

VlnPlot(gCap,features=c("Lrg1"),group.by = "grouping",pt.size = 0,cols =  colors.sample,sort = F)+geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=1.5))+FontSize(x.text = 8,y.text = 8,y.title = 0,x.title=0,main = 0)+NoLegend()
ggsave("SupplFig2_F_gCap.pdf", width = 7, height = 4.5, units = c("cm"), dpi = 200)



#------------------------------------------------------------------------------------------
# LRG1 Expression
EC_order <- c("SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","PV EC", "Arterial EC")
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
sub <- subset(aggr,subset = cellstate %in% EC_order)
DefaultAssay(sub) <-"RNA"

sub$Condition <- as.character(sub$condition)
sub$celltype <- gsub(" ","",sub$cellstate)
sub$id <- paste0(sub$celltype,"_",sub$age,".",sub$group)
select.cells <- colnames(subset(sub,subset = Lrg1 >= 1))
sub1 <- subset(sub, cells = select.cells) 

Lrg1 <- data.frame(table(sub1$id))
table(sub1$id)

Lrg1 <-full_join(data.frame(table(sub$id)),data.frame(table(sub1$id)),by="Var1") %>% mutate(percentage=Freq.y/Freq.x*100)
sub1 <- sub %>% NormalizeData() %>% ScaleData()
Idents(sub1) <- "id"
sub1 <- AverageExpression(sub1, return.seurat = T,assays = "RNA")
sub1 <- data.frame(Var1=names(sub1@assays[["RNA"]]@scale.data["Lrg1",]),Lrg1_expr=sub1@assays[["RNA"]]@scale.data["Lrg1",]) %>% inner_join(Lrg1) %>% mutate(id=Var1) %>% dplyr::select(id,Lrg1_expr,percentage)
sub1[is.na(sub1)] <- 0

#Split ID in groupe and Condition
sub1$group <- gsub("\\..*","",sub1$id)
sub1$Condition <- factor(gsub(".*\\.","",sub1$id),levels = rev(c("PBS","D14","D28","D60")))
groups = c("SVEC_Young","SVEC_Old","Lrg1+gCap_Young","Lrg1+gCap_Old","gCap_Young","gCap_Old","Lrg1+aCap_Young","Lrg1+aCap_Old","aCap_Young","aCap_Old","PVEC_Young","PVEC_Old","ArterialEC_Young","ArterialEC_Old")

sub1$group <- factor(sub1$group,levels = groups)

ggplot(sub1,aes(x=group,y=Condition,size=percentage)) + 
  geom_point(aes(col=Lrg1_expr),alpha=1) +
  scale_colour_gradientn(colours = c("#67001F", "#B2182B" ,"#D6604D" ,"#F4A582", "#FDDBC7", "#F7F7F7" ,"#D1E5F0", "#92C5DE" ,"#4393C3"),values = scales::rescale(c(max(sub1$Lrg1_expr),2,1.5,1,0.5, 0, -0.5, -1,min(sub1$Lrg1_expr))))+
  labs(y="time point")+theme_classic()+rotate_x_text(45)+ylab("")+xlab("")
ggsave("SupplFig2_E.pdf", width = 13, height = 5.5, units = c("cm"), dpi = 100)




#------------------------------------------------------------------------------------------
#Differential Expression Analyses using DESeq2 on pseudo-bulks for Lrg1+ PCEC subpopulations

#For Lrg1+ gCap and Lrg1+ aCap
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap") & group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

young <- subset(x = aggr, subset = age == c("Young"))
old <- subset(x = aggr, subset = age == c("Old"))

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

cts <- cts.split$`gCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")

write.xlsx(DEG.psdblk.gCap.mice,file = "~/DEG.pseudobulks/Lrg1pos_gCap_pseudobulk.xlsx")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% filter(padj<0.05) 

cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.aCap.mice,file = "~/DEG.pseudobulks/Lrg1pos_aCap_pseudobulk.xlsx")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% filter(padj<0.05) 


#For SVEC (merge SVEC and PVEC from PBS to have enough cells for pseudobulks)
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SVEC <-subset(x = aggr, subset = cellstate %in% c("SV EC")& group %ni% c("PBS"))
PVEC <-subset(x = aggr, subset = cellstate %in% c("PV EC") & group %in% c("PBS"))
aggr <- merge(SVEC,PVEC)
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.SVEC.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.SVEC.mice,file = "~/DEG.pseudobulks/SVEC_pseudobulk.xlsx")
DEG.psdblk.SVEC.mice <- DEG.psdblk.SVEC.mice %>% filter(padj<0.05) 


#For PV EC
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("PV EC")& group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.PVEC.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.PVEC.mice,file = "~/DEG.pseudobulks/PVEC_pseudobulk.xlsx")
DEG.psdblk.PVEC.mice <- DEG.psdblk.PVEC.mice %>% filter(padj<0.05)


#--------------------------------------------------------------------------------------------------------------------#
#--------------- Perform Functional enrichment and NicheNet analysis on PCEC subpopulations signatures---------------#
#--------------------------------------------------------------------------------------------------------------------#

# From the differential analyses tables generated above, Functional enrichment analyses were performed using IPA

# Enrichment analysis outputs from IPA
CP.order <- c("Pulmonary Fibrosis Idiopathic Signaling Pathway","Wound Healing Signaling Pathway","ILK Signaling","HIF1Î± Signaling","Actin Cytoskeleton Signaling","EIF2 Signaling","Pulmonary Healing Signaling Pathway","ID1 Signaling Pathway","mTOR Signaling","VEGF Signaling","Apelin Endothelial Signaling Pathway","FAK Signaling","Glycolysis I","TGF-Î² Signaling","NF-ÎºB Signaling","IL-8 Signaling","WNT/Î²-catenin Signaling","Hypoxia Signaling in the Cardiovascular System")
aCap <- read.xlsx("~/aCap_IPA.xlsx",sheet="CP") %>% filter(Ingenuity.Canonical.Pathways %in% c(CP.order)) %>% dplyr::select(Ingenuity.Canonical.Pathways,"-log(p-value)","z-score") %>% mutate(subpop="aCap")
gCap <- read.xlsx("~/gCap_IPA.xlsx",sheet="CP") %>% filter(Ingenuity.Canonical.Pathways %in% c(CP.order)) %>% dplyr::select(Ingenuity.Canonical.Pathways,"-log(p-value)","z-score") %>% mutate(subpop="gCap")
SVEC <- read.xlsx("~/SVEC_IPA.xlsx",sheet="CP") %>% filter(Ingenuity.Canonical.Pathways %in% c(CP.order)) %>% dplyr::select(Ingenuity.Canonical.Pathways,"-log(p-value)","z-score") %>% mutate(subpop="SV EC")
PVEC <- read.xlsx("~/PVEC_IPA.xlsx",sheet="CP") %>% filter(Ingenuity.Canonical.Pathways %in% c(CP.order)) %>% dplyr::select(Ingenuity.Canonical.Pathways,"-log(p-value)","z-score") %>% mutate(subpop="PV EC")

CP <- Reduce(rbind,list(aCap,gCap,SVEC,PVEC))
colnames(CP) <- c("CP","pval","z-score","subpop")
CP$pval[CP$pval>6] <- 6
CP$CP <- factor(CP$CP,levels = rev(CP.order))

ggplot(CP, aes(pval, CP, color = subpop)) +
  geom_point() +ggtitle("") + ylab("") + scale_color_manual(values = c("SV EC"="#A6D854","gCap"="#fa3737","aCap"="#FFD92F","PV EC"="#f2ae5a"))+
  theme(axis.text.y = element_text(size = rel(.75)),axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())+theme_bw()+NoLegend()+ geom_vline(xintercept = -log10(0.05), linetype="dotted", color = "black", size=0.5)
ggsave("Fig3_E.pdf", width = 12, height = 8, units = c("cm"), dpi = 150)


DF.order <- c("Invasion of cells","Vasculogenesis","Cell movement of endothelial cells","Leukocyte migration","Cell-cell contact", "Morphology of vessel",
              "Sprouting","Angiogenesis", "Cell survival","Metabolism of protein","Formation of cytoskeleton",
              "Synthesis of protein","Inflammatory response","Chemotaxis",
              "Quantity of connective tissue cells","Quantity of metal ion", "Growth of connective tissue")

aCap <- read.xlsx("~/aCap_IPA.xlsx",sheet="DF") %>% filter(Diseases.or.Functions.Annotation %in% c(DF.order)) %>% dplyr::select(Diseases.or.Functions.Annotation,"p-value") %>% mutate(subpop="aCap")
gCap <- read.xlsx("~/gCap_IPA.xlsx",sheet="DF") %>% filter(Diseases.or.Functions.Annotation %in% c(DF.order)) %>% dplyr::select(Diseases.or.Functions.Annotation,"p-value") %>% mutate(subpop="gCap")
SVEC <- read.xlsx("~/SVEC_IPA.xlsx",sheet="DF") %>% filter(Diseases.or.Functions.Annotation %in% c(DF.order)) %>% dplyr::select(Diseases.or.Functions.Annotation,"p-value") %>% mutate(subpop="SV EC")
PVEC <- read.xlsx("~/PVEC_IPA.xlsx",sheet="DF") %>% filter(Diseases.or.Functions.Annotation %in% c(DF.order)) %>% dplyr::select(Diseases.or.Functions.Annotation,"p-value") %>% mutate(subpop="PV EC")

DF <- Reduce(rbind,list(aCap,gCap,SVEC,PVEC))
colnames(DF) <- c("DF","pval","subpop")
DF$pval <- -log10(DF$pval)
DF$pval[DF$pval>10] <- 10
DF$DF <- factor(DF$DF,levels = rev(DF.order))

ggplot(DF, aes(pval, DF, color = subpop)) +
  geom_point() +ggtitle("") + ylab("") + scale_color_manual(values = c("SV EC"="#A6D854","gCap"="#fa3737","aCap"="#FFD92F","PV EC"="#f2ae5a"))+
  theme(axis.text.y = element_text(size = rel(.75)),axis.ticks.y = element_blank(),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())+theme_bw()+NoLegend()+ geom_vline(xintercept = -log10(0.05), linetype="dotted", color = "black", size=0.5)
ggsave("Fig3_F.pdf", width = 14, height = 8, units = c("cm"), dpi = 150)


#------------------------------------------------------------------------------------------
#Expression of Hypoxia markers
#Loading dataset
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
markers_hypox <- c("Aldoa","Ankrd37","Apln","App","Bhlhe40","Bnip3","Dnah11","Ero1l","Gadd45b","Gapdh","Gpi1","Hif1a","Idh2","Ldha","Pdgfb","Pkm","Plod1","Serpine1")
SUB <-subset(x = aggr, subset = cellstate %in% c("aCap","gCap","Lrg1+ aCap","Lrg1+ gCap","SV EC"))
EC_order <- c("Lrg1+ aCap","aCap","Lrg1+ gCap","gCap","SV EC")
SUB$cellstate <- factor(SUB$cellstate,levels = rev(EC_order))
SUB <- SUB %>% NormalizeData() 
DotPlot(SUB, features = markers_hypox,group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+NoLegend()
ggsave("Fig3_G_hypoxic_signature.pdf", width = 10, height = 4, units = c("cm"), dpi = 150)

#Expression of Pro-angiogenic markers
angiogenic.features <- c("Dll4","Esm1","Hspg2","Plxnd1","Rgcc","Aplnr","Smad1","Sparc","Col4a1","Col4a2","Flt1","Xbp1","Btg1","Rhoj","Cxcl12","Sox17","Adamts1","Adam15","Vegfa","Ackr3","Lrg1","Nrp1","Kdr","Hspb1","Emp2","Cd34","Serpine1","Unc5b","Hif1a","Tgfbr2")
plot <- DotPlot(SUB, features = rev(angiogenic.features),group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+NoLegend()
plot
ggsave("Fig3_G_angiogenic_markers.pdf", width = 15, height = 4, units = c("cm"), dpi = 150)

gCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ gCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
aCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ aCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
SVEC.angiogenic.features <- plot[["data"]] %>% filter(id=="SV EC"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()


#------------------------------------------------------------------------------------------
#NicheNet Analysis
library(nichenetr)
library(circlize)

#DATA PREPARATION
#-------------------------------------------------------------------------------

#Load Seurat Object
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
select.cells <- names(Idents(aggr)[which(Idents(aggr)%in%c("Prolif. DC","cDC1","cDC2","pDC","Mature DC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Platelets"))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Myeloid cells')
select.cells <- names(Idents(aggr)[which(Idents(aggr)%in%c("Prolif. T cells","CD4 T cells","Reg T cells", "T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes"))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Lymphocytes')
aggr$celltype <- Idents(aggr)


#Read in NicheNet's ligand-target prior model, ligand-receptor network and weighted integrated networks
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks[[1]] <- weighted_networks[[1]] %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks[[2]] <- weighted_networks[[2]] %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to")) %>% drop_na()
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# Convert the NicheNet network gene symbols from human to mouse based on one-to-one orthology
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

#Overview of ligand expression in each cell population
ligands <- data.frame(ligand=sort(unique(colnames(ligand_target_matrix))))
exp_ligands <- ligands[which(ligands$ligand%in%rownames(aggr)),]
receptors <- data.frame(receptors=sort(unique(lr_network$to)))
Lrg1_network <- sig_network %>% filter(from=="LRG1" | to=="LRG1")
ligand_exp_matrix <- aggr@assays[["RNA"]][match(exp_ligands,rownames(aggr@assays[["RNA"]])),]

# STEP1 : Define expressed genes in sender and receiver population 
#-------------------------------------------------------------------------------
## Define receiver cells
receiver = "Lrg1+ gCap"
expressed_genes_receiver = get_expressed_genes(receiver, aggr, pct = 0.25,assay_oi ="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cells
sender_celltypes <- c( "Prolif. AM","AM1","AM2","AM3","IM","LCM","Fibroblasts","Pericytes","SV EC","Prolif. EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Arterial EC","PV EC","Lymphatic EC","Mesothelial cells","AT1","AT2","Multiciliated cells","Mast cells","Neutrophils","Lymphocytes","Myeloid cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, aggr, 0.25,assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# STEP2 : Define the gene set of interest and a background of gene
#-------------------------------------------------------------------------------
#geneset= pro angiogenic features expressed in Lrg1+ gCap

SUB <-subset(x = aggr, idents = c("aCap","gCap","Lrg1+ aCap","Lrg1+ gCap","SV EC"))
EC_order <- c("Lrg1+ aCap","aCap","Lrg1+ gCap","gCap","SV EC")
SUB$cellstate <- factor(SUB$cellstate,levels = rev(EC_order))
SUB <- SUB %>% NormalizeData() %>% ScaleData()
#Expression of Pro-angiogenic markers
angiogenic.features <- c("Hspg2","Plxnd1","Rgcc","Aplnr","Smad1","Sparc","Col4a1","Col4a2","Flt1","Xbp1","Btg1","Rhoj","Cxcl12","Sox17","Adamts1","Adam15","Vegfa","Ackr3","Lrg1","Nrp1","Kdr","Hspb1","Emp2","Cd34","Serpine1","Unc5b","Hif1a","Tgfbr2")
plot <- DotPlot(SUB, features = angiogenic.features,group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)
plot

gCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ gCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
aCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ aCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
SVEC.angiogenic.features <- plot[["data"]] %>% filter(id=="SV EC"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()

#STEP3 : Define a set of potential ligands 
#-------------------------------------------------------------------------------
#ligands that are expressed by the "sender/niche" cell population and bind a receptor expressed by the "receiver/target" population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expr_rec_gCap <- expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#STEP4 : Perform NicheNet's ligand activity analysis on the gene set of interest 
#-------------------------------------------------------------------------------
#INPUT :
#geneset = gene set of interest
#background_expressed_genes = all expressed genes in receiver cells
#ligand_target_matrix = ligand target matrix denoting regulatory potential scores
#potential_ligands = expressed ligands in sender cells

#OUTPUT : statistical metrics showing how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes
#pearson corr coeff > 0.1 = good candidate
ligand_activities = predict_ligand_activities(geneset = gCap.angiogenic.features, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
write.table(ligand_activities,file = "ligand_activity_Lrg1gCap.txt",row.names = F,sep = "\t")

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(14, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(14, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
Idents(aggr) <- "celltype"
Sub <- subset(aggr,idents = sender_celltypes)
Sub <- Sub %>% NormalizeData()
DotPlot(Sub, features = best_upstream_ligands %>% rev(), cols = c("#FFFFB2","#B10026")) + RotatedAxis()


#Same for Lrg1+ aCap
# STEP1 : Define expressed genes in sender and receiver population 
#-------------------------------------------------------------------------------
## Define receiver cells
receiver = "Lrg1+ aCap"
expressed_genes_receiver = get_expressed_genes(receiver, aggr, pct = 0.25,assay_oi ="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cells
sender_celltypes <- c( "Prolif. AM","AM1","AM2","AM3","IM","LCM","Fibroblasts","Pericytes","SV EC","Prolif. EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Arterial EC","PV EC","Lymphatic EC","Mesothelial cells","AT1","AT2","Multiciliated cells","Mast cells","Neutrophils","Lymphocytes","Myeloid cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, aggr, 0.25,assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# STEP2 : Define the gene set of interest and a background of gene
#-------------------------------------------------------------------------------
#geneset= pro angiogenic features expressed in Lrg1+ gCap

aCap.angiogenic.features

#STEP3 : Define a set of potential ligands 
#-------------------------------------------------------------------------------
#ligands that are expressed by the "sender/niche" cell population and bind a receptor expressed by the "receiver/target" population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expr_rec_aCap <- expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#STEP4 : Perform NicheNet's ligand activity analysis on the gene set of interest 
#-------------------------------------------------------------------------------
#INPUT :
#geneset = gene set of interest
#background_expressed_genes = all expressed genes in receiver cells
#ligand_target_matrix = ligand target matrix denoting regulatory potential scores
#potential_ligands = expressed ligands in sender cells

#OUTPUT : statistical metrics showing how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes
#pearson corr coeff > 0.1 = good candidate
ligand_activities = predict_ligand_activities(geneset = aCap.angiogenic.features, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
write.table(ligand_activities,file = "ligand_activity_Lrg1aCap.txt",row.names = F,sep = "\t")

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(14, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(14, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
Idents(aggr) <- "celltype"
Sub <- subset(aggr,idents = sender_celltypes)
Sub <- Sub %>% NormalizeData()
DotPlot(Sub, features = best_upstream_ligands %>% rev(), cols = c("#FFFFB2","#B10026")) + RotatedAxis()


#Same for SV EC
# STEP1 : Define expressed genes in sender and receiver population 
#-------------------------------------------------------------------------------
## Define receiver cells
receiver = "SV EC"
expressed_genes_receiver = get_expressed_genes(receiver, aggr, pct = 0.25,assay_oi ="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cells
sender_celltypes <- c( "Prolif. AM","AM1","AM2","AM3","IM","LCM","Fibroblasts","Pericytes","SV EC","Prolif. EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Arterial EC","PV EC","Lymphatic EC","Mesothelial cells","AT1","AT2","Multiciliated cells","Mast cells","Neutrophils","Lymphocytes","Myeloid cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, aggr, 0.25,assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# STEP2 : Define the gene set of interest and a background of gene
#-------------------------------------------------------------------------------
#geneset= pro angiogenic features expressed in SV EC

SVEC.angiogenic.features

#STEP3 : Define a set of potential ligands 
#-------------------------------------------------------------------------------
#ligands that are expressed by the "sender/niche" cell population and bind a receptor expressed by the "receiver/target" population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expr_rec_SVEC <- expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#STEP4 : Perform NicheNet's ligand activity analysis on the gene set of interest 
#-------------------------------------------------------------------------------
#INPUT :
#geneset = gene set of interest
#background_expressed_genes = all expressed genes in receiver cells
#ligand_target_matrix = ligand target matrix denoting regulatory potential scores
#potential_ligands = expressed ligands in sender cells

#OUTPUT : statistical metrics showing how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes
#pearson corr coeff > 0.1 = good candidate
ligand_activities = predict_ligand_activities(geneset = SVEC.angiogenic.features, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
write.table(ligand_activities,file = "ligand_activity_SVEC.txt",row.names = F,sep = "\t")

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(14, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(14, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
Idents(aggr) <- "celltype"
Sub <- subset(aggr,idents = sender_celltypes)
Sub <- Sub %>% NormalizeData()
DotPlot(Sub, features = best_upstream_ligands %>% rev(), cols = c("#FFFFB2","#B10026")) + RotatedAxis()


#Add best ligand activities for Lrg1+gCap receiver cells
ligand_activities_Lrg1gCap <- read.table( "ligand_activity_Lrg1gCap.txt",header = T) 
ligand_activities_Lrg1aCap <- read.table( "ligand_activity_Lrg1aCap.txt",header = T) 
ligand_activities_SVEC <- read.table( "ligand_activity_SVEC.txt",header = T) 

best_upstream_ligands <- unique(c(ligand_activities_SVEC %>% top_n(14,pearson) %>% pull(test_ligand),ligand_activities_Lrg1aCap %>% top_n(14,pearson) %>% pull(test_ligand),ligand_activities_Lrg1gCap %>% top_n(14,pearson) %>% pull(test_ligand)))
#

#STEP5 : Infer target genes and receptors of top-ranked ligands and visualize in a heatmap
#-------------------------------------------------------------------------------
#ACTIVE TARGET GENE INFERENCE
#Which genes of the geneset are in the top 200 predicted targets of the top potential ligands
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = angiogenic.features, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
#matrix processing for the heatmap
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.1)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% sort() %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network


#RECEPTORS OF TOP-RANKED LIGANDS
# get the ligand-receptor network of the top-ranked ligands
expressed_receptors <- unique(c(expr_rec_aCap,expr_rec_gCap,expr_rec_SVEC))
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
# convert to a matrix
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")

order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
order_receptors_large <- order_receptors

#RECEPTORS OF TOP-RANKED LIGANDS, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
# get only bona fide ligand-receptor interactions
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

# get the weights of the bona fide ligand-receptor interactions as used in the NicheNet model
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
# convert to a matrix
lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_position="right",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict


#SUMMARY VISUALIZATIONS OF THE NICHENET ANALYSIS
#Figure to combine : 1) ligand activity; 2) ligand expression dotplot; 3) ligand logFC; 4) inference of ligand - target gene activity 
# ligand activity heatmap
vis_ligand_pearson = Reduce(function(x, y, ...) merge(x, y,all = TRUE,by="test_ligand"), 
                            list(ligand_activities_Lrg1gCap  %>% mutate(Lrg1.gCap=pearson )%>% dplyr::select(test_ligand,Lrg1.gCap),
                                 ligand_activities_Lrg1aCap  %>% mutate(Lrg1.aCap=pearson )%>% dplyr::select(test_ligand,Lrg1.aCap),
                                 ligand_activities_SVEC  %>% mutate(SVEC=pearson )%>% dplyr::select(test_ligand,SVEC))) %>% filter(test_ligand %in% best_upstream_ligands) %>% column_to_rownames(.,"test_ligand") %>% as.matrix()
vis_ligand_pearson[is.na(vis_ligand_pearson)] <- min(vis_ligand_pearson[!is.na(vis_ligand_pearson)])
vis_ligand_pearson <- vis_ligand_pearson[match(rev(best_upstream_ligands),rownames(vis_ligand_pearson)),]

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson


# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
aggr$celltype <- factor(aggr$celltype,levels = sender_celltypes)
aggr <- aggr %>% NormalizeData() %>% ScaleData()
rotated_dotplot = DotPlot(aggr %>% subset(celltype %in% sender_celltypes),group.by = "celltype", features = order_ligands_adapted , cols = "RdBu",assay = "RNA") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab("")+xlab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6,  ncol(vis_ligand_target)))
legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank(),axis.title.x = element_text()) 
ggsave("SupplFig3_A_LigandActivity.pdf", width = 3.5, height = 10, units = c("cm"), dpi = 150)
rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("") + xlab("") + scale_y_discrete(position = "right")
ggsave("SupplFig3_A_LigandExpression.pdf", width = 12, height = 10, units = c("cm"), dpi = 150)
p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab("")+xlab("")
ggsave("SupplFig3_A_PredictedTarget.pdf", width = 10, height = 10, units = c("cm"), dpi = 150)

# Checking for best upstream regulators + bona fide receptors + top predicted genes expression value, fraction of expressing cells and differential expression
#------------------------------------------------------------------------------- 

# ligand expression Seurat dotplot
sub <- subset(aggr,idents = c("Lrg1+ gCap","Lrg1+ aCap","SV EC")) 
sub <- sub %>% NormalizeData() %>% ScaleData()

vis_ligand_receptor_network_final <- vis_ligand_receptor_network[,c("Apoe","Hgf")]%>% as.data.frame()
vis_ligand_receptor_network_final <- t(apply(vis_ligand_receptor_network_final, 1, function(x) replace(x, x != max(x, na.rm = TRUE), 0))) %>% as.data.frame()
vis_ligand_receptor_network_final <- vis_ligand_receptor_network_final[c(rownames(vis_ligand_receptor_network_strict),"Scarb1"),]
a <- rbind(vis_ligand_receptor_network_strict,0)
rownames(a) <- c(rownames(vis_ligand_receptor_network_strict),"Scarb1")
vis_ligand_receptor_network_final <- cbind(vis_ligand_receptor_network_final,a)
vis_ligand_receptor_network_final<-vis_ligand_receptor_network_final[match(c("Fgfr1","Fgfr3","Sdc4","Flt1","Kdr","Nrp1","Nrp2","Il1r1","Acvrl1","Tgfbr2","Tgfbr3","Scarb1","Itgb1","Itga5","Bmpr2","Ackr3","Plxnd1","Tnfrsf1a","Ltbr","Amfr","Tnfrsf11b"),rownames(vis_ligand_receptor_network_final)),
                                                                     match(rev(best_upstream_ligands[which(best_upstream_ligands%in%colnames(vis_ligand_receptor_network_final))]),colnames(vis_ligand_receptor_network_final))]


p_ligand_receptor_network_final = vis_ligand_receptor_network_final %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_position="right",legend_title = "Prior interaction potential")
p_ligand_receptor_network_final+ ylab("")+xlab("")+NoLegend()
ggsave("Fig3_H.pdf", width = 12, height = 12, units = c("cm"), dpi = 150)

DotPlot(sub, features = rownames(vis_ligand_receptor_network_final),scale=F, cols = c("white","#730d2b"),dot.scale = 4)  + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))+ theme(legend.position = "right", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("") + xlab("") + scale_y_discrete(position = "left")+ scale_x_discrete(position = "top")+NoLegend()
ggsave("Fig3_H_ReceptorsExpr.pdf", width = 12, height = 4, units = c("cm"), dpi = 150)


#------------------------------------------------------------------------------- 
#-.- FOLLOW-UP ANALYSIS : Circos plot visualization to show active ligand-target links between interacting cells -.-#

best_upstream_ligands <- colnames(vis_ligand_receptor_network_final)
best_upstream_receptors <- rownames(vis_ligand_receptor_network_final)

Mast_cell_specific_ligands = c("Hgf")
Epithelial_specific_ligands = c("Sema3e","Bmp4")
Neutrophil_specific_ligands = c("Il1b")
EC_specific_ligands = c("Cxcl12", "Vegfa","Tnfsf10")
Macrophages_specific_ligands =c("Apoe","Tnf","Adam17","Il1a")
Mesenchymal_specific_ligands =c("Fgf2","Pgf","Col18a1","Tgfb3")
general_ligands = setdiff(best_upstream_ligands,c(EC_specific_ligands,Epithelial_specific_ligands,Macrophages_specific_ligands,Mesenchymal_specific_ligands,Mast_cell_specific_ligands,Neutrophil_specific_ligands)) %>% print()

ligand_type_indication_df = tibble(
  ligand_type = c(rep("EC-specific", times = EC_specific_ligands %>% length()),
                  rep("Mesenchymal-specific", times = Mesenchymal_specific_ligands %>% length()),
                  rep("Epithelial-specific", times = Epithelial_specific_ligands  %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("Macrophages-specific", times = Macrophages_specific_ligands  %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("Mast_cell-specific", times = Mast_cell_specific_ligands  %>% length())),
  ligand = c(EC_specific_ligands,Mesenchymal_specific_ligands,Epithelial_specific_ligands, general_ligands,Macrophages_specific_ligands, Neutrophil_specific_ligands,Mast_cell_specific_ligands))

#-.- Visualize ligand-receptor interactions of the prioritized ligands in a circos plot -.-#

#RECEPTORS OF TOP-RANKED LIGANDS


grid_col_ligand =c("General" = "#109e52","Mesenchymal-specific" = "#7d7a78","Macrophages-specific" = "#1D91C0","EC-specific" = "#eda43e","Neutrophil-specific" = "#4300de","Epithelial-specific" = "#2e2f5c","Mast_cell-specific" = "#9315ed")
grid_col_target =c("geneset" = "tomato")

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% dplyr::rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "receptor") %>% inner_join(ligand_type_indication_df)
grid_col_ligand 
grid_col_receptor =c( "receptor" = "#a37b7a")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
#show only links with a weight higher than a predefined cutoff: links belonging to the 66% of lowest scores were removed.
cutoff_include_all_ligands = circos_links$weight %>% quantile(0.33)
circos_links = circos_links %>% filter(weight>cutoff_include_all_ligands)
links_circle = circos_links %>% dplyr::select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

#Prepare the circos visualization: order ligands and receptors
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(EC_specific_ligands,Mesenchymal_specific_ligands,Epithelial_specific_ligands,general_ligands,Macrophages_specific_ligands,Neutrophil_specific_ligands,Mast_cell_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)

order = c(ligand_order,receptor_order)

#Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "EC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mesenchymal-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Epithelial-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophages-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mast_cell-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)


pdf("Fig3_I.pdf", width = 6, height = 6)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,link.sort = TRUE, order = order, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()
dev.off()


#------------------------------------------------------------------------------------------
# Differentially expressed receptors and key genes of Angiogenic pathways
library(DESeq2)
genelist <- c("Fgfr1","Fgfr3","Sdc4","Flt1","Kdr","Nrp1","Nrp2","Il1r1","Acvrl1","Tgfbr2","Tgfbr3","Scarb1","Itgb1","Itga5","Bmpr2","Ackr3","Plxnd1","Tnfrsf1a","Ltbr","Amfr","Tnfrsf11b")
deg_receptors <- read.xlsx("~/Supplementary_table_3.xlsx") %>% filter(celltype %in% c("gCap","aCap","Venous EC")  & gene%in%genelist & baseMean > 3)

aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap", "aCap","Lrg1+ gCap", "gCap","PV EC") & group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  t(x)})

cts <- cts.split[["gCap"]]
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds_gCap <- DESeq(dds,test="Wald")

cts <- cts.split[["aCap"]]
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds_aCap <- DESeq(dds,test="Wald")

cts <- cts.split[["Venous EC"]]
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds_PVEC <- DESeq(dds,test="Wald")


DEG_gCap <- read.xlsx("~/DEG.pseudobulks/gCap_pseudobulk.xlsx")
DEG_aCap <- read.xlsx("~/DEG.pseudobulks/aCap_pseudobulk.xlsx")
DEG_PVEC <- read.xlsx("~/DEG.pseudobulks/PVEC_pseudobulk.xlsx")


for (gene in unique(deg_receptors$gene)) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('PBS', sample), 'PBS', 'Bleo'))
  cts$condition <- factor(cts$condition,levels = c("PBS","Bleo") )
  p <- ggboxplot(cts, x="condition",y= gene,fill = "condition",ylab = gene, palette = c("PBS"=gg_color_hue(2)[2],"Bleo"=gg_color_hue(2)[1]),outlier.shape = NA)+geom_jitter(shape=16, position=position_jitter(0.2),size=1.5)
  stat.test <- tibble(group1=c("PBS"),group2=c("Bleo"),p.adj=c(round(DEG_gCap$padj[which(DEG_gCap$gene==gene)],digits = 3)),y.position=max(counts(dds_gCap, normalized=TRUE)[gene,])+1)
  p + stat_pvalue_manual(stat.test, size = 2.6,label = "p.adj = {p.adj}")+FontSize(x.text = 10,y.text = 10,x.title = 0,y.title = 11)+NoLegend()
  ggsave(paste0("SupplFig3_B_",gene,"_boxplot.pdf"), width = 3.5, height = 5.8, units = c("cm"), dpi = 200)
}

genes <- c("Mmp14",
           "Eng","Acvrl1","Smad6","Id1",
           "Tgfbr1","Smad1","Smad2","Smad5")
for (gene in genes) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('PBS', sample), 'PBS', 'Bleo'))
  cts$condition <- factor(cts$condition,levels = c("PBS","Bleo") )
  p <- ggboxplot(cts, x="condition",y= gene,fill = "condition",ylab = gene, palette = c("PBS"=gg_color_hue(2)[2],"Bleo"=gg_color_hue(2)[1]),outlier.shape = NA)+geom_jitter(shape=16, position=position_jitter(0.2),size=1.5)
  stat.test <- tibble(group1=c("PBS"),group2=c("Bleo"),p.adj=c(round(DEG_gCap$padj[which(DEG_gCap$gene==gene)],digits = 3)),y.position=max(counts(dds_gCap, normalized=TRUE)[gene,])+1)
  p + stat_pvalue_manual(stat.test, size = 2.6,label = "p.adj = {p.adj}")+FontSize(x.text = 10,y.text = 10,x.title = 0,y.title = 11)+NoLegend()
  ggsave(paste0("SupplFig3_C_",gene,"_boxplot.pdf"), width = 3.5, height = 5.8, units = c("cm"), dpi = 200)
}

