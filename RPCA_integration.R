#author : Marin Truchi

source("/home/truchi/Starter_pack.R")

### ======================================================================================================== ###
#INTEGRATION WITH reciprocal PCA (RPCA)
# goal : samples integration where cells in different biological states are less likely to ‘align’ after integration
# from https://satijalab.org/seurat/articles/integration_rpca.html
#install prerequisite : BiocManager::install("glmGamPoi")
library(SeuratData)

## With SCTransform
#load unnormalized data
aggr <- readRDS("/data/truchi_data/Rdata_objects/before_integration.rds")
# split the dataset into a list of seurat objects 
aggr.list <- SplitObject(aggr, split.by = "sample")
rm(aggr)

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
aggr.list <- lapply(X = aggr.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = aggr.list, nfeatures = 3000)
aggr.list <- PrepSCTIntegration(object.list = aggr.list, anchor.features = features)
aggr.list <- lapply(X = aggr.list, FUN = RunPCA, features = features)

#Perform integration
#k.anchor determines the level of integration. Increasing this parameter  will assist the alignment
aggr <- FindIntegrationAnchors(object.list = aggr.list, normalization.method = "SCT",
                               anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 5)
rm(aggr.list)
# this command creates an 'integrated' data assay
aggr <- IntegrateData(anchorset = aggr, normalization.method = "SCT", dims = 1:50)
DefaultAssay(aggr) <- "integrated"

aggr <- RunPCA(aggr, verbose = FALSE,npcs = 80)
aggr <- RunUMAP(aggr, reduction = "pca", dims = 1:80)
aggr <- FindNeighbors(aggr, reduction = "pca", dims = 1:80)
aggr <- FindClusters(aggr, resolution = 1)

# Visualization
p1 <- DimPlot(object = aggr,group.by = "seurat_clusters", label=TRUE)+ ggtitle("Clusters")+NoAxes()+NoLegend()
p2 <- DimPlot(object = aggr, group.by = "condition")+ ggtitle("Condition")+NoAxes()+NoLegend()
p3 <- DimPlot(object = aggr, group.by = "orig.ident",cols = group.colors)+ ggtitle("Sample")+NoAxes()+NoLegend()
plot_grid(p1,p2,p3,ncol = 3)

DefaultAssay(aggr) <- "RNA"
saveRDS(aggr,file = "/data/truchi_data/Rdata_objects/after_integration.rds")

VlnPlot(aggr,  features = c("percent_dropouts"),pt.size = 0.1,sort = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("nUMIs"),pt.size = 0.1,sort = T,log = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("nGenes"),pt.size = 0.1,sort = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("percent_ribo"),pt.size = 0.1,sort = T,log = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("percent_mito"),pt.size = 0.1,sort = T)+ geom_boxplot()
FeaturePlot(object = aggr, features = c("nUMIs"),ncol = 1,max.cutoff = 5000)

genes <- unlist(findnoisygenes.mm(aggr)[4:5])
genes <- rownames(aggr)[which(rownames(aggr)%ni%genes)]

aggr <- aggr %>% NormalizeData() %>% ScaleData()
a <- FindAllMarkers(aggr,features = genes,logfc.threshold = 0.8,min.pct = 0.8,only.pos = T)
b <- FindMarkers(aggr,ident.1 = 3,min.pct = 0.5,logfc.threshold = 0.5,only.pos = T)
DimPlot(aggr,label = T,cells.highlight = rownames(aggr@meta.data %>% filter(seurat_clusters =="18")))

#Remove unwanted cells and rerun the integration
aggr <- readRDS("/data/truchi_data/Rdata_objects/BLEO_FINAL_RPCA.rds")
aggr <- subset(aggr,idents=c(6,20,22,24,29,30,34,43),invert =T)
aggr[['integrated']] <- NULL
aggr[['SCT']] <- NULL
DefaultAssay(aggr) <- "RNA"
aggr.list <- SplitObject(aggr, split.by = "sample")
aggr.list <- lapply(X = aggr.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = aggr.list, nfeatures = 3000)
aggr.list <- PrepSCTIntegration(object.list = aggr.list, anchor.features = features)
aggr.list <- lapply(X = aggr.list, FUN = RunPCA, features = features)
aggr <- FindIntegrationAnchors(object.list = aggr.list, normalization.method = "SCT",
                               anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 5)
aggr <- IntegrateData(anchorset = aggr, normalization.method = "SCT", dims = 1:50)
DefaultAssay(aggr) <- "integrated"

aggr <- RunPCA(aggr, verbose = FALSE,npcs = 80)
aggr <- RunUMAP(aggr, reduction = "pca", dims = 1:80)
aggr <- FindNeighbors(aggr, reduction = "pca", dims = 1:80)
aggr <- FindClusters(aggr, resolution = 1)

# Visualization
p1 <- DimPlot(object = aggr,group.by = "seurat_clusters", label=TRUE)+ ggtitle("Clusters")+NoAxes()+NoLegend()
p2 <- DimPlot(object = aggr, group.by = "condition")+ ggtitle("Condition")+NoAxes()+NoLegend()
plot_grid(p1,p2,ncol = 2)


DefaultAssay(aggr) <- "RNA"
aggr[['SCT']] <- NULL

VlnPlot(aggr,  features = c("percent_dropouts"),pt.size = 0.1,sort = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("nUMIs"),pt.size = 0.1,sort = T,log = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("nGenes"),pt.size = 0.1,sort = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("percent_ribo"),pt.size = 0.1,sort = T,log = T)+ geom_boxplot()
VlnPlot(aggr,  features = c("percent_mito"),pt.size = 0.1,sort = T)+ geom_boxplot()
FeaturePlot(object = aggr, features = c("Cox4i2"),ncol = 1,max.cutoff = 4)

genes <- unlist(findnoisygenes.mm(aggr)[4:5])
genes <- rownames(aggr)[which(rownames(aggr)%ni%genes)]

aggr <- aggr %>% NormalizeData() %>% ScaleData()
a <- FindAllMarkers(aggr,features = genes,logfc.threshold = 0.8,min.pct = 0.8,only.pos = T)
a <- FindMarkers(aggr,ident.1 = 34,ident.2=c(),min.pct = 0.5,logfc.threshold = 0.5,only.pos = F)
FeaturePlot(object = aggr, features = c("Ccl17"),ncol = 1,max.cutoff = 4)
DimPlot(aggr,label = T,cells.highlight = rownames(aggr@meta.data %>% filter(seurat_clusters =="40")))


#---------------------------------------------------------------------------
#LABELLING

aggr$clusters <- Idents(aggr)
aggr$clusters <- factor(aggr$clusters, levels =c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40"))
Idents(aggr) <- "clusters"

new.cluster.ids <- c("gCap","B cells","NK cells","CD8 T cells","CD4 T cells","NKT1","Neutrophils","Macrophage","cMonocytes","gCap","Macrophage","ncMonocytes","NKT2","Reg T cells","Neutrophils","T helper cells","Macrophage","aCap","Fibroblasts","B cells","Venous EC","B cells","DC","Mast cells","Gamma delta T cells","DC","DC",
                     "Prolif. AM","ILC2","Prolif. EC","AT2","Prolif. T cells","Arterial EC","Platelets","Neutrophil-like Monocytes","cMonocytes","Plasmocytes","Mesothelial cells","Lymphatic EC","LCM","B cells")
names(x = new.cluster.ids) <- levels(x = aggr)
aggr <- RenameIdents(object = aggr, new.cluster.ids)
DimPlot(object = aggr,  reduction = 'umap', label = T, pt.size = 0.5)
aggr$celltype <- Idents(aggr)
DefaultAssay(aggr) <- "integrated"

#-----------------------------------------------------------------
#Subclustering of Fibroblast


SUB <-subset(x = aggr, idents = c("Fibroblasts"))
DimPlot(object = SUB, reduction = 'umap',label = T)
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.1)
DimPlot(object = SUB, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = SUB, reduction = 'umap',group.by = "condition",label=T)+NoLegend()
DimPlot(aggr,label = T,cells.highlight = rownames(SUB@meta.data %>% filter(seurat_clusters =="4")))
DimPlot(object = SUB, reduction = 'umap', group.by = "orig.ident",cols = sample.colors)

DefaultAssay(SUB) <- "RNA" 
SUB <- SUB %>% NormalizeData() %>% ScaleData()
s <- FindAllMarkers(SUB,logfc.threshold = 0.5,min.pct = 0.5,only.pos = T)
FeaturePlot(SUB,features = "Cox4i2",max.cutoff = 4)

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(3,4))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Pericytes')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(5))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'AT2')
DimPlot(object = aggr, reduction = 'umap',label = T)

#-----------------------------------------------------------------
#Subclustering of AT2

SUB <-subset(x = aggr, idents = c("AT2"))
DimPlot(object = SUB, reduction = 'umap',label = T)
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.1)
DimPlot(object = SUB, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = SUB, reduction = 'umap',group.by = "condition",label=T)+NoLegend()
DimPlot(aggr,label = T,cells.highlight = rownames(SUB@meta.data %>% filter(seurat_clusters =="3")))
DimPlot(object = SUB, reduction = 'umap', group.by = "orig.ident",cols = sample.colors)

DefaultAssay(SUB) <- "RNA" 
SUB <- SUB %>% NormalizeData() %>% ScaleData()
s <- FindAllMarkers(SUB,logfc.threshold = 0.5,min.pct = 0.5,only.pos = T)
FeaturePlot(SUB,features = "Sftpc",max.cutoff = 4)

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(2))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'AT1')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(1))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Multiciliated cells')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(3))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Fibroblasts')

DimPlot(object = aggr, reduction = 'umap',label = T)

#-----------------------------------------------------------------
#Subclustering of DC

SUB <-subset(x = aggr, idents = c("DC"))
DimPlot(object = SUB, reduction = 'umap',label = T)
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.1)
DimPlot(object = SUB, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = SUB, reduction = 'umap',group.by = "condition",label=T)+NoLegend()
DimPlot(aggr,label = T,cells.highlight = rownames(SUB@meta.data %>% filter(seurat_clusters =="1")))
DimPlot(object = SUB, reduction = 'umap', group.by = "orig.ident",cols = sample.colors)

DefaultAssay(SUB) <- "RNA" 
SUB <- SUB %>% NormalizeData() %>% ScaleData()
s <- FindAllMarkers(SUB,logfc.threshold = 0.5,min.pct = 0.5,only.pos = T)
FeaturePlot(SUB,features = "F13a1",max.cutoff = 4)

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(0))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'cDC2')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(1))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'cMonocytes')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(2))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'cDC1')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(3))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Mature DC')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(4))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'pDC')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(5))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Prolif. DC')

DimPlot(object = aggr, reduction = 'umap',label = T)


#-----------------------------------------------------------------
#Subclustering of Macrophage

SUB <-subset(x = aggr, idents = c("Macrophage"))
DimPlot(object = SUB, reduction = 'umap',label = T)
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.25)
DimPlot(object = SUB, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = SUB, reduction = 'umap',group.by = "condition",label=T)+NoLegend()
DimPlot(aggr,label = T,cells.highlight = rownames(SUB@meta.data %>% filter(seurat_clusters =="1")))
DimPlot(object = SUB, reduction = 'umap', group.by = "orig.ident",cols = sample.colors)

DefaultAssay(SUB) <- "RNA" 
SUB <- SUB %>% NormalizeData() %>% ScaleData()
s <- FindAllMarkers(SUB,logfc.threshold = 0.5,min.pct = 0.5,only.pos = T)
FeaturePlot(SUB,features = "Spp1",max.cutoff = 4)

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(0,2,5,6))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'AM1')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(1))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'AM2')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(3))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'AM3')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(4))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'IM')

DimPlot(object = aggr, reduction = 'umap',label = T)

#-------------------------------------------------------------------------------------
#Final Dataset
colors.Bleo<- c("Prolif. AM"="#a6d8f5","AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0","IM"="#225EA8","LCM"="#03254f",
                "Prolif. DC"="#ffd1ec","cDC1"="#e670b4","cDC2"="#9c119c","Mature DC"="#960c5a","pDC"="#de439a","cMonocytes"="#66119c","ncMonocytes"="#282b5c","Neutrophils"="#4300de","Neutrophil-like Monocytes"="#38015c","Mast cells"="#9315ed","Platelets"="#de1212",
                "Prolif. T cells"="#F7FCB9","CD4 T cells"="#D9F0A3","Reg T cells"="#ADDD8E","T helper cells"="#78C679","Gamma delta T cells"="#41AB5D","ILC2"="#739154","CD8 T cells"="#006837" ,"NKT1"="#004529","NKT2"="#01736b","NK cells"="#46c798","B cells"="#ebd513","Plasmocytes"="#ccbf47",
                "Prolif. EC"="#fee0d2","gCap"="#911414","aCap"="#752740","Venous EC"="#f2ae5a","Arterial EC"="#e37e12","Lymphatic EC"="#7d3f15","Fibroblasts"="#d3d5db","Pericytes"="#806c6c",
                "Mesothelial cells"="#120c04","AT1"="#8184b8","AT2"="#3e447a","Multiciliated cells"="#1486f7")

celltype_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                    "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                    "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                    "Prolif. EC","gCap","aCap","Venous EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                    "Mesothelial cells","AT1","AT2","Multiciliated cells")
aggr$celltype <- factor(Idents(aggr), levels = celltype_order)
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.5,cols = colors.Bleo,repel = T,group.by = "celltype")
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.1,cols = colors.Bleo,repel = T,ncol = 2,split.by = "orig.ident")
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.5,repel = T,group.by = "condition")
aggr$mice <- paste(substr(aggr$mice,1,nchar(aggr$mice)-5),aggr$condition,sep = "_")
Idents(aggr) <- aggr$celltype

sample_order <- c( "Young_D14.Bleo","Young_D14.PBS","Old_D14.Bleo","Old_D14.PBS","Young_D28.Bleo","Young_D28.PBS","Old_D28.Bleo","Old_D28.PBS","Young_D60.Bleo","Young_D60.PBS","Old_D60.Bleo","Old_D60.PBS" )
aggr$sample <- factor(aggr$sample,levels = sample_order)

age <- group <- c()
for (i in 1:nrow(aggr@meta.data)) {
  age[i] <- strsplit(aggr$orig.ident[i]," ")[[1]][1] 
  group[i] <- ifelse(aggr$condition[i]=="PBS",aggr$condition[i],strsplit(aggr$orig.ident[i]," ")[[1]][2])
}
aggr$age <- age 
aggr$group <- factor(group,levels = c("PBS","D14","D28","D60")) 

saveRDS(aggr,file = "/data/truchi_data/Rdata_objects/BLEO_FINAL_RPCA.rds")

DefaultAssay(aggr) <- "RNA" 
aggr <- aggr %>% NormalizeData() %>% ScaleData()
bleo.markers <- FindAllMarkers(object = aggr, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.8)
a <- read.table("/home/truchi/gene_annot_MOUSE.txt",sep = "\t",header = T) %>% mutate(gene=genes) %>% dplyr::select(-genes)
bleo.markers <- merge(bleo.markers,a,by="gene")
bleo.markers <- bleo.markers %>% dplyr::select(gene,gene_biotype,avg_log2FC,p_val,p_val_adj,pct.1,pct.2,cluster)
bleo.markers$cluster <- factor(bleo.markers$cluster,levels = celltype_order)
a <- bleo.markers %>% group_by(cluster)%>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)
write.xlsx(a,"/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Integrated_Analysis_RPCA_Final/Top_table_BLEO_liger.xlsx")

