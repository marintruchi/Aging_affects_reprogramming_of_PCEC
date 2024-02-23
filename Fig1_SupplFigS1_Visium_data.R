# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to analyse Visium data and produce the plot of Figure 1 and Suppl. Figure S1
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")

#------------------------------------------------------------------------------------------------#
#--------------- Load Spatial count tables and integrate them for spots annotation --------------#
#------------------------------------------------------------------------------------------------#
#Code from :
#https://satijalab.org/seurat/articles/spatial_vignette.html

# ========================= Individual analysis =========================#  
# S1B1 
S1B1 <- Load10X_Spatial(data.dir = "~/S1B1")
S1B1[['sample']] <- "Young.D14_2.down.left"
S1B1[['group']] <- "Young.D14"
S1B1[['replicate']] <- "2.down.left"

mito.genes <- grep(pattern = "^mt-", x = rownames(S1B1@assays$Spatial), value = TRUE)
dropouts <- Matrix::colSums(S1B1@assays$Spatial@data == 0)/nrow(S1B1@assays$Spatial)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(S1B1@assays$Spatial), value = TRUE)
percent.mito <- Matrix::colSums(S1B1@assays$Spatial[mito.genes, ])/Matrix::colSums(S1B1@assays$Spatial)
percent.ribo <- Matrix::colSums(S1B1@assays$Spatial[ribo.genes, ])/Matrix::colSums(S1B1@assays$Spatial)
S1B1[['percent.mito']] <- percent.mito
S1B1[['percent.ribo']] <- percent.ribo
S1B1[['dropouts']] <- dropouts
VlnPlot(S1B1, features = c("nFeature_Spatial", "nCount_Spatial","dropouts","percent.mito","percent.ribo"), ncol=5, cols = "lightsteelblue3")

S1B1 <- SCTransform(S1B1, assay = "Spatial", verbose = FALSE)
S1B1 <- RunPCA(S1B1, assay = "SCT", verbose = FALSE)
S1B1 <- RunUMAP(S1B1, reduction = "pca", dims = 1:50)
S1B1 <- FindNeighbors(S1B1, reduction = "pca", dims = 1:50,k.param = 10)
S1B1 <- FindClusters(S1B1, verbose = FALSE, resolution = 1)

p1 <- DimPlot(S1B1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S1B1, label = TRUE, label.size = 3, pt.size.factor = 2)
plot_grid(p1,p2)


# S1C1
S1C1 <- Load10X_Spatial(data.dir = "data/S1C1")
S1C1[['sample']] <- "Old.D14_1.up.right"
S1C1[['group']] <- "Old.D14"
S1C1[['replicate']] <- "1.up.right"

mito.genes <- grep(pattern = "^mt-", x = rownames(S1C1@assays$Spatial), value = TRUE)
dropouts <- Matrix::colSums(S1C1@assays$Spatial@data == 0)/nrow(S1C1@assays$Spatial)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(S1C1@assays$Spatial), value = TRUE)
percent.mito <- Matrix::colSums(S1C1@assays$Spatial[mito.genes, ])/Matrix::colSums(S1C1@assays$Spatial)
percent.ribo <- Matrix::colSums(S1C1@assays$Spatial[ribo.genes, ])/Matrix::colSums(S1C1@assays$Spatial)
S1C1[['percent.mito']] <- percent.mito
S1C1[['percent.ribo']] <- percent.ribo
S1C1[['dropouts']] <- dropouts
VlnPlot(S1C1, features = c("nFeature_Spatial", "nCount_Spatial","dropouts","percent.mito","percent.ribo"), ncol=5, cols = "lightsteelblue3")

S1C1 <- SCTransform(S1C1, assay = "Spatial", verbose = FALSE)
S1C1 <- RunPCA(S1C1, assay = "SCT", verbose = FALSE)
S1C1 <- RunUMAP(S1C1, reduction = "pca", dims = 1:50)
S1C1 <- FindNeighbors(S1C1, reduction = "pca", dims = 1:50,k.param = 10)
S1C1 <- FindClusters(S1C1, verbose = FALSE, resolution = 1)

p1 <- DimPlot(S1C1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S1C1, label = TRUE, label.size = 3, pt.size.factor = 2)
plot_grid(p1,p2)


# S1D1
S1D1 <- Load10X_Spatial(data.dir = "data/S1D1")
S1D1[['sample']] <- "Old.D14_2.up.left"
S1D1[['group']] <- "Old.D14"
S1D1[['replicate']] <- "2.up.left"

mito.genes <- grep(pattern = "^mt-", x = rownames(S1D1@assays$Spatial), value = TRUE)
dropouts <- Matrix::colSums(S1D1@assays$Spatial@data == 0)/nrow(S1D1@assays$Spatial)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(S1D1@assays$Spatial), value = TRUE)
percent.mito <- Matrix::colSums(S1D1@assays$Spatial[mito.genes, ])/Matrix::colSums(S1D1@assays$Spatial)
percent.ribo <- Matrix::colSums(S1D1@assays$Spatial[ribo.genes, ])/Matrix::colSums(S1D1@assays$Spatial)
S1D1[['percent.mito']] <- percent.mito
S1D1[['percent.ribo']] <- percent.ribo
S1D1[['dropouts']] <- dropouts
VlnPlot(S1D1, features = c("nFeature_Spatial", "nCount_Spatial","dropouts","percent.mito","percent.ribo"), ncol=5, cols = "lightsteelblue3")

S1D1 <- SCTransform(S1D1, assay = "Spatial", verbose = FALSE)
S1D1 <- RunPCA(S1D1, assay = "SCT", verbose = FALSE)
S1D1 <- RunUMAP(S1D1, reduction = "pca", dims = 1:50)
S1D1 <- FindNeighbors(S1D1, reduction = "pca", dims = 1:50,k.param = 10)
S1D1 <- FindClusters(S1D1, verbose = FALSE, resolution = 1)

p1 <- DimPlot(S1D1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S1D1, label = TRUE, label.size = 3, pt.size.factor = 2)
plot_grid(p1,p2)


# S2A1
S2A1 <- Load10X_Spatial(data.dir = "data/S2A1")
S2A1[['sample']] <- "Young.D28_1.up.left"
S2A1[['group']] <- "Young.D28"
S2A1[['replicate']] <- "1.up.left"

mito.genes <- grep(pattern = "^mt-", x = rownames(S2A1@assays$Spatial), value = TRUE)
dropouts <- Matrix::colSums(S2A1@assays$Spatial@data == 0)/nrow(S2A1@assays$Spatial)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(S2A1@assays$Spatial), value = TRUE)
percent.mito <- Matrix::colSums(S2A1@assays$Spatial[mito.genes, ])/Matrix::colSums(S2A1@assays$Spatial)
percent.ribo <- Matrix::colSums(S2A1@assays$Spatial[ribo.genes, ])/Matrix::colSums(S2A1@assays$Spatial)
S2A1[['percent.mito']] <- percent.mito
S2A1[['percent.ribo']] <- percent.ribo
S2A1[['dropouts']] <- dropouts
VlnPlot(S2A1, features = c("nFeature_Spatial", "nCount_Spatial","dropouts","percent.mito","percent.ribo"), ncol=5, cols = "lightsteelblue3")

S2A1 <- SCTransform(S2A1, assay = "Spatial", verbose = FALSE)
S2A1 <- RunPCA(S2A1, assay = "SCT", verbose = FALSE)
S2A1 <- RunUMAP(S2A1, reduction = "pca", dims = 1:50)
S2A1 <- FindNeighbors(S2A1, reduction = "pca", dims = 1:50,k.param = 10)
S2A1 <- FindClusters(S2A1, verbose = FALSE, resolution = 1)

p1 <- DimPlot(S2A1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S2A1, label = TRUE, label.size = 3, pt.size.factor = 2)
plot_grid(p1,p2)


# S2C1
S2C1 <- Load10X_Spatial(data.dir = "data/S2C1")
S2C1[['sample']] <- "Old.D28_1.up.right"
S2C1[['group']] <- "Old.D28"
S2C1[['replicate']] <- "1.up.right"

mito.genes <- grep(pattern = "^mt-", x = rownames(S2C1@assays$Spatial), value = TRUE)
dropouts <- Matrix::colSums(S2C1@assays$Spatial@data == 0)/nrow(S2C1@assays$Spatial)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(S2C1@assays$Spatial), value = TRUE)
percent.mito <- Matrix::colSums(S2C1@assays$Spatial[mito.genes, ])/Matrix::colSums(S2C1@assays$Spatial)
percent.ribo <- Matrix::colSums(S2C1@assays$Spatial[ribo.genes, ])/Matrix::colSums(S2C1@assays$Spatial)
S2C1[['percent.mito']] <- percent.mito
S2C1[['percent.ribo']] <- percent.ribo
S2C1[['dropouts']] <- dropouts
VlnPlot(S2C1, features = c("nFeature_Spatial", "nCount_Spatial","dropouts","percent.mito","percent.ribo"), ncol=5, cols = "lightsteelblue3")

S2C1 <- SCTransform(S2C1, assay = "Spatial", verbose = FALSE)
S2C1 <- RunPCA(S2C1, assay = "SCT", verbose = FALSE)
S2C1 <- RunUMAP(S2C1, reduction = "pca", dims = 1:50)
S2C1 <- FindNeighbors(S2C1, reduction = "pca", dims = 1:50,k.param = 10)
S2C1 <- FindClusters(S2C1, verbose = FALSE, resolution = 1)

p1 <- DimPlot(S2C1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S2C1, label = TRUE, label.size = 3, pt.size.factor = 1)
plot_grid(p1,p2)


# S2D1
S2D1 <- Load10X_Spatial(data.dir = "data/S2D1")
S2D1[['sample']] <- "Old.D28_2.up.left"
S2D1[['group']] <- "Old.D28"
S2D1[['replicate']] <- "2.up.left"

mito.genes <- grep(pattern = "^mt-", x = rownames(S2D1@assays$Spatial), value = TRUE)
dropouts <- Matrix::colSums(S2D1@assays$Spatial@data == 0)/nrow(S2D1@assays$Spatial)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(S2D1@assays$Spatial), value = TRUE)
percent.mito <- Matrix::colSums(S2D1@assays$Spatial[mito.genes, ])/Matrix::colSums(S2D1@assays$Spatial)
percent.ribo <- Matrix::colSums(S2D1@assays$Spatial[ribo.genes, ])/Matrix::colSums(S2D1@assays$Spatial)
S2D1[['percent.mito']] <- percent.mito
S2D1[['percent.ribo']] <- percent.ribo
S2D1[['dropouts']] <- dropouts
VlnPlot(S2D1, features = c("nFeature_Spatial", "nCount_Spatial","dropouts","percent.mito","percent.ribo"), ncol=5, cols = "lightsteelblue3")

S2D1 <- SCTransform(S2D1, assay = "Spatial", verbose = FALSE)
S2D1 <- RunPCA(S2D1, assay = "SCT", verbose = FALSE)
S2D1 <- RunUMAP(S2D1, reduction = "pca", dims = 1:50)
S2D1 <- FindNeighbors(S2D1, reduction = "pca", dims = 1:50,k.param = 10)
S2D1 <- FindClusters(S2D1, verbose = FALSE, resolution = 1)

p1 <- DimPlot(S2D1, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(S2D1, label = TRUE, label.size = 3, pt.size.factor = 1)
plot_grid(p1,p2)

visium.list <- list("S1B1"=S1B1,"S1C1"=S1C1,"S1D1"=S1D1,"S2A1"=S2A1,"S2C1"=S2C1,"S2D1"=S2D1)
visium.list <- lapply(X = visium.list, FUN = SCTransform,assay = "Spatial",variable.features.n = 10000)
merged <- Reduce(merge,visium.list)
VariableFeatures(merged) <- Reduce(intersect,list(visium.list[["S1B1"]]@assays[["SCT"]]@var.features, visium.list[["S1C1"]]@assays[["SCT"]]@var.features,visium.list[["S1D1"]]@assays[["SCT"]]@var.features, visium.list[["S2A1"]]@assays[["SCT"]]@var.features,visium.list[["S2C1"]]@assays[["SCT"]]@var.features,visium.list[["S2D1"]]@assays[["SCT"]]@var.features))
DefaultAssay(merged) <- "SCT"

visium.list <- SplitObject(merged, split.by = "sample")
visium.list <- lapply(X = visium.list, FUN = RunPCA,features = VariableFeatures(merged))

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = visium.list, nfeatures = length(VariableFeatures(merged)))
visium.list <- PrepSCTIntegration(object.list = visium.list, anchor.features = features)

#Perform integration
#k.anchor determines the level of integration. Increasing this parameter  will assist the alignment
merged <- FindIntegrationAnchors(object.list = visium.list, normalization.method = "SCT",
                                 anchor.features = features, dims = 1:50, reduction = "rpca", k.anchor = 5)

# this command creates an 'integrated' data assay
merged <- IntegrateData(anchorset = merged, normalization.method = "SCT", dims = 1:50)
DefaultAssay(merged) <- "integrated"

merged <- RunPCA(merged, verbose = FALSE,npcs = 50)
merged <- RunUMAP(merged, reduction = "pca", dims = 1:50)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:50)
merged <- FindClusters(merged, resolution = 1)
DimPlot(merged, reduction = "umap", group.by = c("ident", "sample"),label = T)

DefaultAssay(merged) <- "Spatial"
merged <- ScaleData(merged)
markers <- FindAllMarkers(merged,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
DoHeatmap(merged,features = markers %>% group_by(cluster) %>% top_n(5,avg_log2FC) %>% pull(gene))
VlnPlot(merged,features = c("Mfap4"),pt.size = 0,sort = T)
images=c("slice1","slice1.1.1","slice1.2.2","slice1.3.3","slice1.4.4","slice1.5.5")
SpatialFeaturePlot(merged, features = c("Jchain"),combine = F,images=images,ncol = 3)
FeaturePlot(merged,features = "Sftpc",max.cutoff = 4)+scale_color_viridis()

#------------------------------------------------------------------------------------------
#LABELLING
merged$clusters <- Idents(merged)
merged$clusters <- factor(merged$clusters, levels =c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))
Idents(merged) <- "clusters"
new.cluster.ids <- c("low counts Alveoli","Alveoli","Bronchioles","Interstitium","Alveoli","Ig+ Alveoli","Fibrosis","Macrophages","Fibrosis","low counts Alveoli","Blood vessel","Interstitium","Ig+ Alveoli","Bronchia","Mesothelial","Blood","Connective tissue","SMC layer","Bronchioles","IFN+ Alveoli")

names(x = new.cluster.ids) <- levels(x = merged)
merged <- RenameIdents(object = merged, new.cluster.ids)
merged$int_tissue <- Idents(merged)
set3 <- brewer.pal(n = 12, name = "Set3")
colors.tissue <- c("Fibrosis"="#de5f5f","Interstitium"="#f79c4d","Macrophages"="#41B6C4","Blood vessel"=set3[12],"Bronchia"=set3[7],"Bronchioles"="#75c26e","Alveoli"="#3e447a","Ig+ Alveoli"="#6e578a","IFN+ Alveoli"="#ed95ea","low counts Alveoli" = "#c7c3e3","Blood"=set3[10],"Mesothelial"="#8c6253","Connective tissue"="black","SMC layer"="#806c6c")
DimPlot(object = merged,  reduction = 'umap', label = F, pt.size = 0.2,group.by = "int_tissue",cols = colors.tissue)+NoLegend()+NoAxes()+ggtitle("")
ggsave("SupplFig1_B.pdf", width = 8, height = 8, units = c("cm"), dpi = 100)

DimPlot(object = merged,  reduction = 'umap', label = F, pt.size = 0.2,group.by = "sample",cols = c("Young.D14_2.down.left"="#d1eced","Old.D14_1.up.right"="#9199db","Old.D14_2.up.left"="#727ddb","Young.D28_1.up.left"="#77a7a8","Old.D28_1.up.right"="#59328c","Old.D28_2.up.left"="#3f1e6b"))+NoLegend()+NoAxes()+ggtitle("")
ggsave("SupplFig1_Bbis.pdf", width = 8, height = 8, units = c("cm"), dpi = 100)

tissue.order <- c("Fibrosis","Interstitium","Macrophages","Blood vessel","Bronchia","Bronchioles","Alveoli","Ig+ Alveoli","IFN+ Alveoli","low counts Alveoli","Blood","SMC layer","Mesothelial","Connective tissue" )
slice.order <- c("Young.D14_2.down.left","Old.D14_1.up.right","Old.D14_2.up.left","Young.D28_1.up.left","Old.D28_1.up.right","Old.D28_2.up.left")
merged$int_tissue <- factor(merged$int_tissue,levels = tissue.order)
merged$sample <- factor(merged$sample,levels = slice.order)
saveRDS(merged,file = "/data/truchi_data/Rdata_objects/BLEO_visium_integrated.rds")
VlnPlot(merged,features=c("nCount_Spatial"),group.by = "int_tissue",pt.size = 0,cols =  colors.tissue,sort = F,y.max=20000)+FontSize(x.text = 8,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.text.x = element_blank())
ggsave("SupplFig1_C.pdf", width = 10, height = 6, units = c("cm"), dpi = 200)

#------------------------------------------------------------------------------------------
#Spatial Dimplots
set3 <- brewer.pal(n = 12, name = "Set3")
SpatialDimPlot(merged, label = F, label.size = 3, pt.size.factor = 1.6,images=c("slice1"),cols = colors.tissue)+NoLegend()
ggsave("Fig1_D_youngD14.pdf", width = 5, height = 5, units = c("cm"), dpi = 300)
SpatialDimPlot(merged, label = F, label.size = 3, pt.size.factor = 1.6,images=c("slice1.3.3"),cols = colors.tissue)+NoLegend()
ggsave("Fig1_D_youngD28.pdf", width = 5, height = 5, units = c("cm"), dpi = 300)

SpatialDimPlot(merged, label = F, label.size = 3, pt.size.factor = 1.6,images=c("slice1.1.1"),cols = colors.tissue)+NoLegend()
ggsave("Fig1_D_oldD14a.pdf", width = 5, height = 5, units = c("cm"), dpi = 300)
SpatialDimPlot(merged, label = F, label.size = 3, pt.size.factor = 1.6,images=c("slice1.4.4"),cols = colors.tissue)+NoLegend()
ggsave("Fig1_D_oldD28a.pdf", width = 5, height = 5, units = c("cm"), dpi = 300)

SpatialDimPlot(merged, label = F, label.size = 3, pt.size.factor = 1.6,images=c("slice1.2.2"),cols = colors.tissue)+NoLegend()
ggsave("Fig1_D_oldD14b.pdf", width = 5, height = 5, units = c("cm"), dpi = 300)
SpatialDimPlot(merged, label = F, label.size = 3, pt.size.factor = 1.6,images=c("slice1.5.5"),cols = colors.tissue)+NoLegend()
ggsave("Fig1_D_oldD28b.pdf", width = 5, height = 5, units = c("cm"), dpi = 300)

#------------------------------------------------------------------------------------------
#Tissue proportions
Idents(merged) <- "int_tissue"
metadata <- merged@meta.data 
metadata = as.data.frame.matrix(table(metadata$sample, metadata$int_tissue))
metadata = as.data.frame(metadata / rowSums(metadata))*100
metadata$sample = rownames(metadata)
metadata$sample =factor(metadata$sample, levels = rev(slice.order))
metadata = gather(metadata, int_tissue, percentage, names(metadata)[1]:names(metadata)[ncol(metadata)-1])
metadata$int_tissue <- factor(metadata$int_tissue, levels =rev(tissue.order))

ggplot(metadata, aes(x =percentage , y = sample, fill = int_tissue)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = colors.tissue) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 270),
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +NoLegend()
ggsave("Fig1_D_tissue_proportions.pdf", width = 10, height = 5, units = c("cm"), dpi = 300)

#------------------------------------------------------------------------------------------
#Heatmap of marker genes 
merged$group.cs <- paste(merged$sample,merged$int_tissue,sep = "-")
Idents(merged) <- "group.cs"
average <- AverageExpression(merged, return.seurat = T)

DefaultAssay(merged) <- "integrated"
merged$int_tissue <- factor(merged$int_tissue,levels =tissue.order )
Idents(merged) <- "int_tissue"
markers <- FindAllMarkers(merged,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
markers$cluster <- factor(markers$cluster,levels=tissue.order)
genes.to.plot <- markers %>% group_by(cluster) %>% filter(gene%ni%c("Hba-a1","Hba-a2","Hbb-bs","Hbb-bt")) %>% slice_max(order_by = avg_log2FC,n =5) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)

# Create Supplementary_table_1 : markers of annotated spots
markers %>% write.xlsx("Supplementary_table_1.xlsx")


data <- GetAssayData(average,assay = "integrated",slot = "data")
data <- data[match(unique(c(genes.to.plot$gene)),rownames(data)),]
annotation <- data.frame(group=colnames(data),sample=gsub("-.*","",colnames(data)),age=gsub("\\..*","",colnames(data)),tissue=gsub(".*-","",colnames(data)))
annotation$sample <- factor(annotation$sample,levels = slice.order)
annotation$tissue <- factor(annotation$tissue,levels = tissue.order )
annotation <-arrange(annotation,tissue,sample)
data <- data[,match(annotation$group,colnames(data))]
annotation <- annotation %>% column_to_rownames(.,"group")
data <- flexible_normalization(data)
data[data < -2] <- -2
data[data > 2] <- 2

pdf("SupplFig1_D.pdf", width = 10, height = 10)
pheatmap(data,annotation_col  = annotation,show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = F,angle_col =90,fontsize_col = 8,legend=T,
         annotation_colors = list(tissue= colors.tissue,age=c("Young"="#dba75e","Old"="#917dd1"),sample=c("Young.D14_2.down.left"="#d1eced","Old.D14_1.up.right"="#9199db","Old.D14_2.up.left"="#727ddb","Young.D28_1.up.left"="#77a7a8","Old.D28_1.up.right"="#59328c","Old.D28_2.up.left"="#3f1e6b")))
dev.off()


merged <- subset(merged,idents=c("low counts Alveoli","Alveoli","Ig+ Alveoli")) 
merged$group.cs <- paste(merged$sample,merged$int_tissue,sep = "-")
Idents(merged) <- "group.cs"

average <- AverageExpression(merged, return.seurat = T)
Idents(merged) <- "int_tissue"
DefaultAssay(merged) <- "integrated"
merged$int_tissue <- factor(merged$int_tissue,levels =tissue.order )
genes.to.plot <- c("Ager","Hopx", "Sftpb", "Sftpd", "Mfap4","Inmt", "Emp2", "Thbd","Igkc","Ighm")
data <- GetAssayData(average,assay = "integrated",slot = "data")

data <- data[match(unique(c(genes.to.plot)),rownames(data)),]
annotation <- data.frame(group=colnames(data),sample=gsub(".*[\\.]([^.]+)[_].*", "\\1", colnames(data)),age=gsub("\\..*","",colnames(data)),tissue=gsub(".*-","",colnames(data)))
annotation$sample <- factor(annotation$sample,levels = c("D14","D28"))
annotation$tissue <- factor(annotation$tissue,levels = tissue.order )
annotation <-arrange(annotation,tissue,sample)
data <- data[,match(annotation$group,colnames(data))]
annotation <- annotation %>% column_to_rownames(.,"group")
data <- flexible_normalization(data)
data[data < -2] <- -2
data[data > 2] <- 2
pheatmap(data,annotation_col  = annotation,show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = F,angle_col =90,fontsize_col = 8,legend=T,gaps_col = c(6,12),
         annotation_colors = list(tissue= colors.tissue,age=c("Young"="#dba75e","Old"="#917dd1"),sample=c("D14"="#de5f5f","D28"="#964141")))



#------------------------------------------------------------------------------------------
#Differential analysis between normal and damaged alveolar areas
DefaultAssay(merged) <- "Spatial"
merged <- NormalizeData(merged)
DEG <- FindMarkers(merged,ident.1 = "Alveoli",ident.2 = "low counts Alveoli",min.pct = 0.25,logfc.threshold = 0)
DEG$gene <- ifelse(rownames(DEG) %in% c("Scgb1a1","Scgb3a1","Scgb3a2","Bpifa1","Cyp2f2","Ager","Hopx","Lpcat1","Sftpa1","Cxcl15","Sftpc","Sftpb","Sfta2","Mfap4","Sema3c","Inmt","Npnt","Cldn5","Thbd","Emp2"), rownames(DEG), "")

DEG$threshold = as.factor(ifelse(DEG$p_val_adj > 0.05, 'not significant', 
                                 ifelse(DEG$p_val_adj < 0.05 & DEG$avg_log2FC <= -0.5, 'down-regulated',
                                        ifelse(DEG$p_val_adj < 0.05 & DEG$avg_log2FC >= 0.5, 'up-regulated', 'not significant'))))

ggplot(DEG, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1.25) +
  scale_x_continuous(limits = c(-1.3, 1)) +
  scale_color_manual(values = c("#c7c3e3", 'grey', "#3e447a")) +
  ggtitle("") +
  xlab("log2 fold change") + ylab("-log10 padj") +
  #geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") +
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(data=filter(DEG, p_val_adj<0.05), aes(label = gene, color = threshold), size = 3.5)
ggsave("SupplFig1_E.pdf", width = 14, height = 14, units = c("cm"), dpi = 100)


#Differential analysis between normal and damaged alveolar areas
DEG <- FindMarkers(merged,ident.1 = "Alveoli",ident.2 = "Ig+ Alveoli",min.pct = 0.25,logfc.threshold = 0)
DEG$gene <- ifelse(rownames(DEG) %in% c("Ighm","Igkc","Jchain","Iglc2","Ighg2b","Iglc1","Igha","Ighg2c","Mzb1","Iglc3","Ager","Hopx","Lpcat1","Sftpa1","Cxcl15","Sftpc","Sftpb","Sfta2","Mfap4","Sema3c","Inmt","Npnt","Cldn5","Thbd","Emp2"), rownames(DEG), "")

DEG$threshold = as.factor(ifelse(DEG$p_val_adj > 0.05, 'not significant', 
                                 ifelse(DEG$p_val_adj < 0.05 & DEG$avg_log2FC <= -0.5, 'down-regulated',
                                        ifelse(DEG$p_val_adj < 0.05 & DEG$avg_log2FC >= 0.5, 'up-regulated', 'not significant'))))

ggplot(DEG, aes(x=avg_log2FC, y=-log10(p_val_adj))) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1.25) +
  scale_x_continuous(limits = c(-4.6, 1.3)) +
  scale_color_manual(values = c("#6e578a", 'grey', "#3e447a")) +
  ggtitle("") +
  xlab("log2 fold change") + ylab("-log10 padj") +
  #geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") +
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(data=filter(DEG, p_val_adj<0.05), aes(label = gene, color = threshold), size = 3.5)
ggsave("SupplFig1_F.pdf", width = 14, height = 14, units = c("cm"), dpi = 100)
