# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to prepare data for RNA velocity and produce the plots of figure 6 and supplemental figure S7
# Author : Marin Truchi

# Load required packages and functions
source("prerequisites.R")

library(BUSpaRse)
library(SeuratWrappers)
library(zeallot) # For %<-% that unpacks lists in the Python manner
library(DropletUtils)
library(GGally) # For ggpairs
library(SingleR)
library(scales)
library(plotly)
library(SeuratDisk)

# From : https://bustools.github.io/BUS_notebooks_R/velocity.html
# Run kallisto + bustools for each sample to produce spliced/unspliced matrices 
# ex :
# kb count -i mm_cDNA_introns_index.idx -g tr2g.tsv -x 10xv3 -o kb \
# -c1 cDNA_tx_to_capture.txt -c2 introns_tx_to_capture.txt --workflow lamanno \
# sample_x_L001_R1_001.fastq.gz sample_x_L001_R2_001.fastq.gz \
# sample_x_L002_R1_001.fastq.gz sample_x_L002_R2_001.fastq.gz \
# sample_x_L003_R1_001.fastq.gz sample_x_L003_R2_001.fastq.gz \
# sample_x_L004_R1_001.fastq.gz sample_x_L004_R2_001.fastq.gz


#  For Young D14 Bleomycin sample

# =============== Step 1 : prepare spliced/unspliced matrix =============== #

d <- "~/bleo_YoungD14/output/kb/counts_unfiltered/"
c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = d,spliced_name = "spliced",unspliced_dir = d,unspliced_name = "unspliced")

# How many UMIs are from unspliced transcripts ?
sum(unspliced@x) / (sum(unspliced@x) + sum(spliced@x))

#Matrix dimensions
dim(spliced)
dim(unspliced)

#Most barcodes only have 0 or 1 UMIs detected.
tot_count <- Matrix::colSums(spliced)
summary(tot_count)

#estimate the number of empty droplets
bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)

# Knee plot for filtering empty droplets
# Visualizes the inflection point to filter empty droplets. This function plots 
# different datasets with a different color. Facets can be added after calling this function with `facet_*` functions.

knee_plot <- function(bc_ranks) {
  # purrr pluck shorthand doesn't work on S4Vector DataFrame
  knee_plt <- tibble(rank = map(bc_ranks, ~ .x[["rank"]]), 
                     total = map(bc_ranks, ~ .x[["total"]]),
                     dataset = names(bc_ranks)) %>% 
    unnest(cols = c(rank, total)) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
                  rank_cutoff = map_dbl(bc_ranks, 
                                        ~ max(.x$rank[.x$total >
                                                        metadata(.x)[["inflection"]]])),
                  dataset = names(bc_ranks))
  p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
    geom_line() +
    geom_hline(aes(yintercept = inflection, color = dataset), 
               data = annot, linetype = 2) +
    geom_vline(aes(xintercept = rank_cutoff, color = dataset),
               data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total UMIs")
  return(p)
}
knee_plot(list(spliced = bc_rank, unspliced = bc_uns)) +
  coord_flip()

bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]
# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 0]
sf <- spliced[genes_use, bcs_use]
uf <- unspliced[genes_use, bcs_use]

dim(sf)
dim(uf)

rownames(sf) <- str_remove(rownames(sf), "\\.\\d+")
rownames(uf) <- str_remove(rownames(uf), "\\.\\d+")

#Replace Ensembl ID by Gene symbol
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
genes <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), filters = 'ensembl_gene_id',  values = rownames(sf), mart = ensembl)
sf <- sf[match(genes$ensembl_gene_id,rownames(sf)),]
rownames(sf) <- make.names(genes$external_gene_name,unique=TRUE)
colnames(sf) <- paste0("Young_D14.",colnames(sf))
uf <- uf[match(genes$ensembl_gene_id,rownames(uf)),]
rownames(uf) <- make.names(genes$external_gene_name,unique=TRUE)
colnames(uf) <- paste0("Young_D14.",colnames(uf))

dim(sf)
dim(uf)


# =============== Step 2 : Add metadata from pre-analyzed data =============== #

#Select the subset of interest
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(aggr) <- "cellstate"
sub <- subset(aggr,idents = c("gCap","Lrg1+ gCap","sCap"))
Idents(sub) <- "sample"
YD14 <- subset(sub,idents=c("Young_D14.Bleo","Young_D14.PBS"))
YD14 <- YD14 %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()
YD14 <- RunUMAP(YD14,dims = 1:30)
DimPlot(YD14,label = T,group.by = "cellstate")

# normalize and scale spliced matrice with SCTransform
Seuvelo <- CreateSeuratObject(sf, assay = "spliced") %>% SCTransform(assay = "spliced",variable.features.n = 3000, new.assay.name = "SCT_sf")

YD14$condition <- factor(YD14$condition ,levels = )

YD14$clusters <- factor(YD14$cellstate,levels = c("sCap","Lrg1+ gCap","gCap"))
Idents(YD14) <- "clusters"

cells_to_keep <- YD14@meta.data %>% mutate(barcode=rownames(.))
length(intersect(cells_to_keep$barcode,colnames(Seuvelo)))
cells_to_keep <- cells_to_keep %>% filter(barcode%in%intersect(cells_to_keep$barcode,colnames(Seuvelo)))

# Only keep those cells
Seuvelo <- Seuvelo[, cells_to_keep %>% pull(barcode)]
Seuvelo$clusters <- factor(cells_to_keep$clusters,levels = c("sCap","Lrg1+ gCap","gCap"))
Seuvelo$condition <- factor(cells_to_keep$condition,levels = c("PBS","Bleo"))
# Also only keep relevant cells in the unspliced matrix
uf <- uf[, cells_to_keep %>% pull(barcode)]


# =============== Quality Control =============== #

# normalize and scale unpliced matrice with SCTransform
Seuvelo[["unspliced"]] <- CreateAssayObject(uf)
Seuvelo <- SCTransform(Seuvelo, assay = "unspliced",variable.features.n = 3000, new.assay.name = "SCT_uf")
cols_use <- c("nCount_spliced", "nFeature_spliced", "nCount_unspliced", "nFeature_unspliced")
VlnPlot(Seuvelo, cols_use, pt.size = 0.1, ncol = 2, group.by = "clusters")

# Helper functions for ggpairs
log10_diagonal <- function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + scale_x_log10()
}
log10_points <- function(data, mapping, ...) {
  ggally_points(data, mapping, ...) + scale_x_log10() + scale_y_log10()
}
# How does number of UMI counts in the spliced matrix relate to the number of gene detected in the unspliced matrix?
ggpairs(Seuvelo@meta.data, columns = cols_use,
        upper = list(continuous = "cor"),
        diag = list(continuous = log10_diagonal),
        lower = list(continuous = wrap(log10_points, alpha = 0.1, size=0.3)),
        progress = FALSE)

# =============== UMAP embedding  =============== #

DefaultAssay(Seuvelo) <- "SCT_sf"
Seuvelo <- RunPCA(Seuvelo, verbose = FALSE, npcs = 100)
ElbowPlot(Seuvelo, ndims = 100)
DimPlot(Seuvelo, reduction = "pca",group.by = "clusters", pt.size = 0.5, label = TRUE, repel = TRUE) 

Seuvelo <- RunUMAP(Seuvelo, dims = 1:30) 
DimPlot(Seuvelo, reduction = "umap",group.by = "clusters", pt.size = 0.5, label = TRUE, repel = TRUE) 
DimPlot(Seuvelo, reduction = "umap",group.by = "condition", pt.size = 0.5) 
Idents(Seuvelo) <- "clusters"

#Save seurat object to transfer it in python
Seuvelo[["RNA"]] <- Seuvelo[["spliced"]]
DefaultAssay(Seuvelo) <- "RNA"
SaveH5Seurat(Seuvelo, filename = "Bleo_YD14_CECs_gCap_final.h5Seurat",overwrite = T)
Convert("Bleo_YD14_CECs_gCap_final.h5Seurat", dest = "h5ad",overwrite = T)

# ====================================================================================== # 
#                                                                                        #
#             Run scVelo in python notebook : "scVelo_gCap_Yd14.ipynb"                   #
#                                                                                        #
# ====================================================================================== # 





#  For Old D28 Bleomycin sample

# =============== Step 1 : prepare spliced/unspliced matrix =============== #

d <- "~/bleo_OldD28/output/kb/counts_unfiltered/"
c(spliced, unspliced) %<-% read_velocity_output(spliced_dir = d,spliced_name = "spliced",unspliced_dir = d,unspliced_name = "unspliced")

# How many UMIs are from unspliced transcripts ?
sum(unspliced@x) / (sum(unspliced@x) + sum(spliced@x))

#Matrix dimensions
dim(spliced)
dim(unspliced)

#Most barcodes only have 0 or 1 UMIs detected.
tot_count <- Matrix::colSums(spliced)
summary(tot_count)

#estimate the number of empty droplets
bc_rank <- barcodeRanks(spliced)
bc_uns <- barcodeRanks(unspliced)

# Knee plot for filtering empty droplets
# Visualizes the inflection point to filter empty droplets. This function plots 
# different datasets with a different color. Facets can be added after calling this function with `facet_*` functions.

knee_plot <- function(bc_ranks) {
  # purrr pluck shorthand doesn't work on S4Vector DataFrame
  knee_plt <- tibble(rank = map(bc_ranks, ~ .x[["rank"]]), 
                     total = map(bc_ranks, ~ .x[["total"]]),
                     dataset = names(bc_ranks)) %>% 
    unnest(cols = c(rank, total)) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = map_dbl(bc_ranks, ~ metadata(.x)[["inflection"]]),
                  rank_cutoff = map_dbl(bc_ranks, 
                                        ~ max(.x$rank[.x$total >
                                                        metadata(.x)[["inflection"]]])),
                  dataset = names(bc_ranks))
  p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
    geom_line() +
    geom_hline(aes(yintercept = inflection, color = dataset), 
               data = annot, linetype = 2) +
    geom_vline(aes(xintercept = rank_cutoff, color = dataset),
               data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total UMIs")
  return(p)
}
knee_plot(list(spliced = bc_rank, unspliced = bc_uns)) +
  coord_flip()

bcs_use <- colnames(spliced)[tot_count > metadata(bc_rank)$inflection]
# Remove genes that aren't detected
tot_genes <- Matrix::rowSums(spliced)
genes_use <- rownames(spliced)[tot_genes > 0]
sf <- spliced[genes_use, bcs_use]
uf <- unspliced[genes_use, bcs_use]

dim(sf)
dim(uf)

rownames(sf) <- str_remove(rownames(sf), "\\.\\d+")
rownames(uf) <- str_remove(rownames(uf), "\\.\\d+")

#Replace Ensembl ID by Gene symbol
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
genes <- getBM(attributes=c('external_gene_name','ensembl_gene_id'), filters = 'ensembl_gene_id',  values = rownames(sf), mart = ensembl)
sf <- sf[match(genes$ensembl_gene_id,rownames(sf)),]
rownames(sf) <- make.names(genes$external_gene_name,unique=TRUE)
colnames(sf) <- paste0("Old_D28.",colnames(sf))
uf <- uf[match(genes$ensembl_gene_id,rownames(uf)),]
rownames(uf) <- make.names(genes$external_gene_name,unique=TRUE)
colnames(uf) <- paste0("Old_D28.",colnames(uf))

dim(sf)
dim(uf)


# =============== Step 2 : Add metadata from pre-analyzed data =============== #

#Select the subset of interest
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(aggr) <- "cellstate"
sub <- subset(aggr,idents = c("gCap","Lrg1+ gCap","sCap"))
Idents(sub) <- "sample"
OD28 <- subset(sub,idents=c("Old_D28.Bleo","Old_D28.PBS"))
OD28 <- OD28 %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()
OD28 <- RunUMAP(OD28,dims = 1:30)
DimPlot(OD28,label = T,group.by = "cellstate")

# normalize and scale spliced matrice with SCTransform
Seuvelo <- CreateSeuratObject(sf, assay = "spliced") %>% SCTransform(assay = "spliced",variable.features.n = 3000, new.assay.name = "SCT_sf")

OD28$condition <- factor(OD28$condition ,levels = )

OD28$clusters <- factor(OD28$cellstate,levels = c("sCap","Lrg1+ gCap","gCap"))
Idents(OD28) <- "clusters"

cells_to_keep <- OD28@meta.data %>% mutate(barcode=rownames(.))
length(intersect(cells_to_keep$barcode,colnames(Seuvelo)))
cells_to_keep <- cells_to_keep %>% filter(barcode%in%intersect(cells_to_keep$barcode,colnames(Seuvelo)))

# Only keep those cells
Seuvelo <- Seuvelo[, cells_to_keep %>% pull(barcode)]
Seuvelo$clusters <- factor(cells_to_keep$clusters,levels = c("sCap","Lrg1+ gCap","gCap"))
Seuvelo$condition <- factor(cells_to_keep$condition,levels = c("PBS","Bleo"))
# Also only keep relevant cells in the unspliced matrix
uf <- uf[, cells_to_keep %>% pull(barcode)]


# =============== Quality Control =============== #

# normalize and scale unpliced matrice with SCTransform
Seuvelo[["unspliced"]] <- CreateAssayObject(uf)
Seuvelo <- SCTransform(Seuvelo, assay = "unspliced",variable.features.n = 3000, new.assay.name = "SCT_uf")
cols_use <- c("nCount_spliced", "nFeature_spliced", "nCount_unspliced", "nFeature_unspliced")
VlnPlot(Seuvelo, cols_use, pt.size = 0.1, ncol = 2, group.by = "clusters")

# Helper functions for ggpairs
log10_diagonal <- function(data, mapping, ...) {
  ggally_densityDiag(data, mapping, ...) + scale_x_log10()
}
log10_points <- function(data, mapping, ...) {
  ggally_points(data, mapping, ...) + scale_x_log10() + scale_y_log10()
}
# How does number of UMI counts in the spliced matrix relate to the number of gene detected in the unspliced matrix?
ggpairs(Seuvelo@meta.data, columns = cols_use,
        upper = list(continuous = "cor"),
        diag = list(continuous = log10_diagonal),
        lower = list(continuous = wrap(log10_points, alpha = 0.1, size=0.3)),
        progress = FALSE)

# =============== UMAP embedding  =============== #

DefaultAssay(Seuvelo) <- "SCT_sf"
Seuvelo <- RunPCA(Seuvelo, verbose = FALSE, npcs = 100)
ElbowPlot(Seuvelo, ndims = 100)
DimPlot(Seuvelo, reduction = "pca",group.by = "clusters", pt.size = 0.5, label = TRUE, repel = TRUE) 

Seuvelo <- RunUMAP(Seuvelo, dims = 1:30) 
DimPlot(Seuvelo, reduction = "umap",group.by = "clusters", pt.size = 0.5, label = TRUE, repel = TRUE) 
DimPlot(Seuvelo, reduction = "umap",group.by = "condition", pt.size = 0.5) 
Idents(Seuvelo) <- "clusters"

#Save seurat object to transfer it in python
Seuvelo[["RNA"]] <- Seuvelo[["spliced"]]
DefaultAssay(Seuvelo) <- "RNA"
SaveH5Seurat(Seuvelo, filename = "Bleo_OD28_CECs_gCap_final.h5Seurat",overwrite = T)
Convert("Bleo_OD28_CECs_gCap_final.h5Seurat", dest = "h5ad",overwrite = T)

# ====================================================================================== # 
#                                                                                        #
#             Run scVelo in python notebook : "scVelo_gCap_Od28.ipynb"                   #
#                                                                                        #
# ====================================================================================== # 



## Plot the scVelo outputs 

#Young D14
#Plot genes trends according to latent time in each subpopulation
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
subcolors <- c("sCap"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414")

YD14 <- subset(aggr,subset = cellstate %in% c("gCap","Lrg1+ gCap","sCap") & sample %in% c("Young_D14.Bleo","Young_D14.PBS"))
YD14 <- YD14 %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()
YD14 <- RunUMAP(YD14,dims = 1:30)
DimPlot(YD14,label = T,group.by = "cellstate")

t <- read.table("~/velo_important_genesper_cluster.tsv",sep = "\t",header = T,row.names = 1)
velocyto_genes <- read.table("~/velo_important_genes.tsv",sep = "\t",header = T,row.names = 1) %>% slice_max(fit_likelihood,n=200) %>% dplyr::select(features,fit_likelihood)
velocyto_genes$features <- gsub("\\.","-",velocyto_genes$features)

velo_metadata <- read.table("~/velo_metadata.tsv",sep = "\t",header = T) %>% mutate(barcode=X) %>% dplyr::select(X,barcode,latent_time) %>% column_to_rownames(.,"X")
cells_D14 <- YD14@meta.data %>% dplyr::select(barcode,cellstate) %>% inner_join(velo_metadata) %>% arrange(latent_time) %>% mutate(rank = 1:nrow(.))


# ---- Fig6C density of subpop along trajectory in YD14 ----
ggplot(cells_D14, aes(x=rank, fill=cellstate)) +
  geom_density(alpha=0.75)+scale_fill_manual(values=subcolors)+
  theme_classic()+NoLegend()+NoAxes()


# Define cells initiating and finishing the trajectory
first <- cells_D14 %>% filter(latent_time <= 0.5) %>% pull(barcode)
YD14$half_ltime <- ifelse(YD14$barcode %in% first, "1st_half", "2nd_half" )
Idents(YD14) <- "half_ltime"
common.genes <- intersect(rownames(YD14),velocyto_genes$features)

# Run differential analysis between both groups
mm <- FindMarkers(YD14,ident.1 = "1st_half",ident.2 = "2nd_half",features = common.genes,min.pct = 0,logfc.threshold = 0) %>% mutate(features=rownames(.))
velocyto_genes <- inner_join(velocyto_genes,mm)

# Save the results to run IPA
write.xlsx(velocyto_genes,file = "Velocity_IPA_YD14.xlsx")  


# Plot genes trends according to latent time in each subpopulation 

genes.to.plot <-c("Smad1","Sparc","Col4a1","Ndrg1","Igfbp7",
                  "Ackr3", "Emp1","Cmah","Vwf","Smad7",
                  "Npr3","Hpgd","Clec1a","Bmp6","Itga1")


velo_metadata <- read.table("~/velo_metadata.tsv",sep = "\t",header = T) %>% dplyr::select(X,latent_time) %>% column_to_rownames(.,"X")
genes_counts <- t(GetAssayData(YD14)[genes.to.plot,]) %>% data.frame()
metadata <- merge(YD14@meta.data,genes_counts,by="row.names")%>% column_to_rownames(.,"Row.names")
metadata <- merge(metadata,velo_metadata,by="row.names") %>% column_to_rownames(.,"Row.names")

plist <- list()
for (gene in genes.to.plot) {
  p <- ggplot(metadata, aes_string(x="latent_time", y=gene,  color="cellstate"))+ geom_smooth(method = "gam", aes(fill=cellstate),se=T)+
    scale_color_manual(values=c("sCap"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414"))+
    scale_fill_manual(values = c("sCap"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414")) + theme_classic()+ theme(axis.title.x = element_blank(),axis.text.x = element_blank())
  plist <- list.append(plist,p)
}


#Old D28
#Plot genes trends according to latent time in each subpopulation
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
subcolors <- c("sCap"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414")

OD28 <- subset(aggr,subset = cellstate %in% c("gCap","Lrg1+ gCap","sCap") & sample %in% c("Old_D28.Bleo","Old_D28.PBS"))
OD28 <- OD28 %>% NormalizeData() %>% ScaleData() %>% FindVariableFeatures() %>% RunPCA()
OD28 <- RunUMAP(OD28,dims = 1:30)
DimPlot(OD28,label = T,group.by = "cellstate")

velocyto_genes <- read.table("~/OD28_gCap_velo_important_genes.tsv",sep = "\t",header = T,row.names = 1) %>% slice_max(fit_likelihood,n=200) %>% dplyr::select(features,fit_likelihood)
velocyto_genes$features <- gsub("\\.","-",velocyto_genes$features)

velo_metadata <- read.table("~/OD28_gCap_velo_metadata.tsv",sep = "\t",header = T) %>% mutate(barcode=X) %>% dplyr::select(X,barcode,latent_time) %>% column_to_rownames(.,"X")
cells_D28 <- OD28@meta.data %>% dplyr::select(barcode,cellstate) %>% inner_join(velo_metadata) %>% arrange(latent_time) %>% mutate(rank = 1:nrow(.))

# ---- Fig6D density of subpop along trajectory in OD28 ----
ggplot(cells_D28, aes(x=rank, fill=cellstate)) +
  geom_density(alpha=0.75)+scale_fill_manual(values=subcolors)+
  theme_classic()+NoLegend()+NoAxes()

# Define cells initiating and finishing the trajectory
first <- cells_D28 %>% filter(latent_time <= 0.5) %>% pull(barcode)
OD28$half_ltime <- ifelse(OD28$barcode %in% first, "1st_half", "2nd_half" )
Idents(OD28) <- "half_ltime"
common.genes <- intersect(rownames(OD28),velocyto_genes$features)

# Run differential analysis between both groups
mm <- FindMarkers(OD28,ident.1 = "1st_half",ident.2 = "2nd_half",features = common.genes,min.pct = 0,logfc.threshold = 0) %>% mutate(features=rownames(.))
velocyto_genes <- inner_join(velocyto_genes,mm)

# Save the results to run IPA
write.xlsx(velocyto_genes,file = "Velocity_IPA_OD28.xlsx") 


# ---- SupplFigS7A genes trends according to latent time ----
ggarrange(plotlist=plist,ncol = 5,nrow =3,common.legend = T,legend = "bottom",align = "hv")




# Represent IPA results based on the differential analysis run in both subsets
IPA <- read.xlsx("~/IPA_Velocity.xlsx") %>% arrange(zscore)
IPA$UR <- factor(IPA$UR,levels = IPA$UR)

# ---- SupplFigS7C IPA based on early vs late trajectory signatures ----
ggplot(IPA, aes(x=-log10(pvalue), y=UR, fill=zscore)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.85)+
  geom_vline(xintercept=1.30103,linetype="dashed",color="gray60")+
  scale_fill_gradient2(low = "#024099",mid = "white",high = "#c45b10")+
  theme_minimal()+xlab("-log10(pvalue)")+ylab("enriched pathways")

