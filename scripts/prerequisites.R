# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : List of prerequisite R packages
# Author : Marin Truchi

#----------------------------------------------------------------#
#------------- Load required packages and functions -------------#
#----------------------------------------------------------------#

#---------------- PACKAGES ---------------#

library(Matrix)
library(Seurat)
library(sctransform)
library(tidyverse)
library(ggplot2)
library(BiocGenerics)
library(Matrix.utils)
library(reshape2)
library(cowplot)
library(ggpubr)
library(ggrastr)
library(RColorBrewer)
library(rlist)
library(pheatmap)
library(plyr)
library(dplyr)
library(openxlsx)
library(viridis)  
library(magrittr)
library(matrixStats)
library(data.table)
library(grid)
library(gtable)
library(gridExtra)
library(devtools)
library(ggvenn)
library(enrichR)
library(edgeR)
library(DESeq2)
library(SeuratData)
library(ggrepel)
#--------------- FUNCTIONS ---------------#

#Smooth Norm heatmaps
flexible_normalization <-function(data_in, by_row=TRUE){
  if(by_row){
    row_mean <- apply(data_in,1,mean)
    row_sd   <- apply(data_in,1,sd)
    output <- data_in
    for(i in 1:dim(data_in)[1]){
      output[i,] <- (data_in[i,] - row_mean[i])/row_sd[i]
    }
  }
  #### if by column
  if(!by_row){
    col_mean <- apply(data_in,2,mean)
    col_sd   <- apply(data_in,2,sd)
    output <- data_in
    for(i in 1:dim(data_in)[2]){
      output[,i] <- (data_in[,i] - col_mean[i])/col_sd[i]
    }
  }
  return(output)
}

#GG color Hue
gg_color_hue <-function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#Reverse %in%
'%ni%' <- Negate('%in%')

#Detect noisy/uninterpretable genes in differential analysis
findnoisygenes.mm <- function(x){
  mito.genes <- grep(pattern = "^mt-", x = rownames(x@assays$RNA), value = TRUE)
  ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(x@assays$RNA), value = TRUE)
  mitoribo.genes <- grep(pattern = "^Mrp[sl]", x = rownames(x@assays$RNA), value = TRUE)
  nc.genes <- grep(pattern = "Gm", x = rownames(x@assays$RNA), value = TRUE)
  nc.genes <- nc.genes[str_detect(nc.genes, pattern = "Gm\\d")]
  Rik.genes <- grep(pattern = "Rik", x = rownames(x@assays$RNA), value = TRUE)
  Hemoglobin <- c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt")
  sex.assoc <- c("Xist","Tsix")
  heg <- c("Malat1","Neat1","Lars2")
  stress.genes <- c("Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1","Dnaja1","Hsph1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1", "Btg2", "Dusp1")
  return(list("mito.genes"=mito.genes,"mitoribo.genes"=mitoribo.genes,"ribo.genes"=ribo.genes,"nc.genes"=nc.genes,"Rik.genes"=Rik.genes,"Hemoglobin"=Hemoglobin,"stress.genes"=stress.genes,"sex.assoc"=sex.assoc,"heg"=heg))
}


#Convert Mouse genes to Human genes
convertMouseToHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

convertHumanToMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse , uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

#Code from https://github.com/theislab/2019_Strunz/blob/master/code/Fig1abcdefghi.R
# Define functions for plotting relative frequencies ####

calcFreq <- function(seu, by1, by2) table(as.character(seu@meta.data[,by1]), as.character(seu@meta.data[,by2]))
get_relFreqs <- function(seu, grouping = "batch"){
  freqByGrouping <- as.data.frame.matrix(calcFreq(seu, by1 = 'identifier', by2 = 'celltype'))
  
  relFreqs <- data.frame(freqByGrouping / apply(freqByGrouping, 1, sum), check.names = F)
  relFreqs$id <- rownames(freqByGrouping)
  
  temp <- unique(seu@meta.data[, c("identifier", grouping)])
  map <- temp[, grouping]
  names(map) <- temp$identifier
  relFreqs$group <- map[relFreqs$id]
  
  freqs <- melt(relFreqs, id.vars = c("id", "group"), variable.name = "cluster")
  freqs
}
plot_relFreqs <- function(freqs, clust, cols, order, title = F, save = F, 
                          path = "Plots/rel_frequencies/"){
  freqs$cluster <- paste("Cluster", freqs$cluster)
  freqs$group <- factor(freqs$group, order)
  
  clust = paste("Cluster", clust)
  if(title == F){
    title = clust
  }
  cols <- cols[order]
  p <- ggplot(subset(freqs, cluster %in% clust),aes(x = group, y = value)) + geom_boxplot( fill= cols,outlier.shape = NA)+ 
    geom_jitter(shape=16, position=position_jitter(0.2),size=0.75) + 
    labs(y = "rel. frequency", x = "group", title = title) + guides(fill=FALSE) +
    scale_x_discrete(labels = order)
  
  if(save == T){
    ggsave(paste0(gsub(" ", "_", title), ".pdf"), plot = p, device = "pdf", width = 6, height = 5, units = "in", path = path)
  }
  else{
    plot(p)
  }
}
