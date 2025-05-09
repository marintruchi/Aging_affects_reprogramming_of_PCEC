# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to analyse Visium data and produce the plots of figure 1 and supplementary figure S1
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")
library(STdeconvolve)
library(SpatialExperiment)
library(DropletUtils)

# See STdeconvolve documentation for more details
# https://github.com/JEFworks-Lab/STdeconvolve/issues/37
# https://github.com/JEFworks-Lab/STdeconvolve/blob/devel/docs/combining_datasets.md
# https://jef.works/STdeconvolve/process_bcl_data.html
# https://github.com/JEFworks-Lab/STdeconvolve?tab=readme-ov-file

# Read the data (space ranger outputs have to be organized as described in ?read10xVisium )
S1B1 <- SpatialExperiment::read10xVisium("~/S1B1", type = "sparse", data = "filtered")
colnames(S1B1) <- rownames(S1B1@int_colData@listData[["spatialCoords"]]) <- S1B1@assays@data@listData[["counts"]]@Dimnames[[2]] <- paste0("S1B1.",colnames(S1B1) )
S1C1 <- SpatialExperiment::read10xVisium("~/S1C1", type = "sparse", data = "filtered")
colnames(S1C1) <- rownames(S1C1@int_colData@listData[["spatialCoords"]]) <- S1C1@assays@data@listData[["counts"]]@Dimnames[[2]]<- paste0("S1C1.",colnames(S1C1) )
S1D1 <- SpatialExperiment::read10xVisium("~/S1D1", type = "sparse", data = "filtered")
colnames(S1D1) <- rownames(S1D1@int_colData@listData[["spatialCoords"]]) <- S1D1@assays@data@listData[["counts"]]@Dimnames[[2]]<- paste0("S1D1.",colnames(S1D1) )
S2A1 <- SpatialExperiment::read10xVisium("~/S2A1", type = "sparse", data = "filtered")
colnames(S2A1) <- rownames(S2A1@int_colData@listData[["spatialCoords"]]) <- S2A1@assays@data@listData[["counts"]]@Dimnames[[2]]<- paste0("S2A1.",colnames(S2A1) )
S2C1 <- SpatialExperiment::read10xVisium("~/S2C1", type = "sparse", data = "filtered")
colnames(S2C1) <- rownames(S2C1@int_colData@listData[["spatialCoords"]]) <- S2C1@assays@data@listData[["counts"]]@Dimnames[[2]]<- paste0("S2C1.",colnames(S2C1) )
S2D1 <- SpatialExperiment::read10xVisium("~/S2D1", type = "sparse", data = "filtered")
colnames(S2D1) <- rownames(S2D1@int_colData@listData[["spatialCoords"]]) <- S2D1@assays@data@listData[["counts"]]@Dimnames[[2]]<- paste0("S2D1.",colnames(S2D1) )


# Combine counts and positions for all samples
all_counts <- cbind(as.matrix(counts(S1B1)), as.matrix(counts(S1C1)), as.matrix(counts(S1D1)), as.matrix(counts(S2A1)), as.matrix(counts(S2C1)), as.matrix(counts(S2D1)))
all_pos <- rbind(spatialCoords(S1B1), spatialCoords(S1C1), spatialCoords(S1D1), spatialCoords(S2A1), spatialCoords(S2C1), spatialCoords(S2D1))

# Assign slice identifiers
slices <- c("S1B1","S1C1","S1D1","S2A1","S2C1","S2D1")

S1B1_slice <- rep("S1B1", nrow(spatialCoords(S1B1)))
S1C1_slice <- rep("S1C1", nrow(spatialCoords(S1C1)))
S1D1_slice <- rep("S1D1", nrow(spatialCoords(S1D1)))
S2A1_slice <- rep("S2A1", nrow(spatialCoords(S2A1)))
S2C1_slice <- rep("S2C1", nrow(spatialCoords(S2C1)))
S2D1_slice <- rep("S2D1", nrow(spatialCoords(S2D1)))
all_slice <- c(S1B1_slice, S1C1_slice, S1D1_slice,S2A1_slice,S2C1_slice,S2D1_slice)
names(all_slice) <- rownames(all_pos)

# clean up poor spots and genes
data=data.frame(value=colSums(all_counts))
ggplot(data=data.frame(value=colSums(all_counts)), aes(x=value)) +
  geom_density()+xlim(0,1000)


allClean <- cleanCounts(counts = all_counts,
                        min.reads = 10,
                        min.lib.size = 250,
                        verbose = TRUE,plot =T)

# merge the 6 sections together
all_paths <- list()
all_paths[["S1B1"]] <- as.matrix(t(all_counts[, names(all_slice[all_slice == "S1B1"])]))
all_paths[["S1C1"]] <- as.matrix(t(all_counts[, names(all_slice[all_slice == "S1C1"])]))
all_paths[["S1D1"]] <- as.matrix(t(all_counts[, names(all_slice[all_slice == "S1D1"])]))
all_paths[["S2A1"]] <- as.matrix(t(all_counts[, names(all_slice[all_slice == "S2A1"])]))
all_paths[["S2C1"]] <- as.matrix(t(all_counts[, names(all_slice[all_slice == "S2C1"])]))
all_paths[["S2D1"]] <- as.matrix(t(all_counts[, names(all_slice[all_slice == "S2D1"])]))


# Selecting ODGenes across 6 samples
all_genes <- lapply(all_paths, function(p){
  print(dim(p))
  genes <- colnames(p)
  print(length(genes))
  genes
})
all_genes <- Reduce(intersect, all_genes)
length(all_genes)

# select genes.to.remove with findnoisygenes.mm
findnoisygenes.mm <- function(x){
  mito.genes <- grep(pattern = "^mt-", x = all_genes, value = TRUE)
  ribo.genes <- grep(pattern = "^Rp[sl]", x = all_genes, value = TRUE)
  mitoribo.genes <- grep(pattern = "^Mrp[sl]", x = all_genes, value = TRUE)
  nc.genes <- grep(pattern = "Gm", x = all_genes, value = TRUE)
  nc.genes <- nc.genes[str_detect(nc.genes, pattern = "Gm\\d")]
  Rik.genes <- grep(pattern = "Rik", x = all_genes, value = TRUE)
  Hemoglobin <- c("Hbb-bs", "Hba-a1", "Hba-a2", "Hbb-bt")
  # stress.genes <- c("Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1","Dnaja1","Hsph1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1", "Btg2", "Dusp1")
  return(list("mito.genes"=mito.genes,"mitoribo.genes"=mitoribo.genes,"ribo.genes"=ribo.genes,"nc.genes"=nc.genes,"Rik.genes"=Rik.genes,"Hemoglobin"=Hemoglobin))
}
genes.to.remove <- findnoisygenes.mm(all_genes) %>% unlist()
all_genes <- setdiff(all_genes,genes.to.remove)

all_ODgenes <- lapply(all_paths, function(p) {
  dat <- preprocess(p,
                    extractPos = FALSE,
                    selected.genes = all_genes,
                    nTopGenes = 5,
                    genes.to.remove = NA,
                    removeAbove = 0.95,
                    removeBelow = NA,
                    min.reads = 10,
                    min.lib.size = 1,
                    min.detected = 1,
                    ODgenes = TRUE,
                    od.genes.alpha = 0.05,
                    gam.k = 5,
                    nTopOD = 1000)
  colnames(dat$corpus)
})
unionAllGenes <- Reduce(union, all_ODgenes)
length(unionAllGenes)

allCorpus <- preprocess(t(as.matrix(all_counts)),
                        extractPos = FALSE,
                        selected.genes = unionAllGenes,
                        min.lib.size = 200,
                        min.detected = 10,
                        ODgenes = FALSE  )
allCorpus$pos <- all_pos[rownames(allCorpus$corpus), ]

# filter spots in all_slice based on those kept in the corpus
all_slice_filt <- all_slice[match(rownames(allCorpus$pos),names(all_slice))]
length(all_slice_filt)
allCorpus$slice <- all_slice_filt
dim(allCorpus$corpus)
allCorpus$slm
dim(allCorpus$pos)
length(allCorpus$slice)


# fit LDA models to the data
# very long to run
ldas <- fitLDA(allCorpus$corpus, Ks = seq(5, 20, by = 1),
               perc.rare.thresh = 0.05,
               ncores=20,
               plot=TRUE,
               verbose=TRUE)

optLDAS <- optimalModel(models = ldas, opt = 14)
results <- getBetaTheta(optLDAS, perc.filt = 0.05, betaScale = 1000)
topGenes(results$beta, n=20)

# save the object containing the beta and theta matrices
saveRDS("STdeconvolve.rds")


#Correlation matrix between topics
beta <- as.data.frame(results$beta) %>% t()
max <- data.frame(genes=rownames(beta),max=rowMaxs(beta))
max <- max %>% filter(max>1) %>% pull(genes)
beta <- as.matrix(results$beta)
beta <- beta[,match(max,colnames(beta))]

corMtx <- STdeconvolve::getCorrMtx(m1 = beta,m2 = beta,type = "b")
pairs <- STdeconvolve::lsatPairs(corMtx)

corMtx <- round(corMtx,digits = 3)
corMtx[corMtx<0] <- 0

annotation <- data.frame(row.names = rownames(corMtx),topics=paste0("X",rownames(corMtx)))

# ---- SupplFigS1D Correlation matrix between topics ---- 
pheatmap(as.data.frame(corMtx),annotation_col  = annotation,annotation_row  = annotation,show_colnames = T,show_rownames = T,cluster_cols = T,cluster_rows = T,treeheight_row = 0,legend = F,fontsize_col = 10,fontsize_row = 10,
         color = colorRampPalette(brewer.pal(n = 9, name = "Reds"))(100),annotation_legend = F,angle_col = 0,display_numbers=1,number_color="white")



# Visualize the contribution of genes to each topic
colors <- c("X1"="#E31A1C","X2"="#a1c9ae",  "X3"="#fceadc",  "X4"="#FFFF99",  "X5"="#B2DF8A",  "X6"="#33A02C",  "X7"="grey",  "X8"="#CAB2D6",  "X9"="#A6CEE3",  "X10"="#703210", "X11"="#1F78B4", "X12"="#6A3D9A","X13"="#FDBF6F", "X14"= "#d4583f")

results <- readRDS("STdeconvolve.rds")
deconProp <- results$theta
deconGexp <- results$beta

ps <- lapply(colnames(deconProp), function(celltype) {
  
  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 1))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  
  ## visualize the transcriptional profile
  dat <- data.frame(values = as.vector(log2fc), genes = names(log2fc), order = seq(length(log2fc)))
  # Hide all of the text labels.
  dat$selectedLabels <- ifelse(abs(dat$values)>1 | dat$genes=="Lrg1",dat$genes,"")
  dat$values <- ifelse(dat$values>5 ,5,dat$values)
  dat$text_col <- ifelse(dat$genes=="Lrg1","highlighted","normal")
  dat$hjust <- ifelse(dat$values>0,-0.25,0.8)
  dat$vjust <- ifelse(dat$values>0,-0.8,-0.6)
  dat <- dat %>% filter(values>0 & order<41)
  
  color <- colors[paste0("X",celltype)] %>% as.character()
  plt <- ggplot2::ggplot(data = dat) +
    ggplot2::geom_col(ggplot2::aes(x = order, y = values,
                                   fill = factor(selectedLabels == ""),
                                   color = factor(selectedLabels == "")), width = 1) +
    
    ggplot2::scale_fill_manual(values = c(color,color)) +
    ggplot2::scale_color_manual(values = c(color,color)) +
    
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(-0.8, max(dat$values) + 2.5)) +
    # ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(-2, NA)) +
    
    ggplot2::labs(title = paste0("X", celltype),
                  x = "Gene expression rank",
                  y = "log2(FC)") +
    
    ## placement of gene symbol labels of top genes
    ggplot2::geom_text(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels,color=text_col,vjust = vjust,hjust=hjust), size = 3 ,angle=90) +
    # geom_text_repel(ggplot2::aes(x = order+1, y = values-0.1, label = selectedLabels,color=text_col,vjust = vjust), size = 3,angle=45 ) +
    scale_color_manual(values = c("highlighted"="red", "normal"="black"))+
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size=10, color = "black",angle = 90),
                   axis.text.y = ggplot2::element_text(size=10, color = "black",angle = 90),
                   axis.title.y = ggplot2::element_text(size=10, color = "black"),
                   axis.title.x = ggplot2::element_text(size=10, color = "black",angle = 180),
                   axis.ticks.x = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(size=15),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3, colour = "gray80"),
                   axis.line = ggplot2::element_line(size = 0.5, colour = "black"),
                   legend.position="none")
  plt
})

# ---- SupplFigS1C Top contributing gene for each topic---- 
for (i in c(1,14)) {
  print(ps[[i]])
  ggsave(paste0("topic_",i,"_top_genes.pdf"), width = 14, height = 7, units = c("cm"), dpi = 300)
}

# Save the differentially expressed genes in each topic
ps <- lapply(colnames(deconProp), function(celltype) {
  celltype <- as.numeric(celltype)
  ## highly expressed in cell-type of interest
  highgexp <- names(which(deconGexp[celltype,] > 1))
  ## high log2(fold-change) compared to other deconvolved cell-types
  log2fc <- sort(log2(deconGexp[celltype,highgexp]/colMeans(deconGexp[-celltype,highgexp])), decreasing=TRUE)
  log2fc <- enframe(log2fc,name = "gene",value = "log2FC") # %>% filter(log2fc > 0.4)
  log2fc
})
names(ps) <- paste0("topic.",c(1:14))

# ---- SupplTableS1 Differentially expressed genes in each topic ---- 
write.xlsx(ps,file = "Supplementary_table_S1.xlsx")


# Visualize the repartition of topics per spot in each slice
slices <- c("S1B1","S1C1","S1D1","S2A1","S2C1","S2D1")
viz_topics <- lapply(slices, function(slice) {
  theta1 <- results$theta[rownames(allCorpus$pos)[allCorpus$slice == slice],]
  pos_int1 <- allCorpus$pos[rownames(allCorpus$pos)[allCorpus$slice == slice],] %>% as.data.frame()
  colnames(pos_int1) <- c("y", "x")
  topicCols <- colors[ paste0("X",which(colSums(theta1)>0))]
  
  plot <- vizAllTopics(theta = theta1,
                       pos = pos_int1,
                       r = 12,
                       lwd = 0.01,
                       topicCols =topicCols, 
                       showLegend = TRUE,
                       plotTitle = NA) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1)) +
    ggplot2::geom_rect(data = data.frame(pos_int1),
                       ggplot2::aes(xmin = min(x) - 90, xmax = max(x) + 90,
                                    ymin = min(y) - 90, ymax = max(y) + 90),
                       fill = NA, color = "black", linetype = "solid", size = 0.5) +
    ggplot2::theme(
      plot.background = ggplot2::element_blank()
    ) +
    ggplot2::guides(colour = "none")
  plot
})
names(viz_topics) <- slices

# ---- SupplFigS1B Repartition of topics per spot in each slice---- 
viz_topics[["S2C1"]]


ps <- lapply(colnames(theta1), function(celltype) {
  
  vizTopic(theta = theta1, pos = pos_int1, topic = celltype, plotTitle = paste0("X", celltype),
           size = 0.5, stroke = 1, alpha = 0.5,
           low = "white",
           high = "red") +
    
    ## remove the pixel "Groups", which is the color aesthetic for the pixel borders
    ggplot2::guides(colour = "none")
  
})
gridExtra::grid.arrange(
  grobs = ps,
  layout_matrix = rbind(c(1, 2, 3, 4,5,6,7),
                        c(8,9,10 ,11, 12,13,14))
)



# Plot the mean proportion of each topic across deconvoluted spot per slice
slices <- c("S1B1","S1C1","S1D1","S2A1","S2C1","S2D1")
metadata <- data.frame()
for (slice in slices) {
  theta1 <- results$theta[rownames(allCorpus$pos)[allCorpus$slice == slice],] %>% as.data.frame()
  theta1 <- colMeans(theta1) %>% as.data.frame()
  theta1 <- data.frame(slice=slice,topic=paste0("X",rownames(theta1)),percent=theta1$.)
  metadata <- rbind( metadata,theta1)
}
metadata$slice <- factor(metadata$slice,levels = rev(slices))
metadata$topic <- factor(metadata$topic,levels = rev(c("X1","X14","X13","X4","X5","X6","X2","X3","X10","X9","X11","X7","X8","X12")))

# ---- Fig1F Mean proportion of each topic across deconvoluted spot per slice ---- 
ggplot(metadata, aes(x =percent , y = slice, fill = topic)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = colors) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 270),
        axis.title.y = element_blank(),axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) +NoLegend()



# IPA on differentially expressed genes in each topics

sCP <- c("IL-17A Signaling in Fibroblasts",
         "Pulmonary Fibrosis Idiopathic Signaling Pathway",
         "Wound Healing Signaling Pathway",
         "Signaling by TGF-beta Receptor Complex",
         "Signaling by PDGF",
         "Regulation of Insulin-like Growth Factor (IGF) transport and uptake by IGFBPs",
         "Integrin cell surface interactions",
         "Extracellular matrix organization",
         "Elastic fibre formation",
         "Cell surface interactions at the vascular wall",
         "S100 Family Signaling Pathway" ,
         "Acute Phase Response Signaling",
         "Macrophage Alternative Activation Signaling Pathway",
         "MHC class II antigen presentation",
         "IL-8 Signaling",
         "Antigen Presentation Pathway",
         "Complement cascade",
         "Interferon gamma signaling",
         "Leukocyte Extravasation Signaling",
         "Surfactant metabolism",
         "Actin Cytoskeleton Signaling",
         "NAD Signaling Pathway",
         "Oxidative Phosphorylation",
         "Mitochondrial Fatty Acid Beta-Oxidation",
         "Smooth Muscle Contraction")

topics <- c(1:14)
IPA_topics <- lapply(topics, function(topic) {
  IPA <- read.xlsx("IPA_CP_topics_Visium.xlsx",sheet = topic) %>% filter(CP %in% sCP)
  IPA$zscore <- as.numeric(IPA$zscore)
  IPA$padj <- ifelse(IPA$padj > 10, 10, IPA$padj)
  IPA$padj <- ifelse(IPA$padj < -log10(0.05), NA, IPA$padj)
  IPA$zscore <- ifelse(as.numeric(IPA$zscore) > 4, 4, as.numeric(IPA$zscore))
  IPA$topic <- topic
  IPA
})
IPA_topics <- Reduce(rbind,IPA_topics)
IPA_topics$topic <- factor(as.character(IPA_topics$topic),levels = c(1,14,13,4,2,5,6,3,10,7,9,11,8,12))
IPA_topics$CP <- factor(IPA_topics$CP,levels = rev(sCP))

# ---- Fig1E IPA on differentially expressed genes in each topics ---- 
ggplot(IPA_topics,aes(x=topic,y=CP,size=padj)) + 
  geom_point(aes(col=zscore),alpha=1) +
  scale_colour_gradientn(colours = brewer.pal(n = 9, name = "OrRd"))+
  theme_minimal()+theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.text.x = element_text(size=12))



