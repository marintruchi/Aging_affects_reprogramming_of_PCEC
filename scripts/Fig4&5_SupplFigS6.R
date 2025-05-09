# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to compare PCEC from old and young mice and produce the plots of figure 4 & 5 and Supplementary figure S6
# Author : Marin Truchi

#Load required packages and functions
source("~/prerequisites.R")


#Loading dataset
#For convenience, load directly the seurat object obtained after the processing and integration steps
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
group_order <- c("PBS","D14","D28","D60")

#Colors for plots
colors.group <- c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c")

# Relative proportion boxplot
# Code from https://github.com/theislab/2019_Strunz/blob/master/code/Fig1abcdefghi.R
# Define functions for plotting relative frequencies

aggr$celltype <-aggr$cellstate
aggr$identifier <- aggr$mice

cellstates <- c("sCap","SV EC","Lrg1+ aCap","Lrg1+ gCap","Prolif. EC","aCap","gCap","PV EC","Arterial EC")
# aggr <- subset(aggr,subset = cellstate %in% cellstates)

Idents(aggr) <- "age"
SUB.old <- subset(aggr,idents="Old")
SUB <- subset(aggr,idents="Young")

plot_relfreq_list <- lapply(cellstates, function(cellstate) {
  y.lim <- max(get_relFreqs(aggr, grouping = "group") %>% mutate(value=round(value*100,digits = 2)) %>% filter(cluster==cellstate) %>% pull(value))
  p <- ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group") %>% mutate(value=round(value*100,digits = 2)), cellstate, colors.group, group_order, title = cellstate, save = F)+ylim(0, y.lim+(0.05*y.lim))+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 90)),
                 plot_relFreqs(get_relFreqs(SUB.old, grouping = "group") %>% mutate(value=round(value*100,digits = 2)), cellstate, colors.group, group_order, title = cellstate, save = F)+ylim(0, y.lim+(0.05*y.lim))+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 90)),ncol=1,nrow=2)  
  p
})
names(plot_relfreq_list) <- cellstates

# ---- Fig4A Relative frequencies of EC between Young and Old ----  
ggarrange(plot_relfreq_list[["sCap"]],plot_relfreq_list[["Lrg1+ gCap"]],plot_relfreq_list[["Lrg1+ aCap"]],plot_relfreq_list[["Prolif. EC"]],plot_relfreq_list[["SV EC"]],plot_relfreq_list[["PV EC"]],plot_relfreq_list[["gCap"]],plot_relfreq_list[["aCap"]],plot_relfreq_list[["Arterial EC"]],nrow = 1,ncol = 9)


# Differential Abundance testing
# Code from https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/Markdowns/10_MultiSplComp.html#differential-abundance-between-conditions

DA_edgeR <- function(dataset,group_pair){
  Idents(dataset) <- "group"
  sub <- subset(x=dataset,idents= group_pair)
  abundances <- table(sub$cellstate, sub$mice) 
  abundances <- unclass(abundances) 
  head(abundances)
  
  # Attaching some column metadata.
  extra.info <- table(sub$mice) %>% as.data.frame()
  colnames(extra.info) <- c("mice","lib.size")
  metadata <- sub@meta.data %>% dplyr::select("mice","sample","age","group") %>% unique()
  extra.info <- inner_join(extra.info,metadata)
  extra.info <- extra.info[match(colnames(abundances),extra.info$mice),]
  
  y.ab <- DGEList(abundances, samples=extra.info)
  y.ab
  
  keep <- filterByExpr(y.ab, group=y.ab$samples$group,min.count = 2, min.total.count = 1, large.n = 1, min.prop = 0.1)
  summary(keep)
  y.ab <- y.ab[keep,]
  summary(keep)
  
  design <- model.matrix(~factor(group), y.ab$samples)
  y.ab <- estimateDisp(y.ab, design, trend.method="none")
  summary(y.ab$common.dispersion)
  
  fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
  summary(fit.ab$var.prior)
  
  res <- glmQLFTest(fit.ab, coef=ncol(design))
  summary(decideTests(res))
  
  table <- topTags(res, n = dim(res)[1]) %>% as.data.frame() %>% mutate(test=paste0(group_pair[1],".vs.",group_pair[2]),cellstate=rownames(.))
  table
}

contrasts <- list(c("PBS","D14"),c("PBS","D28"),c("PBS","D60"))
young_DA_2 <- lapply(X=contrasts,FUN = DA_edgeR,dataset=SUB )
old_DA_2 <- lapply(X=contrasts,FUN = DA_edgeR,dataset=SUB.old )
write.xlsx(list("young"=Reduce(rbind,young_DA_2),"old"=Reduce(rbind,old_DA_2)),file = "PCECs_DA.xlsx")



# BLM-induced signature score based on differentially expressed genes between Bleo and PBS
library(DESeq2)
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr$samples <- aggr$mice
Idents(aggr) <- "RNA"

SUB <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap") & group %ni% c("D60"))

DEG.psdblk.gCap.mice <- read.xlsx("~/Lrg1pos_gCap_pseudobulk.xlsx")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 10)

DefaultAssay(aggr) <- "RNA"
cts <- AggregateExpression(SUB,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

#gCap
cts <- cts.split$`gCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")

SUB <- subset(x = aggr, subset = cellstate %in% c("Lrg1+ gCap","gCap"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.gCap.mice  %>% pull(gene)),name="gCap_bleo")
SUB$id <- paste0(SUB$age,"_",SUB$group)
df <- SUB@meta.data
df <-  data.frame(score=df$gCap_bleo1,df$id,age=df$age,group=df$group)
df$age <- factor(df$age,levels = c("Young","Old"))
gcap <- ggboxplot(df, x="group",y= "score",fill = "age",ylab = "",outlier.shape = NA, palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()

# ---- Fig4B BLM-score in gCap ---- 
gcap


# Wilcoxon test on BLM-score
Idents(SUB) <- "id"
SUB@assays[["RNA"]]@data <- rbind(SUB@assays[["RNA"]]@data,t(data.frame(score=SUB$gCap_bleo1)))
DefaultAssay(SUB) <- "RNA"

bleo_gcap <- list(FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0) %>% as.data.frame() %>% mutate(contrast="Young_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D14",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D28",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D28"))
bleo_gcap <- Reduce(rbind,bleo_gcap)
write.xlsx(bleo_gcap,"score_gCap.xlsx")


#Supplementary heatmap of bleomycin-induced markers in gCap 
DEG.psdblk.gCap.mice <- read.xlsx("~/Supplementary_table_S3.xlsx",sheet="gCap")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% top_n(50,-log10(padj))
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.gCap.mice$gene,]  %>% as.data.frame()  
names <- gsub("gCap_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

# ---- SupplFigS6A Heatmap of genes modulated by BLM in gCap ---- 
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)


#aCap
cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.aCap.mice <- read.xlsx("~/Supplementary_table_S3.xlsx",sheet="aCap")
DEG.psdblk.aCap.mice <- DEG.psdblk.aCap.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 5)

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SUB <- subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","aCap"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.aCap.mice  %>% pull(gene)),name="aCap_bleo")
SUB$id <- paste0(SUB$age,"_",SUB$group)
df <- SUB@meta.data
df <-  data.frame(score=df$aCap_bleo1,df$id,age=df$age,group=df$group)
df$age <- factor(df$age,levels = c("Young","Old"))
acap <- ggboxplot(df, x="group",y= "score",fill = "age",ylab = "",outlier.shape = NA, palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()

# ---- Fig4B BLM-score in aCap ---- 
acap

# Wilcoxon test on BLM-score
Idents(SUB) <- "id"
SUB@assays[["RNA"]]@data <- rbind(SUB@assays[["RNA"]]@data,t(data.frame(score=SUB$aCap_bleo1)))
DefaultAssay(SUB) <- "RNA"

bleo_acap <- list(FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0) %>% as.data.frame() %>% mutate(contrast="Young_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D14",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D28",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D28"))
bleo_acap <- Reduce(rbind,bleo_acap)
write.xlsx(bleo_acap,"score_aCap.xlsx")


# Supplementary heatmap of bleomycin-induced markers in aCap 
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.aCap.mice$gene,]  %>% as.data.frame()  
names <- gsub("aCap_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

# ---- SupplFigS6B Heatmap of genes modulated by BLM in aCap---- 
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)



#For sCap (merge sCap and PVEC from PBS to have enough cells for pseudobulks)
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
sCap <-subset(x = aggr, subset = cellstate %in% c("sCap") & group %ni% c("PBS"))
gCap <-subset(x = aggr, subset = cellstate %in% c("gCap") & group %in% c("PBS"))
aggr <- merge(sCap,gCap)
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

DEG.psdblk.sCap.mice <- read.xlsx("~/Supplementary_table_S3.xlsx",sheet="sCap")
DEG.psdblk.sCap.mice <- DEG.psdblk.sCap.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 5)

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SUB <- subset(x = aggr, subset = cellstate %in% c("sCap"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.sCap.mice  %>% pull(gene)),name="sCap_bleo")
SUB$id <- paste0(SUB$age,"_",SUB$group)
df <- SUB@meta.data
df <-  data.frame(score=df$sCap_bleo1,df$id,age=df$age,group=df$group)
df$age <- factor(df$age,levels = c("Young","Old"))
sCap <- ggboxplot(df, x="group",y= "score",fill = "age",ylab = "",outlier.shape = NA, palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()

# ---- SupplFigS6C BLM-score in sCap ---- 
sCap

# Wilcoxon test on BLM-score
Idents(SUB) <- "id"
SUB@assays[["RNA"]]@data <- rbind(SUB@assays[["RNA"]]@data,t(data.frame(score=SUB$sCap_bleo1)))
DefaultAssay(SUB) <- "RNA"

bleo_sCap <- list(FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0) %>% as.data.frame() %>% mutate(contrast="Young_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D14",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D28",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D28"))
bleo_sCap <- Reduce(rbind,bleo_sCap)
write.xlsx(bleo_sCap,"score_sCap.xlsx")


# Supplementary heatmap of bleomycin-induced markers in sCap 
DEG.psdblk.sCap.mice <- DEG.psdblk.sCap.mice %>% top_n(50,-log10(padj))
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.sCap.mice$gene,]  %>% as.data.frame()  
names <- gsub("sCap_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

# ---- SupplFigS6D Heatmap of genes modulated by BLM in sCap ----
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)



# For SVEC (merge SVEC and PVEC from PBS to have enough cells for pseudobulks)
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SVEC <-subset(x = aggr, subset = cellstate %in% c("SVEC") & group %ni% c("PBS"))
PVEC <-subset(x = aggr, subset = cellstate %in% c("PVEC") & group %in% c("PBS"))
aggr <- merge(SVEC,PVEC)
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

DEG.psdblk.SVEC.mice <- read.xlsx("~/Supplementary_table_S3.xlsx",sheet="SV EC")
DEG.psdblk.SVEC.mice <- DEG.psdblk.SVEC.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 1)

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SUB <- subset(x = aggr, subset = cellstate %in% c("SV EC"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.SVEC.mice  %>% pull(gene)),name="SVEC_bleo")
SUB$id <- paste0(SUB$age,"_",SUB$group)
df <- SUB@meta.data
df <-  data.frame(score=df$SVEC_bleo1,df$id,age=df$age,group=df$group)
df$age <- factor(df$age,levels = c("Young","Old"))
SVEC <- ggboxplot(df, x="group",y= "score",fill = "age",ylab = "",outlier.shape = NA, palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()

# ---- SupplFigS6E BLM-score in SV EC ----
SVEC


# Wilcoxon test on BLM-score
Idents(SUB) <- "id"
SUB@assays[["RNA"]]@data <- rbind(SUB@assays[["RNA"]]@data,t(data.frame(score=SUB$SVEC_bleo1)))
DefaultAssay(SUB) <- "RNA"

bleo_SVEC <- list(FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0) %>% as.data.frame() %>% mutate(contrast="Young_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D14",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D28",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D28"))
bleo_SVEC <- Reduce(rbind,bleo_SVEC)
write.xlsx(bleo_SVEC,"score_SVEC.xlsx")


# Supplementary heatmap of bleomycin-induced markers in SVEC 
DEG.psdblk.SVEC.mice <- DEG.psdblk.SVEC.mice %>% filter(gene %ni% unlist(findnoisygenes.mm(aggr))) %>% top_n(50,-log10(padj))
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.SVEC.mice$gene,]  %>% as.data.frame()  
names <- gsub("SVEC_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

# ---- SupplFigS6F Heatmap of genes modulated by BLM in sCap ----
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)


#For PVEC (merge PVEC and PVEC from PBS to have enough cells for pseudobulks)
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("PV EC")& group %ni% c("D60"))
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

DEG.psdblk.PVEC.mice <- read.xlsx("~/Supplementary_table_S3.xlsx",sheet="PV EC")
DEG.psdblk.PVEC.mice <- DEG.psdblk.PVEC.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 5)

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SUB <- subset(x = aggr, subset = cellstate %in% c("PV EC"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.PVEC.mice  %>% pull(gene)),name="PVEC_bleo")
SUB$id <- paste0(SUB$age,"_",SUB$group)
df <- SUB@meta.data
df <-  data.frame(score=df$PVEC_bleo1,df$id,age=df$age,group=df$group)
df$age <- factor(df$age,levels = c("Young","Old"))
PVEC <- ggboxplot(df, x="group",y= "score",fill = "age",ylab = "",outlier.shape = NA, palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()

# ---- SupplFigS6G BLM-score in PV EC ----
PVEC


# Wilcoxon test on BLM-score
Idents(SUB) <- "id"
SUB@assays[["RNA"]]@data <- rbind(SUB@assays[["RNA"]]@data,t(data.frame(score=SUB$PVEC_bleo1)))
DefaultAssay(SUB) <- "RNA"

bleo_PVEC <- list(FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0) %>% as.data.frame() %>% mutate(contrast="Young_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Young_PBS",ident.2 = "Young_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Young_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD28"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_PBS",ident.2 = "Old_D60",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="Old_PBSvsD60"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D14",ident.2 = "Young_D14",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D14"),
                  FindMarkers(SUB,slot = "data",ident.1 = "Old_D28",ident.2 = "Young_D28",features = "score",min.pct = 0,logfc.threshold = 0)%>% as.data.frame() %>% mutate(contrast="OldvsYoung_D28"))
bleo_PVEC <- Reduce(rbind,bleo_PVEC)
write.xlsx(bleo_PVEC,"score_PVEC.xlsx")


# Supplementary heatmap of bleomycin-induced markers in PVEC 
DEG.psdblk.PVEC.mice <- DEG.psdblk.PVEC.mice %>% filter(gene %ni% unlist(findnoisygenes.mm(aggr))) %>% top_n(50,-log10(padj))
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.PVEC.mice$gene,]  %>% as.data.frame()  
names <- gsub("PVEC_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

# ---- SupplFigS6G BLM-score in PV EC ----
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)


# Table of BLM-induced score per cell
table <- write.xlsx(list("gcap"=gcap[["data"]],"acap"=acap[["data"]],"sCap"=sCap[["data"]],"SVEC"=SVEC[["data"]],"PVEC"=PVEC[["data"]]),file = "BLM_induced_score_per_cell.xlsx")



#Differential expression analysis between old and young aCap or gCap cells pseudobulks using DESeq2 in fibrotic condition 

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap","PV EC") & group %ni% c("PBS","D60"))
aggr$samples <- aggr$mice

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
  contam_t <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/FINAL_Figures/Supplementary_table_S2.xlsx") %>% filter(avg_log2FC > 4 & cluster %ni%  c("aCap","gCap","PV EC","Arterial EC","Prolif. EC")) %>% pull(gene)
  # stress.genes <- c("Fosb", "Fos", "Jun", "Junb", "Jund", "Atf3", "Egr1","Dnaja1","Hsph1", "Hspa1a", "Hspa1b", "Hsp90ab1", "Hspa8", "Hspb1", "Ier3", "Ier2", "Btg1", "Btg2", "Dusp1")
  return(list("mito.genes"=mito.genes,"mitoribo.genes"=mitoribo.genes,"ribo.genes"=ribo.genes,"nc.genes"=nc.genes,"Rik.genes"=Rik.genes,"Hemoglobin"=Hemoglobin,"sex.assoc"=sex.assoc,"heg"=heg,"contam_t"=contam_t))
}
unwanted_genes <- unlist(findnoisygenes.mm(aggr))

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

cts <- cts.split$`gCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
dds_gCap <- dds
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="Old_BLM_vs_Young_BLM")

counts <- data.frame(gene = rownames(cts),raw_aggr_counts_mean_Young=cts[,grepl("Young",colnames(cts))] %>% rowMeans(),
                     raw_aggr_counts_mean_Old=cts[,grepl("Old",colnames(cts))] %>% rowMeans())
DEG.psdblk.gCap.mice <- inner_join(DEG.psdblk.gCap.mice,counts) %>% filter(gene %ni% unwanted_genes) %>% filter(raw_aggr_counts_mean_Old > 30 | raw_aggr_counts_mean_Young > 30) %>% 
  dplyr::select("test","gene","log2FoldChange","pvalue","padj","raw_aggr_counts_mean_Old","raw_aggr_counts_mean_Young") %>% arrange(padj)


cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="Old_BLM_vs_Young_BLM")

counts <- data.frame(gene = rownames(cts),raw_aggr_counts_mean_Young=cts[,grepl("Young",colnames(cts))] %>% rowMeans(),
                     raw_aggr_counts_mean_Old=cts[,grepl("Old",colnames(cts))] %>% rowMeans())
DEG.psdblk.aCap.mice <- inner_join(DEG.psdblk.aCap.mice,counts) %>% filter(gene %ni% unwanted_genes) %>% filter(raw_aggr_counts_mean_Old > 20 | raw_aggr_counts_mean_Young > 20) %>% 
  dplyr::select("test","gene","log2FoldChange","pvalue","padj","raw_aggr_counts_mean_Old","raw_aggr_counts_mean_Young") %>% arrange(padj)


# ---- SupplTableS5 ---- 
# Create Supplementary_table_5 : Differential expression analysis (DESeq2) between old and young aCap or gCap cells in fibrotic condition
write.xlsx(list("gCap"=DEG.psdblk.gCap.mice,"aCap"=DEG.psdblk.aCap.mice),file = "Supplementary_table_S5.xlsx")

gene_list <- read.xlsx("~/Supplementary_table_S5.xlsx",sheet="gCap")

upgene_Bleo <- gene_list %>% arrange(dplyr::desc(log2FoldChange)) %>% filter(log2FoldChange > 0.25 & padj<0.05)  %>% na.omit() %>% pull(gene)
dwgene_Bleo <- arrange(gene_list, log2FoldChange) %>% filter(log2FoldChange < -0.25 & padj<0.05)  %>%  na.omit() %>% pull(gene)
genestoshow <- c(upgene_Bleo, dwgene_Bleo)


print(genestoshow)
gene_list$gene <- ifelse(gene_list$gene %in% genestoshow, gene_list$gene, "")
gene_list$threshold = as.factor(ifelse(gene_list$padj > 0.05, 'not significant', 
                                       ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange <= -0.25, 'down-regulated',
                                              ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange >= 0.25, 'up-regulated', 'not significant'))))

# ---- Fig5A VolcanoPlot Old vs Young gCap in Fibrosis ---- 
ggplot(gene_list, aes(x=log2FoldChange, y=-log10(padj))) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1.25) +
  # scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#dba75e', 'grey', '#917dd1')) +
  ggtitle("") +
  xlab("log2 fold change") + ylab("-log10 padj") +
  #geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") +
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(data=filter(gene_list, padj<0.05), aes(label = gene, color = threshold), size = 3.5)


DEG_gCap <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Last_figs/Supplementary_table_S5.xlsx",sheet="gCap")
genes <- c("Cd74","H2-Ab1","Gbp4","Klf2","Klf10","Peg3")

# ---- Fig5B Boxplot ---- 
for (gene in genes) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('Young', sample), 'Young', 'Old')) %>% mutate(time.point=ifelse(grepl('D14', sample), 'D14', 'D28'))
  cts$condition <- factor(cts$condition,levels = c("Young","Old") )
  p <- ggboxplot(cts, x="condition",y= gene,fill = "condition",ylab = gene, palette = c("Young"='#dba75e',"Old"='#917dd1'),outlier.shape = NA)+geom_jitter(aes(shape = time.point), position=position_jitter(0.2),size=1.5)
  ypos <- max(counts(dds_gCap, normalized=TRUE)[gene,])
  stat.test <- tibble(group1=c("Young"),group2=c("Old"),p.adj=c(round(DEG_gCap$padj[which(DEG_gCap$gene==gene)],digits = 5)),y.position=ypos+c(ypos*0.05))
  p + ylim(c(0,ypos+c(ypos*0.1)))+ stat_pvalue_manual(stat.test, size = 2,label = "p.adj = {p.adj}")+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 9)+NoLegend()
}

# Corresponding table 
df.list <- lapply(genes, function(gene) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('Young', sample), 'Young', 'Old')) %>% mutate(time.point=ifelse(grepl('D14', sample), 'D14', 'D28'))
  colnames(cts)[1] <- "Value"
  cts$gene <- gene
  cts
})



# Differential expression analysis between old and young aCap or gCap cells pseudobulks using DESeq2 in physiological condition 

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap","PV EC") & group %in% c("PBS"))
aggr$samples <- aggr$mice

unwanted_genes <- unlist(findnoisygenes.mm(aggr))

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

cts <- cts.split$`gCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
dds_gCap <- dds
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="Old_PBS_vs_Young_PBS")

counts <- data.frame(gene = rownames(cts),raw_aggr_counts_mean_Young=cts[,grepl("Young",colnames(cts))] %>% rowMeans(),
                     raw_aggr_counts_mean_Old=cts[,grepl("Old",colnames(cts))] %>% rowMeans())
DEG.psdblk.gCap.mice <- inner_join(DEG.psdblk.gCap.mice,counts) %>% filter(gene %ni% unwanted_genes) %>% filter(raw_aggr_counts_mean_Old > 5 | raw_aggr_counts_mean_Young > 5) %>% 
  dplyr::select("test","gene","log2FoldChange","pvalue","padj","raw_aggr_counts_mean_Old","raw_aggr_counts_mean_Young") %>% arrange(padj)


cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="Old_PBS_vs_Young_PBS")

counts <- data.frame(gene = rownames(cts),raw_aggr_counts_mean_Young=cts[,grepl("Young",colnames(cts))] %>% rowMeans(),
                     raw_aggr_counts_mean_Old=cts[,grepl("Old",colnames(cts))] %>% rowMeans())
DEG.psdblk.aCap.mice <- inner_join(DEG.psdblk.aCap.mice,counts) %>% filter(gene %ni% unwanted_genes) %>% filter(raw_aggr_counts_mean_Old > 5 | raw_aggr_counts_mean_Young > 5) %>% 
  dplyr::select("test","gene","log2FoldChange","pvalue","padj","raw_aggr_counts_mean_Old","raw_aggr_counts_mean_Young") %>% arrange(padj)


# ---- SupplTableS6 ---- 
# Create Supplementary_table_6 : Differential expression analysis (DESeq2) between old and young aCap or gCap cells in physiological condition
write.xlsx(list("gCap"=DEG.psdblk.gCap.mice,"aCap"=DEG.psdblk.aCap.mice,"PVEC"=DEG.psdblk.PVEC.mice),file = "Supplementary_table_S6.xlsx")


gene_list <- read.xlsx("~/Supplementary_table_S6.xlsx",sheet="gCap")

upgene_Bleo <- gene_list %>% arrange(dplyr::desc(log2FoldChange)) %>% filter(log2FoldChange > 0.25 & padj<0.05)  %>% na.omit() %>% pull(gene)
dwgene_Bleo <- arrange(gene_list, log2FoldChange) %>% filter(log2FoldChange < -0.25 & padj<0.05)  %>%  na.omit() %>% pull(gene)
genestoshow <- c(upgene_Bleo, dwgene_Bleo)
print(genestoshow)

gene_list$gene <- ifelse(gene_list$gene %in% genestoshow, gene_list$gene, "")

gene_list$threshold = as.factor(ifelse(gene_list$padj > 0.05, 'not significant', 
                                       ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange <= -0.25, 'down-regulated',
                                              ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange >= 0.25, 'up-regulated', 'not significant'))))


# ---- Fig5D VolcanoPlot Old vs Young gCap in Fibrosis ---- 
ggplot(gene_list, aes(x=log2FoldChange, y=-log10(padj))) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1) +
  # scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#dba75e', 'grey', '#917dd1')) +
  ggtitle("") +
  xlab("log2 fold change") + ylab("-log10 padj") +
  #geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") +
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(data=filter(gene_list, padj<0.05), aes(label = gene, color = threshold), size = 3.5)+FontSize(x.text =9 ,x.title = 9,y.text = 9,y.title = 9 )
ggsave("Fig5D.pdf", width = 15, height = 12, units = c("cm"), dpi = 200)


DEG_gCap <- read.xlsx("~/Supplementary_table_S6.xlsx",sheet="gCap")
genes <- c("Lrg1","Col15a1","Ntrk2","Aplnr","Vwf","Slc6a2","Prss23","Plat")

# ---- Fig5E Boxplot ---- 
for (gene in genes) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('Young', sample), 'Young', 'Old')) %>% mutate(time.point=ifelse(grepl('D14', sample), 'D14', ifelse(grepl('D28', sample), 'D28', 'D60')))
  cts$condition <- factor(cts$condition,levels = c("Young","Old") )
  p <- ggboxplot(cts, x="condition",y= gene,fill = "condition",ylab = gene, palette = c("Young"='#dba75e',"Old"='#917dd1'),outlier.shape = NA)+geom_jitter(aes(shape = time.point), position=position_jitter(0.2),size=1.5)
  ypos <- max(counts(dds_gCap, normalized=TRUE)[gene,])
  stat.test <- tibble(group1=c("Young"),group2=c("Old"),p.adj=c(round(DEG_gCap$padj[which(DEG_gCap$gene==gene)],digits = 5)),y.position=ypos+c(ypos*0.05))
  p + ylim(c(0,ypos+c(ypos*0.1)))+ stat_pvalue_manual(stat.test, size = 2,label = "p.adj = {p.adj}")+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 9)+NoLegend()
  ggsave(paste0("Fig5_E_",gene,"_boxplot.pdf"), width = 4.5, height = 4.5, units = c("cm"), dpi = 200)
}

df.list <- lapply(genes, function(gene) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('Young', sample), 'Young', 'Old')) %>% mutate(time.point=ifelse(grepl('D14', sample), 'D14', 'D28'))
  colnames(cts)[1] <- "Value"
  cts$gene <- gene
  cts
})
write.xlsx(Reduce(rbind,df.list),file = "Fig5E_data.xlsx")

# IPA based on differentially expressed genes between physiological Old and Young gCap 
old.cp <- c("Th1 and Th2 Activation Pathway","Antigen Presentation Pathway","IL-10 Signaling","Integrin cell surface interactions","IL-4 Signaling","Collagen biosynthesis and modifying enzymes","Pulmonary Fibrosis Idiopathic Signaling Pathway","Platelet Adhesion to exposed collagen","Interferon gamma signaling","S100 Family Signaling Pathway","FGF Signaling","Extracellular matrix organization")
young.cp <- c("Apelin Endothelial Signaling Pathway","CXCR4 Signaling","Thrombin Signaling")

DF_UP <- read.xlsx("~/Pathways_IPA_OldvsYoung_TableS5.xlsx",sheet = "TableS5_UP")
DF_DOWN <- read.xlsx("~/Pathways_IPA_OldvsYoung_TableS5.xlsx",sheet = "TableS5_Down")

DF_UP <- data.frame("DF"=DF_UP$Ingenuity.Canonical.Pathways,pval=DF_UP$padj,score=DF_UP$`z-score`) %>% filter(DF%in%old.cp)
DF_UP$pval <- round(as.numeric(DF_UP$pval),digits=2)
DF_UP <- DF_UP %>% arrange(-pval)
DF_UP$group <- "UP"

DF_DOWN <- data.frame("DF"=DF_DOWN$Ingenuity.Canonical.Pathways,pval=DF_DOWN$padj,score=DF_DOWN$`z-score`) %>% filter(DF%in%young.cp)
DF_DOWN$pval <- round(as.numeric(DF_DOWN$pval),digits=2)
DF_DOWN <- DF_DOWN %>% arrange(-pval)
DF_DOWN$group <- "DOWN"

DF_UP <- rbind(DF_UP,DF_DOWN)
DF_UP$DF <- factor(DF_UP$DF,levels=rev(DF_UP$DF))


# ---- Fig5F IPA on Old vs Young PBS gCap ----
ggplot(data=DF_UP, aes(x=DF, y=pval, fill=group)) +
  geom_bar(stat="identity", color="black",width=0.8)+coord_flip()+theme_classic()+
  scale_fill_manual(values = c("DOWN"='#dba75e',"UP"='#917dd1'))+ylab("-log10(padj)")+xlab("")


# Commonly differentially expressed genes in Bleo vs PBS and Old vs Young PBS
OldvsYoung.gCap.PBS_pos <- read.xlsx("~/Supplementary_table_S6.xlsx",sheet="gCap")%>% filter(padj <0.05 & log2FoldChange>0)
OldvsYoung.gCap.PBS_neg <- read.xlsx("~/Supplementary_table_S6.xlsx",sheet="gCap")%>% filter(padj <0.05 & log2FoldChange<0)
BleovsPBS.gCap_pos <- read.xlsx("~/Supplementary_table_S3.xlsx") %>% filter(celltype=="gCap" & baseMean > 3 & log2FoldChange>0)
BleovsPBS.gCap_neg <- read.xlsx("~/Supplementary_table_S3.xlsx") %>% filter(celltype=="gCap" & baseMean > 3 & log2FoldChange<0)
ggvenn(list("OldvsYoung.gCap.PBS_pos"=OldvsYoung.gCap.PBS_pos$gene,"BleovsPBS.gCap_pos"=BleovsPBS.gCap_pos$gene),fill_color = c(gg_color_hue(2)[1], gg_color_hue(2)[2]),stroke_size=0.4,text_size = 3)+ggtitle("")
ggvenn(list("OldvsYoung.gCap.PBS_neg"=OldvsYoung.gCap.PBS_neg$gene,"BleovsPBS.gCap_neg"=BleovsPBS.gCap_neg$gene),fill_color = c(gg_color_hue(2)[1], gg_color_hue(2)[2]),stroke_size=0.4,text_size = 3)+ggtitle("")

intersect(OldvsYoung.gCap.PBS_pos$gene,BleovsPBS.gCap_pos$gene)
intersect(OldvsYoung.gCap.PBS_neg$gene,BleovsPBS.gCap_neg$gene)
genes_highlight <- c(intersect(OldvsYoung.gCap.PBS_pos$gene,BleovsPBS.gCap_pos$gene),intersect(OldvsYoung.gCap.PBS_neg$gene,BleovsPBS.gCap_neg$gene))

test <- unique(c(OldvsYoung.gCap.PBS_pos$gene,OldvsYoung.gCap.PBS_neg$gene))

OldvsYoung.gCap.PBS <- read.xlsx("~/Supplementary_table_S6.xlsx",sheet="gCap") %>% mutate(OldvsYoung=log2FoldChange)  %>% dplyr::select(gene,OldvsYoung) 
BleovsPBS.gCap <- read.xlsx("~/Supplementary_table_S3.xlsx",sheet="gCap") %>% mutate(BleovsPBS=log2FoldChange)  %>% dplyr::select(gene,BleovsPBS)
a <- inner_join(OldvsYoung.gCap.PBS,BleovsPBS.gCap) %>% filter(gene%in%test)
a$threshold = as.factor(ifelse(a$gene %ni% genes_highlight, 'not significant', 
                               ifelse(a$gene %in% intersect(OldvsYoung.gCap.PBS_neg$gene,BleovsPBS.gCap_neg$gene), 'down-regulated',
                                      ifelse(a$gene %in% intersect(OldvsYoung.gCap.PBS_pos$gene,BleovsPBS.gCap_pos$gene), 'up-regulated', 'not significant'))))
a$gene <- ifelse(a$gene %in% genes_highlight, a$gene, "")

# ---- Fig5G Commonly differentially expressed genes in Bleo vs PBS and Old vs Young PBS ----
ggplot(a, aes(x=OldvsYoung, y=BleovsPBS)) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1) +
  # scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#5088cc', 'grey', '#cc5050')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ggtitle("") + xlim(c(-3,2.6))+
  xlab("Old PBS vs Young PBS") + ylab("BLM vs PBS") +
  theme_minimal() + 
  theme(legend.position = "none") +
  geom_text_repel(aes(label = gene, color = threshold), size = 3)+FontSize(x.text =8 ,x.title = 10,y.text = 8,y.title = 10 )
ggsave("Fig5_G.pdf", width = 12, height = 12, units = c("cm"), dpi = 200)


