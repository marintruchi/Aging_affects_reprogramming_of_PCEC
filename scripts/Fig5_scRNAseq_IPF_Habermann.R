# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to explore PCEC subpopulations in the IPF cell atlas dataset of Habermann et al. and produce the plots of figure 5
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")
library(nichenetr)# to use the function "convert_mouse_to_human_symbols"

#----------------------------------------------------------------------------------------------------------#
#-------------------------- Relevance of PCEC's bleomycin-induced signatures in IPF -----------------------#
#----------------------------------------------------------------------------------------------------------#

#For convenience, load directly the seurat object obtained after the processing and integration steps
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
subcolors <- c("SV EC"="#A6D854","gCap"="#911414","aCap"="#5e1e33","Arterial EC"="#e37e12","PV EC"="#f2ae5a")

# Seurat object of the Habermann et al. dataset reintegrated and reannoted
Habermann <- readRDS("~/Habermann_integrated.rds")
SUB <- subset(Habermann,idents=c("SV EC","Arterial EC","aCap","gCap","PV EC"))
DefaultAssay(SUB) <- "RNA" 

SUB$samples <- paste0(SUB$Diagnosis,".",SUB$orig.ident)
SUB <- SUB %>% NormalizeData() %>% ScaleData()
SUB <- RunUMAP(SUB,dims = 1:50,assay = "integrated")

col.diagnosis <- c("Control"=gg_color_hue(2)[2],"IPF"=gg_color_hue(2)[1])
DimPlot(object = SUB, reduction = 'umap',cols = subcolors )+NoLegend()+NoAxes()
ggsave("Fig5_A.pdf", width = 8, height = 8, units = c("cm"), dpi = 200)

#------------------------------------------------------------------------------------------
#Relative proportion boxplot
#Code from https://github.com/theislab/2019_Strunz/blob/master/code/Fig1abcdefghi.R
# Define functions for plotting relative frequencies ####
SUB$celltype <- factor(SUB$celltype_MT,levels=c('SV EC',"PV EC","Arterial EC","gCap","aCap"))
SUB$identifier <- SUB$orig.ident
Idents(SUB) <- "celltype"

plot_relFreqs(get_relFreqs(SUB, grouping = "Diagnosis"), "SV EC", col.diagnosis, c("Control","IPF"), title = "", save = F)+theme_classic()+ggtitle("")+FontSize(x.text = 0,y.text = 0,x.title = 0,y.title = 0)
ggsave("Fig5_B_SVEC.pdf", width = 3.5, height = 6, units = c("cm"), dpi = 100)
plot_relFreqs(get_relFreqs(SUB, grouping = "Diagnosis"), "PV EC", col.diagnosis, c("Control","IPF"), title = "", save = F)+theme_classic()+ggtitle("")+FontSize(x.text = 0,y.text = 0,x.title = 0,y.title = 0)
ggsave("Fig5_B_PVEC.pdf", width = 3.5, height = 6, units = c("cm"), dpi = 100)
plot_relFreqs(get_relFreqs(SUB, grouping = "Diagnosis"), "gCap", col.diagnosis, c("Control","IPF"), title = "", save = F)+theme_classic()+ggtitle("")+FontSize(x.text = 0,y.text = 0,x.title = 0,y.title = 0)
ggsave("Fig5_B_gCap.pdf", width = 3.5, height = 6, units = c("cm"), dpi = 100)
plot_relFreqs(get_relFreqs(SUB, grouping = "Diagnosis"), "aCap", col.diagnosis, c("Control","IPF"), title = "", save = F)+theme_classic()+ggtitle("")+FontSize(x.text = 0,y.text = 0,x.title = 0,y.title = 0)
ggsave("Fig5_B_aCap.pdf", width = 3.5, height = 6, units = c("cm"), dpi = 100)
plot_relFreqs(get_relFreqs(SUB, grouping = "Diagnosis"), "Arterial EC", col.diagnosis, c("Control","IPF"), title = "", save = F)+theme_classic()+ggtitle("")+FontSize(x.text = 0,y.text = 0,x.title = 0,y.title = 0)
ggsave("Fig5_B_arterialEC.pdf", width = 3.5, height = 6, units = c("cm"), dpi = 100)


#------------------------------------------------------------------------------------------
#Violin plots of Lrg1+ PCEC markers in IPF
SUB <- subset(SUB,idents=c('SV EC','PV EC',"aCap","gCap"))
VlnPlot(SUB,features = "LRG1",group.by = "celltype",split.by = "Diagnosis",cols = col.diagnosis,sort = F,pt.size = 0)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=0.6)) +labs(y = "LRG1")+FontSize(x.text = 12,y.text = 8)+ labs(y = "",x="")+ggtitle("")+NoLegend()
ggsave("Fig5_C_LRG1.pdf", width = 9, height = 7, units = c("cm"), dpi = 100)

VlnPlot(SUB,features = "COL15A1",group.by = "celltype",split.by = "Diagnosis",cols = col.diagnosis,sort = F,pt.size = 0)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=0.6)) +labs(y = "COL15A1")+FontSize(x.text = 12,y.text = 8)+ labs(y = "",x="")+ggtitle("")+NoLegend()
ggsave("Fig5_C_COL15A1.pdf", width = 9, height = 7, units = c("cm"), dpi = 100)

VlnPlot(SUB,features = "SERPINE1",group.by = "celltype",split.by = "Diagnosis",cols = col.diagnosis,sort = F,pt.size = 0)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=0.6)) +labs(y = "SERPINE1")+FontSize(x.text = 12,y.text = 8)+ labs(y = "",x="")+ggtitle("")+NoLegend()
ggsave("Fig5_C_SERPINE1.pdf", width = 9, height = 7, units = c("cm"), dpi = 100)

VlnPlot(SUB,features = "EDNRB",group.by = "celltype",split.by = "Diagnosis",cols = col.diagnosis,sort = F,pt.size = 0)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=0.6)) +labs(y = "EDNRB")+FontSize(x.text = 12,y.text = 8)+ labs(y = "",x="")+ggtitle("")+NoLegend()
ggsave("Fig5_C_EDNRB.pdf", width = 9, height = 7, units = c("cm"), dpi = 100)

VlnPlot(SUB,features = "APLNR",group.by = "celltype",split.by = "Diagnosis",cols = col.diagnosis,sort = F,pt.size = 0)+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=0.6)) +labs(y = "APLNR")+FontSize(x.text = 12,y.text = 8)+ labs(y = "",x="")+ggtitle("")+NoLegend()
ggsave("Fig5_C_APLNR.pdf", width = 9, height = 7, units = c("cm"), dpi = 100)


#------------------------------------------------------------------------------------------
# Compare PCEC signatures induced by bleomycin or by IPF

#Perform differential expression between IPF and Control pseudobulks
cts <- AggregateExpression(SUB,group.by = c("celltype_MT","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))

# split data.frame
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){rownames(x) <- gsub('.*_', '', rownames(x))
t(x)})

# For gCap
cts <- cts.split$`gCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
DEG.psdblk.gCap <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")

# For aCap
cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
DEG.psdblk.aCap <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")

# For PV EC
cts <- cts.split$`PV EC`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
DEG.psdblk.PVEC <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")

# For Arterial EC
cts <- cts.split$`Arterial EC`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
DEG.psdblk.ArterialEC <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")

#For SVEC (merge SVEC and PVEC from controls to have enough cells for pseudobulks)
Habermann$id <- paste0(Habermann$celltype_MT,"-",Habermann$Diagnosis)
SUB <- subset(Habermann,subset= id%in% c("PV EC-Control","SV EC-Control","SV EC-IPF"))
SUB$samples <- paste0(SUB$Diagnosis,".",SUB$orig.ident)
DefaultAssay(SUB) <- "RNA" 

cts <- AggregateExpression(SUB,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
DEG.psdblk.SVEC <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")

# Create Supplementary_table_6 : Differential expression analysis (DESeq2) between IPF patients and healthy control for aCap, gCap, SV EC, PV EC and Arterial EC from the Habermann et al. Dataset
write.xlsx(list("aCap"=DEG.psdblk.aCap,"gCap"=DEG.psdblk.gCap,"SVEC"=DEG.psdblk.SVEC,"PVEC"=DEG.psdblk.PVEC,"ArterialEC"=DEG.psdblk.ArterialEC),file ="Supplementary_table_6.xlsx")


# To compare PCEC signatures induced by bleomycin or by IPF, load differential expression tables obtained for Lrg1+PCEC (see script for Fig3)
SVEC_IPF_1 <- DEG.psdblk.SVEC %>% filter(padj < 0.05) %>% pull(gene) 
SVEC_Bleo_1 <- read.xlsx("~/DEG.pseudobulks/SVEC_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
SVEC_Bleo_1 <- SVEC_Bleo_1[complete.cases(SVEC_Bleo_1)]

intersect_SVEC <- intersect(SVEC_Bleo,SVEC_IPF_1)

SVEC_IPF <- DEG.psdblk.SVEC %>% filter(gene %in% intersect_SVEC) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
SVEC_Bleo <- read.xlsx("~/DEG.pseudobulks/SVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
SVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(SVEC_Bleo$gene),Bleo=SVEC_Bleo$Bleo) %>% filter(gene %in% intersect_SVEC)
SVEC <- merge(SVEC_Bleo,SVEC_IPF)
SVEC$comparison <- sign(SVEC$Bleo) == sign(SVEC$IPF)
SVEC <- SVEC %>% filter(comparison == TRUE)
ggvenn(list(IPF=SVEC_IPF$gene,Bleo=SVEC_Bleo$gene),show_percentage = T,
       fill_color = c("#5e1e33", "#5e1e33"),stroke_size = 0.3, set_name_size = 3,text_size =3)

#Numbers for Venn Diagramm of figure 5E
nrow(SVEC)
length(SVEC_IPF_1)-nrow(SVEC)
length(SVEC_Bleo_1)-nrow(SVEC)

SVECgenes <- c("APLN","APLNR","CCND1","CD34","COL15A1","COL4A1","COL4A2","CXCL12","EBF1","ENPP2","GJA1","HIF1A","IGFBP7","INHBB","KDR","LAMA4","LAMB1","LAMC1","NREP","PDGFB","PDLIM1","RGCC","SMAD1","SOX17","SOX4","SPARC","SPRY1","THBS1","TNFRSF10B","VWA1",
               "ACVRL1","BMPR2","CEBPD","CRYAB","DUSP1","ELN","FOXF1","FGFR3","GATA2","GATA6","HES1","KLF2","KLF4","KLF9","PTGS1","SMAD6","SMAD7","SORBS1","SRGN","TEK","TGFB2","THBD","TIMP3","TM6SF1","VEGFC","VWF",
               "STXBP6","ABI3BP","TMEM100","PTGIS","HDAC9","ARL4D","INMT","ADGRG6","LTBP1","SAMD5","MGLL","MEOX1","FAM13C","PRDM1","FLT4","RASSF2","ANXA3","ANKRD44","EMP3")

SVEC_IPF <- DEG.psdblk.SVEC  %>% mutate(SVEC_IPF=log2FoldChange)%>% dplyr::select(gene,SVEC_IPF)
SVEC_Bleo <- read.xlsx("~/DEG.pseudobulks/SVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
SVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(SVEC_Bleo$gene),SVEC_Bleo=SVEC_Bleo$Bleo) 
SVEC <- merge(SVEC_Bleo,SVEC_IPF)%>% filter(gene%in%SVECgenes)

PVEC_IPF <- DEG.psdblk.PVEC %>% mutate(PVEC_IPF=log2FoldChange)%>% dplyr::select(gene,PVEC_IPF)
PVEC_Bleo <- read.xlsx("~/DEG.pseudobulks/PVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
PVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(PVEC_Bleo$gene),PVEC_Bleo=PVEC_Bleo$Bleo) 
PVEC <- merge(PVEC_Bleo,PVEC_IPF)%>% filter(gene%in%SVECgenes)

gCap_IPF <- DEG.psdblk.gCap  %>% mutate(gCap_IPF=log2FoldChange)%>% dplyr::select(gene,gCap_IPF)
gCap_Bleo <- read.xlsx("~/DEG.pseudobulks/Lrg1pos_gCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
gCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(gCap_Bleo$gene),gCap_Bleo=gCap_Bleo$Bleo) 
gCap <- merge(gCap_Bleo,gCap_IPF)%>% filter(gene%in%SVECgenes)

aCap_IPF <- DEG.psdblk.aCap  %>% mutate(aCap_IPF=log2FoldChange)%>% dplyr::select(gene,aCap_IPF)
aCap_Bleo <- read.xlsx("~/DEG.pseudobulks/Lrg1pos_aCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
aCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(aCap_Bleo$gene),aCap_Bleo=aCap_Bleo$Bleo) 
aCap <- merge(aCap_Bleo,aCap_IPF) %>% filter(gene%in%SVECgenes)

data <- Reduce(f = full_join,list(SVEC,PVEC,gCap,aCap)) %>% column_to_rownames(.,"gene")
cellstate_order <- c("SVEC","PVEC","gCap","aCap")
sample.order <- c("IPF","Bleo")
annotation <- data.frame(group=colnames(data),sample=gsub(".*_","",colnames(data)),celltype=gsub("_.*","",colnames(data)))
annotation$sample <- factor(annotation$sample,levels = sample.order)
annotation$celltype <- factor(annotation$celltype,levels = cellstate_order)
annotation <-arrange(annotation,celltype,sample)
data <- data[,match(annotation$group,colnames(data))]
annotation <- annotation %>% column_to_rownames(.,"group")
data[data < -3] <- -3
data[data > 3] <- 3
data[is.na(data)] <- 0
subcolors <- c("SVEC"="#A6D854","PVEC"="#f2ae5a","gCap"="#911414","aCap"="#5e1e33")
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))
pdf("Fig5_E.pdf", width = 10, height = 2)
pheatmap(t(data),annotation_row  = annotation,show_rownames = F,cluster_cols = T,cluster_rows = F,color = myColor,angle_col =90,fontsize_col = 8,treeheight_col=0,
         annotation_colors = list(celltype= subcolors,sample=c("IPF"="#5e5e5e","Bleo"="#b6b6b8")))
dev.off()

#For Table 2
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- SVEC$gene
df <- listAttributes(mart)
G_list <- getBM(filters= "hgnc_symbol", attributes= c("description","external_gene_name"),values=genes,mart= mart) %>% mutate(gene=external_gene_name)

SVEC <- merge(SVEC,G_list,by.x="gene")
SVEC$description <- gsub("\\[Source.*","",SVEC$description)
SVEC$expression_in_fibrosis <- ifelse(SVEC$SVEC_Bleo>0,"upregulated","downregulated")
SVEC <- SVEC %>%  dplyr::select(gene,description,expression_in_fibrosis)
write.xlsx(SVEC,file = "~/Table2.xlsx")


#aCap signature in bleo and in IPF
aCap_IPF_1 <- DEG.psdblk.aCap %>% filter(padj < 0.05) %>% pull(gene) 
aCap_Bleo_1 <- read.xlsx("~/DEG.pseudobulks/Lrg1pos_aCap_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
aCap_Bleo_1 <- aCap_Bleo_1[complete.cases(aCap_Bleo_1)]

intersect_aCap <- intersect(aCap_Bleo_1,aCap_IPF_1)
aCap_IPF <- DEG.psdblk.aCap %>% filter(gene %in% intersect_aCap) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
aCap_Bleo <- read.xlsx("~/DEG.pseudobulks/Lrg1pos_aCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
aCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(aCap_Bleo$gene),Bleo=aCap_Bleo$Bleo) %>% filter(gene %in% intersect_aCap)
aCap <- merge(aCap_Bleo,aCap_IPF)
aCap$comparison <- sign(aCap$Bleo) == sign(aCap$IPF)
aCap <- aCap %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5F
nrow(aCap)
length(aCap_IPF_1)-nrow(aCap)
length(aCap_Bleo_1)-nrow(aCap)

#gCap signature in bleo and in IPF
gCap_IPF_1 <- DEG.psdblk.gCap %>% filter(padj < 0.05) %>% pull(gene) 
gCap_Bleo_1 <- read.xlsx("~/DEG.pseudobulks/Lrg1pos_gCap_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
gCap_Bleo_1 <- gCap_Bleo_1[complete.cases(gCap_Bleo_1)]

intersect_gCap <- intersect(gCap_Bleo_1,gCap_IPF_1)

gCap_IPF <- DEG.psdblk.gCap %>% filter(gene %in% intersect_gCap) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
gCap_Bleo <- read.xlsx("~/DEG.pseudobulks/Lrg1pos_gCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
gCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(gCap_Bleo$gene),Bleo=gCap_Bleo$Bleo) %>% filter(gene %in% intersect_gCap)
gCap <- merge(gCap_Bleo,gCap_IPF)
gCap$comparison <- sign(gCap$Bleo) == sign(gCap$IPF)
gCap <- gCap %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5F
nrow(gCap)
length(gCap_IPF_1)-nrow(gCap)
length(gCap_Bleo_1)-nrow(gCap)


#PVEC signature in bleo and in IPF
PVEC_IPF_1 <- DEG.psdblk.PVEC %>% filter(padj < 0.05) %>% pull(gene) 
PVEC_Bleo_1 <- read.xlsx("~/DEG.pseudobulks/PVEC_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
PVEC_Bleo_1 <- PVEC_Bleo_1[complete.cases(PVEC_Bleo_1)]

intersect_PVEC <- intersect(PVEC_Bleo_1,PVEC_IPF_1)

PVEC_IPF <- DEG.psdblk.PVEC %>% filter(gene %in% intersect_PVEC) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
PVEC_Bleo <- read.xlsx("~/DEG.pseudobulks/PVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
PVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(PVEC_Bleo$gene),Bleo=PVEC_Bleo$Bleo) %>% filter(gene %in% intersect_PVEC)
PVEC <- merge(PVEC_Bleo,PVEC_IPF)
PVEC$comparison <- sign(PVEC$Bleo) == sign(PVEC$IPF)
PVEC <- PVEC %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5F
nrow(PVEC)
length(PVEC_IPF_1)-nrow(PVEC)
length(PVEC_Bleo_1)-nrow(PVEC)








