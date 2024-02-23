# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to compare PCEC from old and young mice and produce the plots of figure 4 and Supplementary figure S4
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")

#----------------------------------------------------------------------------------------------------------#
#------------------------- PCEC subpopulations monitoring between young and old mice ----------------------#
#----------------------------------------------------------------------------------------------------------#

#Loading dataset
#For convenience, load directly the seurat object obtained after the processing and integration steps
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(aggr) <- "celltype"
group_order <- c("PBS","D14","D28","D60")

#Colors for plots
colors.group <- c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c")

#------------------------------------------------------------------------------------------
#Relative proportion boxplot
#Code from https://github.com/theislab/2019_Strunz/blob/master/code/Fig1abcdefghi.R
# Define functions for plotting relative frequencies ####

Idents(aggr) <- "cellstate"
SUB <-subset(x = aggr, idents = c("aCap","gCap","Lrg1+ aCap","Lrg1+ gCap","SV EC","Prolif. EC","PV EC","Arterial EC"))
subcolors <- c("SV EC"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414","Lrg1+ aCap"="#FFD92F","aCap"="#5e1e33","Prolif. EC"="#fee0d2","PV EC"="#f2ae5a","Arterial EC"="#e36805")

SUB$celltype <- factor(SUB$cellstate,levels = c("SV EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Prolif. EC","PV EC","Arterial EC"))
SUB$identifier <- SUB$mice

Idents(SUB) <- "age"
SUB.old <- subset(SUB,idents="Old")
SUB <- subset(SUB,idents="Young")

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "SV EC", colors.group, group_order, title = "SV EC", save = F)+ylim(0, 0.1)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "SV EC", colors.group, group_order, title = "SV EC", save = F)+ylim(0, 0.1)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_SVEC.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "Lrg1+ gCap", colors.group, group_order, title = "Lrg1+ gCap", save = F)+ylim(0, 0.42)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "Lrg1+ gCap", colors.group, group_order, title = "Lrg1+ gCap", save = F)+ylim(0, 0.42)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),),ncol=1,nrow=2)
ggsave("Fig4_A_Lrg1gCap.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "gCap", colors.group, group_order, title = "gCap", save = F)+ylim(0, 0.9)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "gCap", colors.group, group_order, title = "gCap", save = F)+ylim(0, 0.9)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_gCap.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "Lrg1+ aCap", colors.group, group_order, title = "Lrg1+ aCap", save = F)+ylim(0, 0.1)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "Lrg1+ aCap", colors.group, group_order, title = "Lrg1+ aCap", save = F)+ylim(0, 0.1)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_Lrg1aCap.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "aCap", colors.group, group_order, title = "aCap", save = F)+ylim(0, 0.16)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "aCap", colors.group, group_order, title = "aCap", save = F)+ylim(0, 0.16)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_aCap.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "Prolif. EC", colors.group, group_order, title = "Prolif. EC", save = F)+ylim(0, 0.11)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "Prolif. EC", colors.group, group_order, title = "Prolif. EC", save = F)+ylim(0, 0.11)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_ProlifEC.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "PV EC", colors.group, group_order, title = "PV EC", save = F)+ylim(0, 0.16)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "PV EC", colors.group, group_order, title = "PV EC", save = F)+ylim(0, 0.16)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_PVEC.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group"), "Arterial EC", colors.group, group_order, title = "Arterial EC", save = F)+ylim(0, 0.11)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),
          plot_relFreqs(get_relFreqs(SUB.old, grouping = "group"), "Arterial EC", colors.group, group_order, title = "Arterial EC", save = F)+ylim(0, 0.11)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank()),ncol=1,nrow=2)
ggsave("Fig4_A_ArterialEC.pdf", width = 3, height = 8, units = c("cm"), dpi = 200)

#------------------------------------------------------------------------------------------
#Differentially expressed genes between Bleo and PBS
library(DESeq2)
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap","PV EC") & group %ni% c("D60"))
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
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
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 10)

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
SUB <- subset(x = aggr, subset = cellstate %in% c("Lrg1+ gCap","gCap"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.gCap.mice$gene),name="gCap_bleo")
# VlnPlot(SUB %>% subset(subset=cellstate%in% c("gCap","Lrg1+gCap")),features = "gCap_bleo1",group.by = "group",split.by = "age",sort = T)

metadata <- SUB@meta.data %>% filter(cellstate %in% c("gCap","Lrg1+ gCap")) %>% dplyr::select(gCap_bleo1,mice) %>% group_by(mice) %>% summarise_at(vars(gCap_bleo1), list(gCap_bleo1 = mean))
metadata <- metadata %>% mutate(sample = gsub("\\..*","",mice ),condition=ifelse(grepl('PBS', mice), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition))
metadata$age <- factor(metadata$age,levels = c("Young","Old"))
ggboxplot(metadata, x="group",y= "gCap_bleo1",fill = "age",ylab = "",add = "dotplot",add.params = list(dotsize=0.8), palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()
ggsave("Fig4_B_gCap.pdf", width = 7, height = 7, units = c("cm"), dpi = 200)

#Supplementary heatmap of bleomycin-induced markers in gCap 
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% top_n(50,-log10(padj))
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.gCap.mice$gene,]  %>% as.data.frame()  
names <- gsub("gCap_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

pdf("SupplFig4_A_heatmap_DEG_gCap.pdf", width = 8, height = 4)
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)
dev.off()


#aCap
cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
DEG.psdblk.aCap.mice <- DEG.psdblk.aCap.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 5)

SUB <- subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","aCap"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.aCap.mice$gene),name="gCap_bleo")

metadata <- SUB@meta.data %>% filter(cellstate %in% c("aCap","Lrg1+ aCap")) %>% dplyr::select(gCap_bleo1,mice) %>% group_by(mice) %>% summarise_at(vars(gCap_bleo1), list(gCap_bleo1 = mean))
metadata <- metadata %>% mutate(sample = gsub("\\..*","",mice ),condition=ifelse(grepl('PBS', mice), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition))
metadata$age <- factor(metadata$age,levels = c("Young","Old"))
ggboxplot(metadata, x="group",y= "gCap_bleo1",fill = "age",ylab = "",add = "dotplot",add.params = list(dotsize=0.8), palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()
ggsave("Fig4_B_aCap.pdf", width = 7, height = 7, units = c("cm"), dpi = 200)

#Supplementary heatmap of bleomycin-induced markers in aCap 
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.aCap.mice$gene,]  %>% as.data.frame()  
names <- gsub("aCap_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))
pdf("SupplFig4_B_heatmap_DEG_aCap.pdf", width = 8, height = 4)
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)
dev.off()

#(Pulmonary) Venous EC
cts <- cts.split$`Venous EC`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.PVEC.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
DEG.psdblk.PVEC.mice <- DEG.psdblk.PVEC.mice %>% filter(padj<0.05 & log2FoldChange > 0.5 & baseMean > 5)

SUB <- subset(x = aggr, subset = cellstate %in% c("PV EC"))
SUB <- SUB %>%  NormalizeData() %>% ScaleData()
SUB <- AddModuleScore(SUB,features = list(DEG.psdblk.PVEC.mice$gene),name="gCap_bleo")

metadata <- SUB@meta.data %>% filter(cellstate %in% c("PV EC")) %>% dplyr::select(gCap_bleo1,mice) %>% group_by(mice) %>% summarise_at(vars(gCap_bleo1), list(gCap_bleo1 = mean))
metadata <- metadata %>% mutate(sample = gsub("\\..*","",mice ),condition=ifelse(grepl('PBS', mice), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition))
metadata$age <- factor(metadata$age,levels = c("Young","Old"))
ggboxplot(metadata, x="group",y= "gCap_bleo1",fill = "age",ylab = "",add = "dotplot",add.params = list(dotsize=0.8), palette = c("Young"="#dba75e","Old"="#917dd1"))+FontSize(x.text = 10,y.text = 10,x.title = 0)+NoLegend()
ggsave("SupplFig4_C.pdf", width = 7, height = 7, units = c("cm"), dpi = 200)

#Supplementary heatmap of bleomycin-induced markers in PVEC
data <- data.frame(counts(dds, normalized=TRUE))[DEG.psdblk.PVEC.mice$gene,]  %>% as.data.frame()  
names <- gsub("Venous.EC_","",colnames(data))
annotation <- data.frame(row.names = colnames(data)) %>% mutate(sample = gsub("\\..*","",names ),condition=ifelse(grepl('PBS', names), 'PBS', 'Bleo'),age=gsub("_.*","",sample ),group=ifelse(condition=="Bleo",gsub(".*_","",sample),condition)) %>% dplyr::select(group,age)
data <- flexible_normalization(data)
data[data>3] <- 3

annot.cols <- list("group"=c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c"),"age"=c("Young"="#dba75e","Old"="#917dd1"))
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))
# pdf("SupplFig4_B_heatmap_DEG_aCap.pdf", width = 8, height = 4)
pheatmap(t(data),annotation_row = annotation,annotation_colors = annot.cols,cluster_rows = T,color = myColor,fontsize_col = 8,treeheight_col = 0,breaks = Breaks,show_rownames = F,angle_col = 90)
# dev.off()


#------------------------------------------------------------------------------------------
#Expression of Lrg1 in spatial transcriptomics data
#https://satijalab.org/seurat/articles/spatial_vignette.html

set3 <- brewer.pal(n = 12, name = "Set3")
colors.tissue <- c("Fibrosis"="#de5f5f","Interstitium"="#f79c4d","Macrophages"="#41B6C4","Blood vessel"=set3[12],"Bronchia"=set3[7],"Bronchioles"="#75c26e","Alveoli"="#3e447a","Ig+ Alveoli"="#6e578a","IFN+ Alveoli"="#ed95ea","low counts Alveoli" = "#c7c3e3","Blood"=set3[10],"Mesothelial"="#8c6253","Connective tissue"="black","SMC layer"="#806c6c")
tissue.order <- c("Fibrosis","Interstitium","Macrophages","Blood vessel","Bronchia","Bronchioles","Alveoli","Ig+ Alveoli","IFN+ Alveoli","low counts Alveoli","Blood","SMC layer","Mesothelial","Connective tissue" )
slice.order <- c("Young.D14_2.down.left","Old.D14_1.up.right","Old.D14_2.up.left","Young.D28_1.up.left","Old.D28_1.up.right","Old.D28_2.up.left")

# For convenience, load the seurat object of Visium data
merged <- readRDS("~/BLEO_visium_integrated.rds")
merged$group.cs <- paste(merged$sample,merged$int_tissue,sep = "-")
Idents(merged) <- "group.cs"
sub <- merged
DefaultAssay(sub) <- "Spatial"
sub <- NormalizeData(sub)

sub$id <- sub$group.cs
select.cells <- colnames(subset(sub,subset = Lrg1 >= 1))
sub1 <- subset(sub, cells = select.cells) 
Lrg1 <- data.frame(table(sub1$id))
table(sub1$id)
Lrg1 <- data.frame(table(sub$id))
Lrg1 <-full_join(data.frame(table(sub$id)),data.frame(table(sub1$id)),by="Var1") %>% mutate(percentage=Freq.y/Freq.x*100)
sub1 <- sub 
Idents(sub1) <- "id"
sub1 <- AverageExpression(sub1, return.seurat = T,assays = "Spatial")
sub1 <- data.frame(Var1=names(sub1@assays[["Spatial"]]@scale.data["Lrg1",]),Lrg1_expr=sub1@assays[["Spatial"]]@scale.data["Lrg1",]) %>% inner_join(Lrg1) %>% mutate(id=Var1) %>% dplyr::select(id,Lrg1_expr,percentage)
sub1[is.na(sub1)] <- 0

#Split ID in groupe and Condition
sub1$tissue <- factor(gsub(".*-","",sub1$id),levels = rev(tissue.order))
sub1$group <- factor(gsub("-.*","",sub1$id),levels = slice.order)

p1 <- ggplot(sub1,aes(x=group,y=tissue,size=percentage)) + 
  geom_point(aes(col=Lrg1_expr),alpha=1) +
  scale_colour_gradientn(colours = c("#B2182B" ,"#D6604D" ,"#F4A582", "#FDDBC7", "#F7F7F7" ,"#D1E5F0", "#92C5DE" ,"#4393C3"),values = scales::rescale(c(max(sub1$Lrg1_expr),1.5,1,0.5, 0, -1.25, -2.5,min(sub1$Lrg1_expr))))+
  labs(y="time point")+theme_classic()+rotate_x_text(45)+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank())
p1+NoLegend()
ggsave("Fig4_C.pdf", width = 7, height = 8, units = c("cm"), dpi = 200)


#------------------------------------------------------------------------------------------
#Heatmap of PCEC pathological and physiological genes in spatial transcriptomics data

merged <- readRDS("~/BLEO_visium_integrated.rds")
set3 <- brewer.pal(n = 12, name = "Set3")
tissue.order <- c("Fibrosis","Interstitium","Macrophages","Blood vessel","Bronchia","Bronchioles","Alveoli","Ig+ Alveoli","IFN+ Alveoli","low counts Alveoli","Blood","SMC layer","Mesothelial","Connective tissue" )
slice.order <- c("Young.D14_2.down.left","Old.D14_1.up.right","Old.D14_2.up.left","Young.D28_1.up.left","Old.D28_1.up.right","Old.D28_2.up.left")
colors.tissue <- c("Fibrosis"="#de5f5f","Interstitium"="#f79c4d","Macrophages"="#41B6C4","Blood vessel"=set3[12],"Bronchia"=set3[7],"Bronchioles"="#75c26e","Alveoli"="#3e447a","Ig+ Alveoli"="#6e578a","IFN+ Alveoli"="#ed95ea","low counts Alveoli" = "#c7c3e3","Blood"=set3[10],"Mesothelial"="#8c6253","Connective tissue"="black","SMC layer"="#806c6c")
tissue.order <- c("Fibrosis","Interstitium","Macrophages","Blood vessel","Bronchia","Bronchioles","Alveoli","Ig+ Alveoli","IFN+ Alveoli","low counts Alveoli","Blood","SMC layer","Mesothelial","Connective tissue" )
slice.order <- c("Young.D14_2.down.left","Old.D14_1.up.right","Old.D14_2.up.left","Young.D28_1.up.left","Old.D28_1.up.right","Old.D28_2.up.left")

merged$group.cs <- paste(merged$sample,merged$int_tissue,sep = "-")
Idents(merged) <- "group.cs"
average <- AverageExpression(merged, return.seurat = T)

merged$int_tissue <- factor(merged$int_tissue,levels =tissue.order )
Idents(merged) <- "int_tissue"

data <- GetAssayData(average,assay = "integrated",slot = "data")
EC.genes <- unique(c(read.xlsx("~/Supplementary_table_3.xlsx")%>% filter(celltype=="gCap" & log2FoldChange > 1) %>% top_n(30,-log10(padj)) %>% pull(gene),read.xlsx("~/Supplementary_table_3.xlsx")%>% filter(celltype=="aCap" & log2FoldChange >0.5) %>% top_n(20,-log10(padj))%>% pull(gene)))
EC.genes<- c("Sema3c","Ptprb","Calcrl","Cldn5","Epas1","Cd93","Hpgd","Tmem100","Itga1","Ramp2","Cd36","Adgrf5","Tspan7","Cdh5","Scn7a","Bmpr2","Cavin2",EC.genes)
EC.genes <- c(EC.genes[which(EC.genes %in% rownames(data))])

data <- data[match(unique(c(EC.genes)),rownames(data)),]
annotation <- data.frame(group=colnames(data),sample=gsub("-.*","",colnames(data)),age=gsub("\\..*","",colnames(data)),tissue=gsub(".*-","",colnames(data)))
annotation$sample <- factor(annotation$sample,levels =slice.order)
annotation$tissue <- factor(annotation$tissue,levels = tissue.order )
annotation <-arrange(annotation,tissue,sample)
data <- data[,match(annotation$group,colnames(data))]
annotation <- annotation %>% column_to_rownames(.,"group")
data <- flexible_normalization(data)
data[data < -2] <- -2
data[data > 2] <- 2

pdf("SupplFig4_D_heatmap_Spatial_markers.pdf", width = 10, height = 6)
pheatmap(data,annotation_col  = annotation,show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = F,angle_col =90,fontsize_col = 8,legend=TRUE,gaps_row = c(17),
         annotation_colors = list(tissue= colors.tissue,age=c("Young"="#dba75e","Old"="#917dd1"),sample=c("Young.D14_2.down.left"="#d1eced","Old.D14_1.up.right"="#9199db","Old.D14_2.up.left"="#727ddb","Young.D28_1.up.left"="#77a7a8","Old.D28_1.up.right"="#59328c","Old.D28_2.up.left"="#3f1e6b")))
dev.off()


#------------------------------------------------------------------------------------------
#Differential expression analysis between old and young aCap or gCap cells pseudobulks using DESeq2 in fibrotic condition 

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap","PV EC") & group %ni% c("PBS","D60"))
aggr$samples <- aggr$mice

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
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Old_vs_Young")

cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Old_vs_Young")

cts <- cts.split$`Venous EC`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.PVEC.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Old_vs_Young")

# Create Supplementary_table_4 : Differential expression analysis (DESeq2) between old and young aCap or gCap cells in fibrotic condition
write.xlsx(list("OldvsYoung.aCap"=DEG.psdblk.aCap.mice,"OldvsYoung.gCap"=DEG.psdblk.gCap.mice,"OldvsYoung.PVEC"=DEG.psdblk.PVEC.mice),file = "Supplementary_table_4.xlsx")

gene_list <- read.xlsx("Supplementary_table_4.xlsx",sheet="OldvsYoung.gCap")
# remove unexpressed or non informative genes
gene_list <-gene_list %>% filter(baseMean > 5 & gene %ni% unlist(findnoisygenes.mm(aggr)))

upgene_Bleo <- gene_list %>% arrange(dplyr::desc(log2FoldChange)) %>% filter(log2FoldChange > 0.25 & padj<0.05)  %>% na.omit() %>% pull(gene)
dwgene_Bleo <- arrange(gene_list, log2FoldChange) %>% filter(log2FoldChange < -0.25 & padj<0.05)  %>%  na.omit() %>% pull(gene)
genestoshow <- c(upgene_Bleo, dwgene_Bleo)
print(genestoshow)

gene_list$gene <- ifelse(gene_list$gene %in% genestoshow, gene_list$gene, "")

gene_list$threshold = as.factor(ifelse(gene_list$padj > 0.05, 'not significant', 
                                       ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange <= -0.25, 'down-regulated',
                                              ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange >= 0.25, 'up-regulated', 'not significant'))))

ggplot(gene_list, aes(x=log2FoldChange, y=-log10(padj))) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1.25) +
  # scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#dba75e', 'grey', '#917dd1')) +
  ggtitle("") +
  xlab("log2 fold change") + ylab("-log10 padj") +
  #geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") +
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(data=filter(gene_list, padj<0.05), aes(label = gene, color = threshold), size = 2.8)
ggsave("Fig4_F.pdf", width = 10, height = 8, units = c("cm"), dpi = 200)

ggvenn(list("dwgene_Bleo"=dwgene_Bleo,"dwgene_PBS"=dwgene_PBS),fill_color = c(gg_color_hue(2)[1], gg_color_hue(2)[2]),stroke_size=0.4,text_size = 3)+ggtitle("")
ggvenn(list("upgene_Bleo"=upgene_Bleo,"upgene_PBS"=upgene_PBS),fill_color = c(gg_color_hue(2)[1], gg_color_hue(2)[2]),stroke_size=0.4,text_size = 3)

intersect(dwgene_Bleo,dwgene_PBS)
intersect(upgene_Bleo,upgene_PBS)


#------------------------------------------------------------------------------------------
#Differential expression analysis between old and young aCap or gCap cells pseudobulks using DESeq2 in physiological condition 

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap","PV EC") & group %in% c("PBS"))
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

cts <- cts.split$`gCap`
df_PBS <- cts
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
colData_PBS <- colData

dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Old_vs_Young")
dds_PBS <- dds
DEG.psdblk.gCap.mice_PBS <- DEG.psdblk.gCap.mice

cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Old_vs_Young")

cts <- cts.split$`Venous EC`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Young', samples), 'Young', 'Old')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Young', 'Old'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds)
resultsNames(dds)
DEG.psdblk.PVEC.mice <- results(dds, name = "condition_Old_vs_Young")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Old_vs_Young")

# Create Supplementary_table_5 : Differential expression analysis (DESeq2) between old and young aCap or gCap cells in fibrotic condition
write.xlsx(list("OldvsYoung.aCap"=DEG.psdblk.aCap.mice,"OldvsYoung.gCap"=DEG.psdblk.gCap.mice,"OldvsYoung.PVEC"=DEG.psdblk.PVEC.mice),file = "Supplementary_table_5.xlsx")


gene_list <- read.xlsx("Supplementary_table_5.xlsx",sheet="OldvsYoung.gCap")
# remove unexpressed or non informative genes
gene_list <-gene_list %>% filter(baseMean > 5 & gene %ni% unlist(findnoisygenes.mm(aggr)))

upgene_PBS <- gene_list %>% arrange(dplyr::desc(log2FoldChange)) %>% filter(log2FoldChange > 0.25 & padj<0.05)  %>% na.omit()%>% pull(gene)
dwgene_PBS <- arrange(gene_list, log2FoldChange) %>% filter(log2FoldChange < -0.25 & padj<0.05)  %>%  na.omit() %>% pull(gene)
genestoshow <- c(upgene_PBS , dwgene_PBS )

Bleo.genes.gCap.table<- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Figures/DEG_gCap_pseudobulk.xlsx")
Bleo.genes.gCap<-genestoshow[which(genestoshow %in% Bleo.genes.gCap.table$gene)]

gene_list$gene <- ifelse(gene_list$gene %in% genestoshow, gene_list$gene, "")
gene_list$threshold = as.factor(ifelse(gene_list$padj > 0.05, 'not significant', 
                                       ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange <= -0.25, 'down-regulated',
                                              ifelse(gene_list$padj < 0.05 & gene_list$log2FoldChange >= 0.25, 'up-regulated', 'not significant'))))

ggplot(gene_list, aes(x=log2FoldChange, y=-log10(padj))) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 1.25) +
  # scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#dba75e', 'grey', '#917dd1')) +
  ggtitle("") +
  xlab("log2 fold change") + ylab("-log10 padj") +
  #geom_hline(yintercept = -log10(0.05), colour="#990000", linetype="dashed") + geom_vline(xintercept = 0.5, colour="#990000", linetype="dashed") + geom_vline(xintercept = -0.5, colour="#990000", linetype="dashed") +
  theme_classic() + 
  theme(legend.position = "none") +
  geom_text_repel(data=filter(gene_list, padj<0.05), aes(label = gene, color = threshold), size = 2.8)
ggsave("Fig4_G.pdf", width = 10, height = 8, units = c("cm"), dpi = 200)






#------------------------------------------------------------------------------------------
#Commonly differentially expressed genes in Bleo vs PBS young and Old vs Young PBS
OldvsYoung.gCap.PBS_pos <- read.xlsx("~/Supplementary_table_5.xlsx",sheet="OldvsYoung.gCap")%>% filter(baseMean > 3 & padj <0.05 & log2FoldChange>0)
OldvsYoung.gCap.PBS_neg <- read.xlsx("~/Supplementary_table_5.xlsx",sheet="OldvsYoung.gCap")%>% filter(baseMean > 3 & padj <0.05 & log2FoldChange<0)
BleovsPBS.gCap_pos <- read.xlsx("~/Supplementary_table_3.xlsx") %>% filter(celltype=="gCap" & baseMean > 3 & log2FoldChange>0)
BleovsPBS.gCap_neg <- read.xlsx("~/Supplementary_table_3.xlsx") %>% filter(celltype=="gCap" & baseMean > 3 & log2FoldChange<0)
ggvenn(list("OldvsYoung.gCap.PBS_pos"=OldvsYoung.gCap.PBS_pos$gene,"BleovsPBS.gCap_pos"=BleovsPBS.gCap_pos$gene),fill_color = c(gg_color_hue(2)[1], gg_color_hue(2)[2]),stroke_size=0.4,text_size = 3)+ggtitle("")
ggvenn(list("OldvsYoung.gCap.PBS_neg"=OldvsYoung.gCap.PBS_neg$gene,"BleovsPBS.gCap_neg"=BleovsPBS.gCap_neg$gene),fill_color = c(gg_color_hue(2)[1], gg_color_hue(2)[2]),stroke_size=0.4,text_size = 3)+ggtitle("")

intersect(OldvsYoung.gCap.PBS_pos$gene,BleovsPBS.gCap_pos$gene)
intersect(OldvsYoung.gCap.PBS_neg$gene,BleovsPBS.gCap_neg$gene)
genes_highlight <- c(intersect(OldvsYoung.gCap.PBS_pos$gene,BleovsPBS.gCap_pos$gene),intersect(OldvsYoung.gCap.PBS_neg$gene,BleovsPBS.gCap_neg$gene))

test <- unique(c(OldvsYoung.gCap.PBS_pos$gene,OldvsYoung.gCap.PBS_neg$gene))

OldvsYoung.gCap.PBS <- read.xlsx("~/Supplementary_table_5.xlsx",sheet="OldvsYoung.gCap") %>% mutate(OldvsYoung=log2FoldChange)  %>% dplyr::select(gene,OldvsYoung) 
BleovsPBS.gCap <- read.xlsx("~/DEG.pseudobulks/gCap_pseudobulk.xlsx") %>% mutate(BleovsPBS=log2FoldChange)  %>% dplyr::select(gene,BleovsPBS)
a <- inner_join(OldvsYoung.gCap.PBS,BleovsPBS.gCap) %>% filter(gene%in%test)
a$threshold = as.factor(ifelse(a$gene %ni% genes_highlight, 'not significant', 
                               ifelse(a$gene %in% intersect(OldvsYoung.gCap.PBS_neg$gene,BleovsPBS.gCap_neg$gene), 'down-regulated',
                                      ifelse(a$gene %in% intersect(OldvsYoung.gCap.PBS_pos$gene,BleovsPBS.gCap_pos$gene), 'up-regulated', 'not significant'))))
a$gene <- ifelse(a$gene %in% genes_highlight, a$gene, "")
ggplot(a, aes(x=OldvsYoung, y=BleovsPBS)) +
  ggrastr::geom_point_rast(aes(color = threshold), shape = 20, alpha = 1, size = 2) +
  # scale_x_continuous(limits = c(-2, 2)) +
  scale_color_manual(values = c('#5088cc', 'grey', '#cc5050')) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  ggtitle("") + xlim(c(-3,2.6))+
  xlab("Old vs Young") + ylab("Bleo vs PBS") +
  theme_minimal() + 
  theme(legend.position = "none") +
  geom_text_repel(aes(label = gene, color = threshold), size = 3)
ggsave("Fig4_I.pdf", width = 18, height = 8, units = c("cm"), dpi = 200)


DF_list <- c("Angiogenesis","Cell movement","Migration of cells","Sprouting","Apoptosis","Cell viability","Activation of antigen presenting cells","Recruitment of leukocytes","Inflammatory response")
DF_PBS <- read.xlsx("~/Old_Young_gCap_IPA.xlsx",sheet = "DF_PBS")

DF_PBS <- data.frame("DF"=DF_PBS$Diseases.or.Functions.Annotation,pval=-log10(DF_PBS$`p-value`),score=DF_PBS$`Activation.z-score`) %>% filter(DF%in%DF_list)
DF_PBS$score <- round(as.numeric(DF_PBS$score),digits=2)
DF_PBS <- DF_PBS[order(abs(DF_PBS$score)),]
DF_PBS$DF <- factor(DF_PBS$DF,levels = DF_PBS$DF) 

ggplot(data=DF_PBS, aes(x=DF, y=score)) +
  geom_bar(stat="identity", color="black", fill=c(rep("#917dd1",4),"#dba75e",rep("#917dd1",4)),width=0.8)+coord_flip()+theme_classic()
ggsave("Fig4_H.pdf", width = 9, height = 8, units = c("cm"), dpi = 100)
