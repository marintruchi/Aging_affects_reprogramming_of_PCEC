# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to explore the scRNA-seq dataset and produce the plots of figure 2
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")

#--------------------------------------------------------------------------------------------#
#------------------------------- Explore the scRNA-seq dataset ------------------------------#
#--------------------------------------------------------------------------------------------#

#Loading dataset
#For convenience, load directly the seurat object obtained after the processing and integration steps 
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(aggr) <- "celltype"
sample_order <- c( "Young_D14.Bleo","Young_D14.PBS","Old_D14.Bleo","Old_D14.PBS","Young_D28.Bleo","Young_D28.PBS","Old_D28.Bleo","Old_D28.PBS","Young_D60.Bleo","Young_D60.PBS","Old_D60.Bleo","Old_D60.PBS" )
group_order <- c("PBS","D14","D28","D60")
celltype_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                    "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                    "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                    "Prolif. EC","gCap","aCap","Venous EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                    "Mesothelial cells","AT1","AT2","Multiciliated cells")

#Colors for plots
colors.group <- c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c")
colors.run=c("Young D14"="#FDB863","Young D28"="#E08214","Young D60"="#B35806","Old D14"="#B2ABD2","Old D28"="#8073AC","Old D60"="#542788")
colors.Bleo<- c("Prolif. AM"="#a6d8f5","AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0","IM"="#225EA8","LCM"="#0e3378",
                "Prolif. DC"="#ffd1ec","cDC1"="#e670b4","cDC2"="#9c119c","Mature DC"="#960c5a","pDC"="#de439a","cMonocytes"="#66119c","ncMonocytes"="#282b5c","Neutrophils"="#4300de","Neutrophil-like Monocytes"="#38015c","Mast cells"="#9315ed","Platelets"="#de1212",
                "Prolif. T cells"="#F7FCB9","CD4 T cells"="#D9F0A3","Reg T cells"="#ADDD8E","T helper cells"="#78C679","Gamma delta T cells"="#41AB5D","ILC2"="#739154","CD8 T cells"="#006837" ,"NKT1"="#004529","NKT2"="#01736b","NK cells"="#46c798","B cells"="#ebd513","Plasmocytes"="#ccbf47",
                "Prolif. EC"="#fee0d2","gCap"="#911414","aCap"="#5e1e33","Venous EC"="#f2ae5a","Arterial EC"="#e36805","Lymphatic EC"="#7d3f15","Fibroblasts"="#d3d5db","Pericytes"="#806c6c",
                "Mesothelial cells"="#120c04","AT1"="#8184b8","AT2"="#3e447a","Multiciliated cells"="#1486f7")
colors.sample<- c("Young_D14.PBS"="#a6d2ed","Old_D14.PBS"="#a6d2ed","Young_D28.PBS"="#a6d2ed","Old_D28.PBS"="#a6d2ed", "Young_D60.PBS"="#a6d2ed","Old_D60.PBS"="#a6d2ed",
                  "Young_D14.Bleo"="#e62727","Old_D14.Bleo"="#8f1414","Young_D28.Bleo"="#e05c5c","Old_D28.Bleo"="#913c3c","Young_D60.Bleo"="#e08d8d","Old_D60.Bleo"="#946666")

aggr_norm <- aggr %>% NormalizeData()
markers <- FindAllMarkers(aggr_norm,min.pct = 0.2, logfc.threshold = 0.25,only.pos = T)

# Create Supplementary_table_2 : markers of annotated clusters 
markers %>% write.xlsx(file = "Supplementary_table_2.xlsx")

#UMAPs
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.1,cols = colors.Bleo,repel = T,group.by = "celltype")+NoAxes()+NoLegend()+ggtitle("")
ggsave("Fig2_B_subpop.pdf", width = 15, height = 15, units = c("cm"), dpi = 100)

DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.1,group.by = "group",cols = colors.group)+NoAxes()+NoLegend()+ggtitle("")
ggsave("Fig2_B_group.pdf", width = 15, height = 15, units = c("cm"), dpi = 100)


aggr$epytllec <- factor(aggr$celltype,levels = rev(levels(aggr$celltype)))
aggr <- NormalizeData(aggr)
p1 <- DotPlot(aggr,assay = "RNA",scale=T,group.by = "epytllec",features = c("Mki67","Mrc1","Ctsk","Trem2","C1qc","Prg4","Mgl2","Xcr1","Cd209a","Ccl22","Siglech","F13a1","Cd300e","Adgre4","Mmp9","Fcer1a","Ppbp",
                                                                            "Cd4","Ikzf2","Cxcr6","Cd163l1","Il1rl1","Cd8b1","Ly6c2","Ccl5","Gzma","Cd79a","Jchain",
                                                                            "Sema3c","Ednrb","Slc6a2","Adgrg6","Ccl21a","Col1a1","Cox4i2","Upk3b","Aqp5","Sfta2","Ccdc153"),cols = "RdBu",dot.scale = 4) + RotatedAxis() +FontSize(x.text = 12,y.text = 12) +theme(axis.title.x = element_blank(),axis.title.y = element_blank())
p1+NoLegend()
as_ggplot(ggpubr::get_legend(p1,position = "top"))



#------------------------------------------------------------------------------------------
#Differential Expression Analyses using DESeq2 on pseudo-bulks for main subpopulations
library(DESeq2)
cellstates_BP <- c("cDC1","cDC2","cMonocytes","ncMonocytes","Neutrophils","Mast cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells",
                   "CD8 T cells","NKT1","NKT2","NK cells","B cells","gCap","aCap","Venous EC","Fibroblasts")

aggr <-subset(x = aggr, subset = cellstate %in% cellstates_BP & group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

for (CT in cellstates_BP) {
  cts <- cts.split[[CT]]
  colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
  colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
  dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >=10,]
  print(paste("Testing",CT))
  dds <- DESeq(dds,test="Wald")
  DEG.pseudobulk <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
  write.xlsx(DEG.pseudobulk,file = paste0("~/DEG.pseudobulks/",CT,"_pseudobulk.xlsx"))
}

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
cellstates_BP <- c("AM1","AM2","AM3")
aggr <-subset(x = aggr, subset = cellstate %in% cellstates_BP & group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
DEG.pseudobulk <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.pseudobulk,file = "~/DEG.pseudobulks/AM_pseudobulk.xlsx")

DE.list <- list(read.xlsx("~/DEG.pseudobulks/AM_pseudobulk.xlsx")%>% mutate(celltype="AM"),
                read.xlsx("~/DEG.pseudobulks/Fibroblasts_pseudobulk.xlsx")%>% mutate(celltype="Fibroblasts"),
                read.xlsx("~/DEG.pseudobulks/gCap_pseudobulk.xlsx")%>% mutate(celltype="gCap"),
                read.xlsx("~/DEG.pseudobulks/aCap_pseudobulk.xlsx")%>% mutate(celltype="aCap"),
                read.xlsx("~/DEG.pseudobulks/Venous EC_pseudobulk.xlsx")%>% mutate(celltype="Venous EC"),
                read.xlsx("~/DEG.pseudobulks/cDC2_pseudobulk.xlsx")%>% mutate(celltype="cDC2"),
                read.xlsx("~/DEG.pseudobulks/Neutrophils_pseudobulk.xlsx")%>% mutate(celltype="Neutrophils"),
                read.xlsx("~/DEG.pseudobulks/NKT1_pseudobulk.xlsx")%>% mutate(celltype="NKT1"),
                read.xlsx("~/DEG.pseudobulks/NKT2_pseudobulk.xlsx")%>% mutate(celltype="NKT2"),
                read.xlsx("~/DEG.pseudobulks/CD4 T cells_pseudobulk.xlsx")%>% mutate(celltype="CD4 T cells"),
                read.xlsx("~/DEG.pseudobulks/CD8 T cells_pseudobulk.xlsx")%>% mutate(celltype="CD8 T cells"),
                read.xlsx("~/DEG.pseudobulks/Reg T cells_pseudobulk.xlsx")%>% mutate(celltype="Reg T cells"))

colors.Bleo<- c("AM"="#41B6C4","Fibroblasts"="black","gCap"="#911414","aCap"="#5e1e33","Venous EC"="#f2ae5a","cDC2"="#9c119c","Neutrophils"="#4300de","NKT1"="#004529","NKT2"="#01736b",
                "CD4 T cells"="#D9F0A3","CD8 T cells"="#006837","Reg T cells"="#ADDD8E" )

#bind data together
DE.list<-bind_rows(DE.list)
DE.list<-DE.list %>% group_by(celltype) %>% mutate( sig = ifelse(padj < 0.05 & abs(log2FoldChange) > .25 , "Sig", "NS"))
cols <- data.frame(hexcode=colors.Bleo,celltype=c(names(colors.Bleo))) 
DE.list<-full_join(DE.list, cols, by = "celltype")

#set color to use
DE.list<-mutate(DE.list, col_use = ifelse(sig== "Sig", hexcode, "dark gray"))
DE.list<-mutate(DE.list, updown = ifelse(padj < 0.05 & log2FoldChange >=0.25, "Upregulated", ifelse(padj < 0.05 & log2FoldChange <= -0.25, "Downregulated", "NS")))
DE.list <- DE.list[complete.cases(DE.list), ]

de_res_table<-table(DE.list$updown, DE.list$celltype)
de_res_table<-as.data.frame.matrix(de_res_table)

#add color for sig points
col_for_plot<-as.character(DE.list$col_use)
table(col_for_plot)

DE.list<-mutate(DE.list, alp_for_plot = ifelse(updown == "NS", 0.25, 1))
alp_for_plot<-DE.list$alp_for_plot
x<-unique(DE.list$celltype)
DE.list$celltype <- factor(DE.list$celltype, levels = x)
DE.list<-group_by(DE.list, celltype) 

#get top and bottom 5 DE genes for each type
sigs<-dplyr::filter(DE.list, sig == "Sig")
top_5<-sigs %>% group_by(celltype) %>% slice_max(order_by = log2FoldChange, n = 5)
bot_5<-sigs %>% group_by(celltype) %>% slice_min(order_by = log2FoldChange, n = 5)
five<-dplyr::full_join(top_5, bot_5)
#five <- sigs

DE.list_sig<-dplyr::filter(DE.list, updown != "NS")

# Create Supplementary_table_3 : differentially expressed genes between BLM and PBS 
write.xlsx(DE.list_sig,file = "Supplementary_table_3.xlsx")

DE.list_not<-dplyr::filter(DE.list, updown == "NS")

col_for_plot<-DE.list_sig$hexcode

#Clean out low expressed genes
DE.list <- DE.list %>% filter(baseMean>5)

p1 <- ggplot(DE.list, aes( x = celltype, y = log2FoldChange,)) +
  geom_jitter_rast(data=DE.list_not,color="dark grey", width = 0.3, height = 0.0, alpha = .25, shape = 1, raster.dpi = 100) +
  geom_jitter(data=DE.list_sig,color=col_for_plot, width = 0.3, height = 0.0, shape = 18) + theme_bw() +
  theme(axis.text.x = element_text(size = 9,angle = 45, vjust = 0.5, hjust=0.5), legend.position = "none", axis.line = element_line(colour = "black"), panel.border = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #geom_text(data = five, label = five$gene, fontface = "italic", size = 2.3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p1
ggsave("Fig2_C_DEGlog2FC.pdf", width = 15, height = 8, units = c("cm"), dpi = 300)

a <- table(DE.list_sig$sig,DE.list_sig$celltype) %>% as.data.frame.matrix() %>% t()
a <- data.frame(celltype=factor(rownames(a),levels = names(colors.Bleo)),DEG=a[,1])
a[a>1500] <- 1500
p2 <- ggplot(a, aes(x=celltype, y=DEG, fill=celltype)) + 
  geom_bar(stat="identity")+ scale_fill_manual(values=colors.Bleo)+theme_classic()+ theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank())+NoLegend()
p2
ggsave("Fig2_C_DEGgeombar.pdf", width = 15, height = 3.5, units = c("cm"), dpi = 300)


#------------------------------------------------------------------------------------------
#Relative proportion boxplot of Alveolar Macrophages subpopullations
#Code from https://github.com/theislab/2019_Strunz/blob/master/code/Fig1abcdefghi.R
# Define functions for plotting relative frequencies ####

aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Macro <- subset(aggr,idents=c("Prolif. AM","AM1", "AM2", "AM3","IM","LCM"))
celltype_order <- c("Prolif. AM","AM1", "AM2", "AM3","IM","LCM")
colors.Bleo <- c("Prolif. AM"="#a6d8f5","AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0","IM"="#225EA8","LCM"="#0e3378")
Macro <- RunUMAP(Macro,dims = 1:80)
plot <- DimPlot(object = Macro,  reduction = 'umap', pt.size = 0.1,cols = colors.Bleo,label = T,repel = T,group.by = "celltype")+NoAxes()+ggtitle("")
plot+NoLegend()

Idents(Macro) <- "age"
Macro$celltype <- factor(Macro$cellstate,levels = c("Prolif. AM","AM1", "AM2", "AM3","IM","LCM"))
Macro$identifier <- Macro$mice
Macro.old <- subset(Macro,idents="Old")
Macro <- subset(Macro,idents="Young")

ggarrange(plot_relFreqs(get_relFreqs(Macro, grouping = "group"), "AM1", colors.group, group_order, title = "AM1", save = F)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()),
          plot_relFreqs(get_relFreqs(Macro.old, grouping = "group"), "AM1", colors.group, group_order, title = "AM1", save = F)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()),ncol=1,nrow=2)
ggsave("Fig2_D_AM1.pdf", width = 4, height = 10, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(Macro, grouping = "group"), "AM2", colors.group, group_order, title = "AM2", save = F)+ylim(0, 0.75)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()),
          plot_relFreqs(get_relFreqs(Macro.old, grouping = "group"), "AM2", colors.group, group_order, title = "AM2", save = F)+ylim(0, 0.75)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()),ncol=1,nrow=2)
ggsave("Fig2_D_AM2.pdf", width = 4, height = 10, units = c("cm"), dpi = 200)

ggarrange(plot_relFreqs(get_relFreqs(Macro, grouping = "group"), "AM3", colors.group, group_order, title = "AM3", save = F)+ylim(0, 0.50)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()),
          plot_relFreqs(get_relFreqs(Macro.old, grouping = "group"), "AM3", colors.group, group_order, title = "AM3", save = F)+ylim(0, 0.50)+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()),ncol=1,nrow=2)
ggsave("Fig2_D_AM3.pdf", width = 4, height = 10, units = c("cm"), dpi = 200)


