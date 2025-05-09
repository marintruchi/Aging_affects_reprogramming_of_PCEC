# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to explore the Chromium scRNA-seq dataset and produce the plots of figure 2 and supplementary figure S2 & S3
# Author : Marin Truchi

#Load required packages and functions
source("~/prerequisites.R")


#Loading dataset
#For convenience, load directly the seurat object obtained after the processing and integration steps 
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
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


young <- aggr %>% subset(subset=age=="Young")
old <- aggr %>% subset(subset=age=="Old")

# ---- SupplFigS2A VlnPlot QC ---- 
ggarrange(VlnPlot(young,group.by = "orig.ident",cols = colors.run,features = c("nCount_RNA"),pt.size = 0)+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 8,main = 12)+NoLegend(),
          VlnPlot(young,group.by = "orig.ident",cols = colors.run,features = c("nFeature_RNA"),pt.size = 0)+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 8,main = 12)+NoLegend(),
          VlnPlot(young,group.by = "orig.ident",cols = colors.run,features = c("percent_mito"),pt.size = 0)+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 8,main = 12)+NoLegend(),ncol = 3,nrow = 1)

ggarrange(VlnPlot(old,group.by = "orig.ident",cols = colors.run,features = c("nCount_RNA"),pt.size = 0)+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 8,main = 12)+NoLegend(),
          VlnPlot(old,group.by = "orig.ident",cols = colors.run,features = c("nFeature_RNA"),pt.size = 0)+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 8,main = 12)+NoLegend(),
          VlnPlot(old,group.by = "orig.ident",cols = colors.run,features = c("percent_mito"),pt.size = 0)+FontSize(x.text = 8,y.text = 8,x.title = 0,y.title = 8,main = 12)+NoLegend(),ncol = 3,nrow = 1)

# ---- SupplFigS2B Number of cells per sample ---- 
run.order <- c("Young D14","Old D14","Young D28" ,"Old D28","Young D60","Old D60")
df <- table(aggr$mice) %>% as.data.frame()
meta <- aggr@meta.data %>% dplyr::select("mice","orig.ident","group") %>% mutate(Var1=mice) %>% unique() %>% inner_join(df)
meta$orig.ident <- factor(meta$orig.ident,levels = run.order)
meta <- ddply(meta, "orig.ident",transform, Freq_ypos=cumsum(Freq)- 0.5*Freq)
ggplot(meta, aes(x=orig.ident, y=Freq, fill=group)) + labs(title = "Cell Number",)+
  geom_bar(stat="identity", color = "grey30")+ scale_fill_manual(values=colors.group)+theme_classic()+xlab("")+ylab("")+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title = element_text(hjust = 0.5))+NoLegend()+scale_y_continuous(expand = c(0,0))


#Celltypes relative proportions
metadata <- aggr@meta.data 
metadata$cellstate <- factor(metadata$celltype, levels = celltype_order)
metadata = as.data.frame.matrix(table(metadata$orig.ident, metadata$cellstate))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$group = rownames(metadata)
metadata$group =factor(metadata$group, levels = run.order)
metadata = gather(metadata, cellstate, percentage, 'Prolif. AM':'Multiciliated cells')
metadata$cellstate <- factor(metadata$cellstate, levels =celltype_order)

# ---- SupplFigS2C celltypes relative proportions ---- 
ggplot(metadata, aes(x = group, y = percentage, fill = cellstate)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = colors.Bleo) +
  ggtitle("CEC relative proportions")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.49, 'cm'),
        legend.title = element_text(size = 0))+ggtitle("")+ guides(fill=guide_legend(ncol =1))+scale_y_continuous(expand = c(0,0))+RotatedAxis()


# ---- Fig2B UMAP all populations----  
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.1,cols = colors.Bleo,repel = T,group.by = "celltype")+NoAxes()+NoLegend()+ggtitle("")

# ---- SupplFigS2D UMAP cells colored per group ----  
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.1,group.by = "group",cols = colors.group)+NoAxes()+NoLegend()+ggtitle("")

# Dot plot of best marker per population
aggr$epytllec <- factor(aggr$celltype,levels = rev(levels(aggr$celltype)))
aggr <- NormalizeData(aggr)
p1 <- DotPlot(aggr,assay = "RNA",scale=T,group.by = "epytllec",features = c("Mki67","Mrc1","Ctsk","Trem2","C1qc","Prg4","Mgl2","Xcr1","Cd209a","Ccl22","Siglech","F13a1","Cd300e","Adgre4","Mmp9","Fcer1a","Ppbp",
                                                                            "Cd4","Ikzf2","Cxcr6","Cd163l1","Il1rl1","Cd8b1","Ly6c2","Ccl5","Gzma","Cd79a","Jchain",
                                                                            "Sema3c","Ednrb","Slc6a2","Adgrg6","Ccl21a","Col1a1","Cox4i2","Upk3b","Aqp5","Sfta2","Ccdc153"),cols = "RdBu",dot.scale = 4) + RotatedAxis() +FontSize(x.text = 12,y.text = 12) +theme(axis.title.x = element_blank(),axis.title.y = element_blank())
p1+NoLegend()
as_ggplot(ggpubr::get_legend(p1,position = "top"))


#Differential Expression Analyses BLM vs PBS using DESeq2 on pseudo-bulks for main subpopulations
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

# Create Supplementary_table_S3 : differentially expressed genes between BLM and PBS 
write.xlsx(DE.list_sig,file = "Supplementary_table_S3.xlsx")

DE.list_not<-dplyr::filter(DE.list, updown == "NS")

col_for_plot<-DE.list_sig$hexcode

#Clean out low expressed genes
DE.list <- DE.list %>% filter(baseMean>5)


# ---- Fig2C Differentially expressed genes BLM vs PBS ---- 
p1 <- ggplot(DE.list, aes( x = celltype, y = log2FoldChange,)) +
  geom_jitter_rast(data=DE.list_not,color="dark grey", width = 0.3, height = 0.0, alpha = .25, shape = 1, raster.dpi = 100) +
  geom_jitter(data=DE.list_sig,color=col_for_plot, width = 0.3, height = 0.0, shape = 18) + theme_bw() +
  theme(axis.text.x = element_text(size = 9,angle = 45, vjust = 0.5, hjust=0.5), legend.position = "none", axis.line = element_line(colour = "black"), panel.border = element_blank()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+ 
  #geom_text(data = five, label = five$gene, fontface = "italic", size = 2.3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
p1

a <- table(DE.list_sig$sig,DE.list_sig$celltype) %>% as.data.frame.matrix() %>% t()
a <- data.frame(celltype=factor(rownames(a),levels = names(colors.Bleo)),DEG=a[,1])
a[a>1500] <- 1500
p2 <- ggplot(a, aes(x=celltype, y=DEG, fill=celltype)) + 
  geom_bar(stat="identity")+ scale_fill_manual(values=colors.Bleo)+theme_classic()+ theme(axis.text.x = element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank())+NoLegend()
p2


# Relative proportion boxplot of Alveolar Macrophages subpopullations
# Code from https://github.com/theislab/2019_Strunz/blob/master/code/Fig1abcdefghi.R
# Define functions for plotting relative frequencies 

#Loading dataset
#For convenience, load directly the seurat object obtained after the processing and integration steps
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
group_order <- c("PBS","D14","D28","D60")

#Colors for plots
colors.group <- c("PBS"="#8dbee0","D14"="#de5f5f","D28"="#964141","D60"="#663c3c")

aggr$identifier <- aggr$mice

aggr <- subset(aggr,subset = cellstate %in% c("AM1","AM2","AM3"))
colors.Bleo <- c("AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0")
aggr <- RunUMAP(aggr,dims = 1:80)
DimPlot(object = aggr,  reduction = 'umap', pt.size = 0.5,cols = colors.Bleo,label = T,repel = T,group.by = "celltype")+NoAxes()+NoLegend()+labs(title = "")
ggsave("UMAP_AM.pdf", width = 10, height = 10, units = c("cm"), dpi = 200)

DefaultAssay(aggr) <- "RNA"
aggr <- aggr %>% NormalizeData() %>% ScaleData()
markers <- FindAllMarkers(aggr,logfc.threshold = 0.5,min.pct = 0.4,only.pos = T)
markers <- markers %>% group_by(cluster) %>% slice_max(avg_log2FC,n=5)
aggr$cellstate <- factor(aggr$cellstate,levels = c("AM1","AM2","AM3"))

# ---- SupplFigS2E Markers of AM subpopulations ---- 
DotPlot(aggr, features = rev(markers$gene),group.by ="cellstate",cols = "RdBu")+FontSize(x.text = 8,y.text = 8)+ylab("")+xlab("")+RotatedAxis()+theme(legend.title=element_blank(),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size = 10))+ coord_flip()




Idents(aggr) <- "age"
SUB.old <- subset(aggr,idents="Old")
SUB <- subset(aggr,idents="Young")

celltypes <- as.character(unique(aggr$cellstate))
plot_relfreq_list <- lapply(celltypes, function(celltype) {
  y.lim <- max(get_relFreqs(aggr, grouping = "group") %>% mutate(value=round(value*100,digits = 2)) %>% filter(cluster==celltype) %>% pull(value))
  p <- ggarrange(plot_relFreqs(get_relFreqs(SUB, grouping = "group") %>% mutate(value=round(value*100,digits = 2)), celltype, colors.group, group_order, title = celltype, save = F)+ylim(0, y.lim+(0.05*y.lim))+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 90)),
                 plot_relFreqs(get_relFreqs(SUB.old, grouping = "group") %>% mutate(value=round(value*100,digits = 2)), celltype, colors.group, group_order, title = celltype, save = F)+ylim(0, y.lim+(0.05*y.lim))+theme_classic()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank(),axis.text.x = element_text(angle = 90)),ncol=1,nrow=2)  
  p
})
names(plot_relfreq_list) <- celltypes

# ---- SupplFigS2F Relative proportions of AM subpopulations ----
ggarrange(plot_relfreq_list[["AM1"]],plot_relfreq_list[["AM2"]],plot_relfreq_list[["AM3"]],nrow = 1,ncol = 3)

# Run differential abundance analysis using the edgeR model
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

# Write the differential abundance table 
contrasts <- list(c("PBS","D14"),c("PBS","D28"),c("PBS","D60"))
young_DA_2 <- lapply(X=contrasts,FUN = DA_edgeR,dataset=SUB )
old_DA_2 <- lapply(X=contrasts,FUN = DA_edgeR,dataset=SUB.old )
write.xlsx(list("young"=Reduce(rbind,young_DA_2),"old"=Reduce(rbind,old_DA_2)),file = "AM_DA.xlsx")




#Subsetting within alveolar capillary endothelial cells 
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
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'sCap')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(10))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Lrg1+ aCap')
DimPlot(object = aggr, reduction = 'umap',label=T)+NoLegend()
aggr$cellstate <- Idents(aggr)

#Subsetting within Venous EC 
DefaultAssay(aggr) <- "integrated"
Idents(aggr) <- "cellstate"
SUB <-subset(x = aggr, idents = c("Venous EC"))
DimPlot(object = SUB, reduction = 'umap',label = T)
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.6)
DimPlot(object = SUB, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = SUB, reduction = 'umap',group.by = "group",label=F,cols = colors.group)+NoLegend()

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(2))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'SV EC')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %ni% c(2))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'PV EC')
DimPlot(object = aggr, reduction = 'umap',label=T)+NoLegend()

aggr$cellstate <- Idents(aggr)
cellstate_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                     "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                     "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                     "Prolif. EC","sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap" ,"SV EC","PV EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                     "Mesothelial cells","AT1","AT2","Multiciliated cells")
aggr$cellstate <- factor(aggr$cellstate,levels = cellstate_order)
Idents(aggr) <- "cellstate"



#Capillary endothelial cellstates markers
Idents(aggr) <- "cellstate"
SUB <-subset(x = aggr, idents = c( "Prolif. EC","sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap" ,"SV EC","PV EC","Arterial EC"))
SUB$cellstate <- factor(SUB$cellstate,levels = c("sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap" , "Prolif. EC","SV EC","PV EC","Arterial EC"))
subcolors <- c("sCap"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414","Lrg1+ aCap"="#FFD92F","aCap"="#5e1e33","Prolif. EC"="#fee0d2","SV EC"="#ccc357" ,"PV EC"="#f2ae5a","Arterial EC"="#e36805")
SUB <- RunUMAP(SUB,dims=1:50)

# ---- Fig2D UMAP EC subpopulations ----
DimPlot(object = SUB,group.by = "cellstate", reduction = 'umap',cols = subcolors,pt.size = 0.2)+NoLegend()+NoAxes()+ggtitle("")


# PCEC markers heatmap
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
markers$cluster <- factor(markers$cluster,levels = c("sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Prolif. EC","SV EC","PV EC","Arterial EC"))
genes.to.plot <- markers %>% filter(gene%ni%ngenes) %>% group_by(cluster) %>% top_n(9,ratio) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)

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

# ---- SupplFigS3A Heatmap of EC subpopulations markers ----
pheatmap(t(data),annotation_row  = annotation,show_rownames = F,cluster_cols = F,cluster_rows = F,angle_col =90,fontsize_col = 10,
         annotation_colors = list(cellstate= subcolors,sample=colors.group))

genes.to.plot <- markers %>% filter(gene%ni%ngenes) %>% filter(cluster%in%c("sCap","Lrg1+ gCap","Lrg1+ aCap")) %>% group_by(cluster) %>% top_n(15,ratio) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)

library('biomaRt')
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- genes.to.plot$gene
df <- listAttributes(mart)
G_list <- getBM(filters= "mgi_symbol", attributes= c("description","external_gene_name"),values=genes,mart= mart) %>% mutate(gene=external_gene_name)

genes.to.plot <- merge(genes.to.plot,G_list,by.x="gene")
genes.to.plot <- genes.to.plot  %>% dplyr::select(cluster ,gene,description) %>% arrange(cluster  ,gene)
genes.to.plot$description <- gsub("\\[Source.*","",genes.to.plot$description)

# ---- Table 1 EC subpopulations markers ----
write.xlsx(genes.to.plot,"table1.xlsx")


# ---- Fig2F Vlnplot EC subpopulations markers ----   
VlnPlot(SUB,features=c("Lrg1"),group.by = "cellstate",pt.size = 0.1,cols =  subcolors,sort = F,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(SUB,features=c("Col15a1"),group.by = "cellstate",pt.size = 0.1,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(SUB,features=c("Vwa1"),group.by = "cellstate",pt.size = 0.1,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(SUB,features=c("Ednrb"),group.by = "cellstate",pt.size = 0.1,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(SUB,features=c("Aplnr"),group.by = "cellstate",pt.size = 0.1,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(SUB,features=c("Ackr1"),group.by = "cellstate",pt.size = 0.1,cols =  subcolors,sort = F,same.y.lims = T,y.max = 4.5)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())


# EC relative proportions
celltype_order <- c("sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Prolif. EC","SV EC","PV EC","Arterial EC")
metadata <- SUB@meta.data %>% filter(cellstate %in% celltype_order)
metadata$cellstate <- factor(metadata$cellstate, levels = celltype_order)
metadata = as.data.frame.matrix(table(metadata$group, metadata$cellstate))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$group = rownames(metadata)
metadata$group =factor(metadata$group, levels = group_order)
metadata = gather(metadata, cellstate, percentage, 'sCap':'Arterial EC')
metadata$cellstate <- factor(metadata$cellstate, levels =celltype_order)

# ---- Fig2E EC subpopulations relative proportions----
ggplot(metadata, aes(x = group, y = percentage, fill = cellstate)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = subcolors) +
  ggtitle("EC relative proportions")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+ggtitle("")+NoLegend()


#Comparison between SVEC markers in human and mouse fibrosis

aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
Idents(aggr) <- "cellstate"
SUB <- subset(aggr,idents=c('SV EC',"PV EC","gCap","Lrg1+ gCap","sCap"))
DefaultAssay(SUB) <- "RNA"
SUB <- SUB %>% NormalizeData() %>% ScaleData()
SUB$cellstate <- factor(SUB$cellstate,levels = c('SV EC',"PV EC","gCap","Lrg1+ gCap","sCap"))

compa.list <- c("Aplnr","Rgcc","Kit","Sox11",
                "Abcb1a","Col15a1","Filip1","Galnt15","Plvap","Ndrg1","Meox1",
                "Gja1","Hspg2","Lamb1","Nr5a2","Spry1","Vwa1",
                "Ackr1","Amigo2","Hdac9","Il1r1","Lrrc1","Mmp16","Nrp2","Prss23","Selp", "Slc6a2")

# ---- Fig2G Capillary, Systemic and Venous markers in mouse ----
DotPlot(SUB, features = compa.list,group.by ="cellstate",cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)+ylab("")+xlab("")+theme(legend.title=element_blank(),legend.key.size = unit(0.3, 'cm'),legend.text = element_text(size = 8))



#Seurat object of the Habermann et al. dataset reintegrated and reannoted
Habermann <- readRDS("/data/truchi_data/Rdata_objects/Habermann_integrated.rds")
Habermann <- subset(Habermann,idents=c('SV EC',"aCap","gCap","PV EC","Arterial EC"))
DefaultAssay(Habermann) <- "integrated"
Habermann <- RunUMAP(Habermann, dims = 1:80)
Habermann <- FindNeighbors(Habermann,  dims = 1:80,k.param = 10)
Habermann <- FindClusters(Habermann, resolution = 1)
DimPlot(object = Habermann, reduction = 'umap',label=T)+NoLegend()
DimPlot(object = Habermann, reduction = 'umap',label=T,group.by = "celltype_MT")+NoLegend()
DimPlot(object = Habermann, reduction = 'umap',group.by = "Diagnosis")+NoLegend()

DefaultAssay(Habermann) <- "RNA" 
Habermann <- Habermann %>% NormalizeData() %>% ScaleData()
mm <- FindAllMarkers(Habermann,min.pct = 0.25,logfc.threshold = 0.5,only.pos = T)
sCapvsSVEC <- FindMarkers(Habermann,ident.1 = "4",ident.2 = "3",min.pct = 0.25,logfc.threshold = 0.5)
FeaturePlot(Habermann,features = c("ACKR1","SELP","BNC2","LRRC1"))
FeaturePlot(Habermann,features = c("KDR","RGCC","BTNL9","APLNR"))
FeaturePlot(Habermann,features = c("VWA1","COL15A1","PLVAP","VWF"))
FeaturePlot(Habermann,features = c("LRG1","SLC6A2","MEOX1","MEOX2"))

select.cells <- names(Idents(Habermann)[which(Idents(Habermann) %in% c(9))])
Idents(Habermann) <- "celltype_MT"
Habermann <- SetIdent(object = Habermann, cells = select.cells, value = 'sCap')

# ---- Fig2I UMAP EC Habermann IPF datataset----
DimPlot(object = Habermann, reduction = 'umap',label=F,cols = subcolors)+NoLegend()+NoAxes()+ggtitle("")

Habermann$celltype_MT <- Idents(Habermann)
Habermann.order <- c("sCap","gCap","aCap",'SV EC',"PV EC","Arterial EC")
Habermann$celltype_MT <- factor(Habermann$celltype_MT,levels =Habermann.order)

# ---- Fig2J VlnPlot EC markers Habermann IPF ----
VlnPlot(Habermann,features=c("LRG1"),group.by = "celltype_MT",pt.size = 0.1,sort = F,y.max = 4.5,cols = subcolors)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(Habermann,features=c("COL15A1"),group.by = "celltype_MT",pt.size = 0.1,sort = F,y.max = 4.5,cols = subcolors)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(Habermann,features=c("VWA1"),group.by = "celltype_MT",pt.size = 0.1,sort = F,y.max = 4.5,cols = subcolors)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(Habermann,features=c("EDNRB"),group.by = "celltype_MT",pt.size = 0.1,sort = F,y.max = 4.5,cols = subcolors)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(Habermann,features=c("APLNR"),group.by = "celltype_MT",pt.size = 0.1,sort = F,y.max = 4.5,cols = subcolors)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())
VlnPlot(Habermann,features=c("ACKR1"),group.by = "celltype_MT",pt.size = 0.1,sort = F,y.max = 4.5,cols = subcolors)+FontSize(x.text = 0,y.text = 8)+NoLegend()+ggtitle("")+theme(axis.title.x = element_blank(),axis.title.y = element_blank())


# Capillary, systemic and venous markers in IPF
Idents(Habermann) <- "celltype_MT"
SUB <- subset(Habermann,idents=c('SV EC',"PV EC","gCap","sCap"))
DefaultAssay(SUB) <- "RNA"
SUB <- SUB %>% NormalizeData() %>% ScaleData()
SUB$celltype_MT <- factor(SUB$celltype_MT,levels = c('SV EC',"PV EC","gCap","sCap"))

compa.list.hs <- c(str_to_upper(c("Aplnr","Rgcc","Kdr","Adgrf5")),
                   "ABCB1",str_to_upper(c("Col15a1","Filip1","Galnt15","Plvap","Ndrg1","Meox1",
                                          "Gja1","Hspg2","Lamb1","Nr5a2","Spry1","Vwa1",
                                          "Ackr1","Cpxm2","Hdac9","Il1r1","Lrrc1","Mmp16","Nrp2","Prss23","Selp", "Znf385d")))

# ---- Fig2K Capillary, Systemic and Venous markers in IPF ----
DotPlot(SUB, features = compa.list.hs,group.by ="celltype_MT",cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)+ylab("")+xlab("")+theme(legend.title=element_blank(),legend.key.size = unit(0.28, 'cm'),legend.text = element_text(size = 8))


SUB <- subset(Habermann,idents=c('Arterial EC','SV EC',"PV EC","gCap","aCap","sCap"))
metadata <- SUB@meta.data 
metadata$celltype_MT <- factor(metadata$celltype_MT, levels = Habermann.order)
metadata = as.data.frame.matrix(table(metadata$Diagnosis, metadata$celltype_MT))
metadata = as.data.frame(metadata / rowSums(metadata))
metadata$group = rownames(metadata)
metadata = gather(metadata, celltype_MT, percentage, 'sCap':'Arterial EC')
metadata$celltype_MT <- factor(metadata$celltype_MT, levels =Habermann.order)

# ---- Fig2I EC proportion Habermann IPF ----
ggplot(metadata, aes(x = group, y = percentage, fill = celltype_MT)) +
  geom_bar(stat = 'identity', color = "grey30") +
  scale_fill_manual(values = subcolors) +
  ggtitle("CEC relative proportions")+
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        axis.line = element_line(colour = "grey50"),
        axis.ticks = element_line(colour = "grey50"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))+ggtitle("")+NoLegend()


# Global list of markers per population
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
colors.Bleo<- c("Prolif. AM"="#a6d8f5","AM1"="#7FCDBB","AM2"="#41B6C4","AM3"="#1D91C0","IM"="#225EA8","LCM"="#0e3378",
                "Prolif. DC"="#ffd1ec","cDC1"="#e670b4","cDC2"="#9c119c","Mature DC"="#960c5a","pDC"="#de439a","cMonocytes"="#66119c","ncMonocytes"="#282b5c","Neutrophils"="#4300de","Neutrophil-like Monocytes"="#38015c","Mast cells"="#9315ed","Platelets"="#de1212",
                "Prolif. T cells"="#F7FCB9","CD4 T cells"="#D9F0A3","Reg T cells"="#ADDD8E","T helper cells"="#78C679","Gamma delta T cells"="#41AB5D","ILC2"="#739154","CD8 T cells"="#006837" ,"NKT1"="#004529","NKT2"="#01736b","NK cells"="#46c798","B cells"="#ebd513","Plasmocytes"="#ccbf47",
                "sCap"="#A6D854","Lrg1+ gCap"="#fa3737","gCap"="#911414","Lrg1+ aCap"="#FFD92F","aCap"="#5e1e33","Prolif. EC"="#fee0d2","SV EC"="#ccc357","PV EC"="#f2ae5a","Arterial EC"="#e36805","Lymphatic EC"="#7d3f15","Fibroblasts"="#d3d5db","Pericytes"="#806c6c",
                "Mesothelial cells"="#120c04","AT1"="#8184b8","AT2"="#3e447a","Multiciliated cells"="#1486f7")
cellstate_order <- c("Prolif. AM","AM1","AM2","AM3","IM","LCM",
                     "Prolif. DC","cDC1","cDC2","Mature DC","pDC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Neutrophils","Mast cells","Platelets",
                     "Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes",
                     "Prolif. EC","sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap" ,"SV EC","PV EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes",
                     "Mesothelial cells","AT1","AT2","Multiciliated cells")


DefaultAssay(aggr) <- "RNA"
aggr <- aggr %>% NormalizeData()
ngenes <- rownames(aggr)[which(rownames(aggr)%ni%unlist(findnoisygenes.mm(aggr)))]
Idents(aggr) <- "cellstate"
markers <- FindAllMarkers(aggr,features = ngenes,logfc.threshold = 0.5,min.pct = 0.4,only.pos = T)

# ---- SupplTableS2 Celltypes markers ---- 
write.xlsx(markers,file = "Supplementary_table_S2.xlsx")


# Lrg1 expression in all celltypes at all time points
aggr$celltype_global <- ifelse(aggr$cellstate%in%c("Prolif. DC","cDC1","cDC2","Mature DC","pDC"),"DCs",
                               ifelse(aggr$cellstate%in%c("cMonocytes","ncMonocytes","Neutrophil-like Monocytes"),"Monocytes",
                                      ifelse(aggr$cellstate%in%c("Prolif. T cells","CD4 T cells","Reg T cells","T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes"),"Lymphocytes",as.character(aggr$cellstate))))   
aggr$celltype_global <-factor(aggr$celltype_global,levels = c("sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","Prolif. EC","aCap","SV EC","PV EC","Arterial EC","Lymphatic EC","Fibroblasts","Pericytes","Mesothelial cells","AT1","AT2","Multiciliated cells",
                                                              "Prolif. AM","AM1","AM2","AM3","IM","LCM","DCs","Monocytes","Mast cells","Neutrophils","Platelets","Lymphocytes"))

Idents(aggr) <- "group"
DefaultAssay(aggr) <- "RNA"
SUB <-subset(x = aggr, subset = group=="PBS")
SUB <- SUB %>% NormalizeData()
SUB1 <-subset(x = aggr, idents = "D14")
SUB1 <- SUB1 %>% NormalizeData()
SUB2 <-subset(x = aggr, idents = "D28")
SUB2 <- SUB2 %>% NormalizeData()
SUB3 <-subset(x = aggr, idents = "D60")
SUB3 <- SUB3 %>% NormalizeData()


# ---- SupplFigS3B Lrg1 expression per condition in all celltypes ----
VlnPlot(SUB,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)

VlnPlot(SUB1,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)

VlnPlot(SUB2,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)

VlnPlot(SUB3,features = "Lrg1",group.by = "celltype_global",cols = colors.Bleo,pt.size = 0,y.max = 5)+NoLegend()+ylab("")+xlab("")+ggtitle("")+
  geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 2,jitter.width=6))+FontSize(x.text = 10,y.text = 7)

# LRG1 Expression per time point, age, and EC subpopulation
EC_order <- c("sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","SV EC","PV EC", "Arterial EC")
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
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
groups = c("sCap_Young","sCap_Old","Lrg1+gCap_Young","Lrg1+gCap_Old","gCap_Young","gCap_Old","Lrg1+aCap_Young","Lrg1+aCap_Old","aCap_Young","aCap_Old","SVEC_Young","SVEC_Old","PVEC_Young","PVEC_Old","ArterialEC_Young","ArterialEC_Old")

sub1$group <- factor(sub1$group,levels = groups)


# ---- SupplFigS3E Lrg1 expression per time point, age, and EC subpopulation ----
ggplot(sub1,aes(x=group,y=Condition,size=percentage)) + 
  geom_point(aes(col=Lrg1_expr),alpha=1) +
  scale_colour_gradientn(colours = c("#67001F", "#B2182B" ,"#D6604D" ,"#F4A582", "#FDDBC7", "#F7F7F7" ,"#D1E5F0", "#92C5DE" ,"#4393C3"),values = scales::rescale(c(max(sub1$Lrg1_expr),2,1.5,1,0.5, 0, -0.5, -1,min(sub1$Lrg1_expr))))+
  labs(y="time point")+theme_classic()+rotate_x_text(45)+ylab("")+xlab("")+theme(legend.title=element_blank(),legend.key.size = unit(0.4, 'cm'),legend.text = element_text(size = 8))


# Lrg1 Expression in capillary EC from Strunz et al.

# Seurat object of the Strunz et al. dataset reintegrated and reannoted
Strunz <- readRDS(file = "~/Strunz_RPCA.rds")
SUB <-subset(x = Strunz, idents = c("9","18"))
DefaultAssay(SUB) <- "integrated" 
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.5)
DimPlot(SUB)
SUB <-subset(x = SUB, idents = c("3","0","2","5"))
SUB <- RunUMAP(SUB, dims = 1:80)
SUB <- FindNeighbors(SUB,  dims = 1:80,k.param = 10)
SUB <- FindClusters(SUB, resolution = 0.5)

# Project annotated dataset on this one
ref <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
Idents(ref) <- "cellstate"
celltype_order <- c("sCap","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","SV EC","PV EC", "Arterial EC", "Lymphatic EC")
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
FeaturePlot(SUB,features = "Gja5")
ngenes <- rownames(SUB)[which(rownames(SUB)%ni%unlist(findnoisygenes.mm(SUB)))]
markers <- FindAllMarkers(SUB,features = ngenes,logfc.threshold = 0.4,min.pct = 0.4,only.pos = T)
markers$ratio <- markers$pct.1/markers$pct.2*markers$avg_log2FC
genes.to.plot <- markers %>% group_by(cluster) %>% top_n(50,ratio) %>% dplyr::arrange(desc(avg_log2FC),.by_group = TRUE)
VlnPlot(SUB,features="nCount_RNA",cols =  colors.Bleo,pt.size = 0)

select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(0))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'gCap')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(1))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'macrovascular EC')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(2))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'aCap')
select.cells <- names(Idents(SUB)[which(Idents(SUB) %in% c(3))])
SUB <- SetIdent(object = SUB, cells = select.cells, value = 'Lymphatic EC')
DimPlot(object = SUB, reduction = 'umap',label = T)

SUB$grouping <- factor(SUB$grouping, levels = c("PBS","d3","d7","d10","d14","d21","d28"))
aCap <- subset(SUB,ident=c('aCap'))
gCap <- subset(SUB,ident=c('gCap'))


# ---- SupplFigS3F Lrg1 Expression in capillary EC from Strunz ----
VlnPlot(gCap,features=c("Lrg1"),group.by = "grouping",pt.size = 0,cols =  colors.sample,sort = F)+FontSize(x.text = 10,y.text = 8)+geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=1.5))+FontSize(x.text = 8,y.text = 8,y.title = 0,x.title=0,main = 0)+NoLegend()
VlnPlot(aCap,features=c("Lrg1"),group.by = "grouping",pt.size = 0,cols =  colors.sample,sort = F)+geom_jitter(shape=10,size=0.1, position=position_jitterdodge(seed = 1, dodge.width = 0.9,jitter.width=1.5))+FontSize(x.text = 8,y.text = 8,y.title = 0,x.title=0,main = 0)+NoLegend()

