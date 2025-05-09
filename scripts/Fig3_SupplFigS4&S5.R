# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to annotate Lrg1+ PCEC functions and produce the plots of figure 3 and Supplementary figure S4 & S5
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")


# Differential Expression Analyses using DESeq2 on pseudo-bulks for Lrg1+ PCEC subpopulations

#For Lrg1+ gCap and Lrg1+ aCap
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap","Lrg1+ gCap","aCap","gCap") & group %ni% c("D60"))
DefaultAssay(aggr) <- "RNA"
aggr$samples <- aggr$mice

young <- subset(x = aggr, subset = age == c("Young"))
old <- subset(x = aggr, subset = age == c("Old"))

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  #rownames(x) <- gsub('.*\\.', '', rownames(x))  
  t(x)})

cts <- cts.split$`gCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.gCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")

write.xlsx(DEG.psdblk.gCap.mice,file = "Lrg1pos_gCap_pseudobulk.xlsx")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% filter(padj<0.05) 

cts <- cts.split$`aCap`
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.aCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.aCap.mice,file = "Lrg1pos_aCap_pseudobulk.xlsx")
DEG.psdblk.gCap.mice <- DEG.psdblk.gCap.mice %>% filter(padj<0.05) 


#For sCap (merge sCap and PVEC from PBS to have enough cells for pseudobulks)
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
sCap <-subset(x = aggr, subset = cellstate %in% c("sCap") & group %ni% c("PBS"))
gCap <-subset(x = aggr, subset = cellstate %in% c("gCap") & group %in% c("PBS"))
aggr <- merge(sCap,gCap)
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DEG.psdblk.sCap.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.sCap.mice,file = "sCap_pseudobulk.xlsx")
DEG.psdblk.sCap.mice <- DEG.psdblk.sCap.mice %>% filter(padj<0.05) 


#For SVEC (merge SVEC and PVEC from PBS to have enough cells for pseudobulks)
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
SVEC <-subset(x = aggr, subset = cellstate %in% c("SV EC")& group %ni% c("PBS"))
PVEC <-subset(x = aggr, subset = cellstate %in% c("PV EC") & group %in% c("PBS"))
aggr <- merge(SVEC,PVEC)
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DEG.psdblk.SVEC.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.SVEC.mice,file = "SVEC_pseudobulk.xlsx")
DEG.psdblk.SVEC.mice <- DEG.psdblk.SVEC.mice %>% filter(padj<0.05) 


#For PV EC
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("PV EC")& group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.PVEC.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.PVEC.mice,file = "PVEC_pseudobulk.xlsx")
DEG.psdblk.PVEC.mice <- DEG.psdblk.PVEC.mice %>% filter(padj<0.05)


#For Arterial EC
aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Arterial EC")& group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
plotDispEsts(dds)
resultsNames(dds)
DESeq2::plotMA(results(dds), ylim=c(-2,2))
DEG.psdblk.ArterialEC.mice <- results(dds, name = "condition_Bleo_vs_PBS")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_Bleo_vs_PBS")
write.xlsx(DEG.psdblk.ArterialEC.mice,file = "ArterialEC_pseudobulk.xlsx")
DEG.psdblk.ArterialEC.mice <- DEG.psdblk.ArterialEC.mice %>% filter(padj<0.05)


# From the differential analyses tables generated above, Functional enrichment analyses were performed using IPA

# Enrichment analysis outputs from IPA
# Upstream Regulators
UR <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/IPA_UR_EC.xlsx")
UR.order <- c("SPEN","MYC","CEBPB","EZH2","HIF1A","NFkB","STAT6","SMAD3",
              "NKX2-3","YAP1","SOX4","FOXM1","NFKBIA",
              "NFE2L2","SOX1","EPAS1","GATA6","FOXO1","KLF2")
UR$cellstate <- factor(UR$cellstate,levels = rev(c("sCap","Lrg1_gCap","Lrg1_aCap","SVEC","PVEC")))
UR$TF <- factor(UR$TF,levels=UR.order)

# ---- SupplFigS4F IPA Upstream Regulators ----
ggplot(UR, aes(y=cellstate, x = TF, color = zscore, size = pval)) + 
  geom_point() + 
  scale_colour_gradient2(low = "#024099",mid = "white",high = "#c45b10")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('')

# Canonical pathways
CP <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/IPA_CP_EC.xlsx")
CP.order <- rev(c("Extracellular matrix organization","TGF-β Signaling","Pulmonary Fibrosis Idiopathic Signaling Pathway",
                  "Activin Inhibin Signaling Pathway","Gap Junction Signaling","Ephrin Receptor Signaling","Smooth Muscle Contraction","RHO GTPase cycle","STAT3 Pathway","Signaling by VEGF",
                  "FAK Signaling","Beta-catenin independent WNT signaling","Interferon gamma signaling"))
CP$cellstate <- factor(CP$cellstate,levels = c("sCap","Lrg1_gCap","Lrg1_aCap","SVEC","PVEC"))
CP$CP <- factor(CP$CP,levels=CP.order)

# ---- Fig3A IPA Canonical pathways ----
ggplot(CP, aes(x=cellstate, y = CP, color = zscore, size = pval)) + 
  geom_point() + 
  scale_colour_gradient2(low = "#024099",mid = "white",high = "#c45b10")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('')


# Diseases and Molecular Functions
DF <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/IPA_DF_EC.xlsx")
DF.order <- rev(c("Sprouting","Branching of cells","Cell movement","Migration of cells","Fibrogenesis","Formation of cellular protrusions",
                  "Cell death of endothelial cells","Microtubule dynamics","Organization of cytoskeleton","Cell survival"))

DF$cellstate <- factor(DF$cellstate,levels = c("sCap","Lrg1_gCap","Lrg1_aCap","SVEC","PVEC"))
DF$DF <- factor(DF$DF,levels=DF.order)

# ---- Fig3B IPA Diseases and Molecular Functions ----
ggplot(DF, aes(x=cellstate, y = DF, color = zscore, size = pval)) + 
  geom_point() + 
  scale_colour_gradient2(low = "#024099",mid = "white",high = "#c45b10")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab('') + xlab('')



#Expression of Hypoxia markers
#Loading dataset
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
markers_hypox <- c("Aldoa","Ankrd37","Apln","App","Bhlhe40","Bnip3","Dnah11","Ero1l","Gadd45b","Gapdh","Gpi1","Hif1a","Idh2","Ldha","Pdgfb","Pkm","Plod1","Serpine1")
SUB <-subset(x = aggr, subset = cellstate %in% c("aCap","gCap","Lrg1+ aCap","Lrg1+ gCap","sCap"))
EC_order <- c("Lrg1+ aCap","aCap","Lrg1+ gCap","gCap","sCap")
SUB$cellstate <- factor(SUB$cellstate,levels = rev(EC_order))
SUB <- SUB %>% NormalizeData() 

# ---- Fig3C Hypoxia markers in PCEC ----
DotPlot(SUB, features = markers_hypox,group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+
  FontSize(x.text = 8,y.text = 8)+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+NoLegend()

#Expression of Pro-angiogenic markers
angiogenic.features <- c("Dll4","Esm1","Hspg2","Plxnd1","Rgcc","Aplnr","Smad1","Sparc","Col4a1","Col4a2","Flt1","Xbp1","Btg1","Rhoj","Cxcl12","Sox17","Adamts1","Adam15","Vegfa","Ackr3","Lrg1","Nrp1","Kdr","Hspb1","Emp2","Cd34","Serpine1","Unc5b","Hif1a","Tgfbr2")

# ---- Fig3D Angiogenic markers in PCEC ----
DotPlot(SUB, features = rev(angiogenic.features),group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+
  FontSize(x.text = 8,y.text = 8)+theme(axis.title.x = element_blank(),axis.title.y = element_blank())+NoLegend()

# NicheNet Analysis
library(nichenetr)
library(circlize)

# DATA PREPARATION

# Select top angiogenic associated genes for NicheNet analysis
gCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ gCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
aCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ aCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
SVEC.angiogenic.features <- plot[["data"]] %>% filter(id=="sCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()

# Load Seurat Object
aggr <- readRDS(file = "~/Truchi_et_al_seuratobj.rds")
DefaultAssay(aggr) <- "RNA"
select.cells <- names(Idents(aggr)[which(Idents(aggr)%in%c("Prolif. DC","cDC1","cDC2","pDC","Mature DC","cMonocytes","ncMonocytes","Neutrophil-like Monocytes","Platelets"))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Myeloid cells')
select.cells <- names(Idents(aggr)[which(Idents(aggr)%in%c("Prolif. T cells","CD4 T cells","Reg T cells", "T helper cells","Gamma delta T cells","ILC2","CD8 T cells","NKT1","NKT2","NK cells","B cells","Plasmocytes"))])
aggr <- SetIdent(object = aggr, cells = select.cells, value = 'Lymphocytes')
aggr$celltype <- Idents(aggr)


# Read in NicheNet's ligand-target prior model, ligand-receptor network and weighted integrated networks
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))

lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
head(lr_network)

weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks[[1]] <- weighted_networks[[1]] %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks[[2]] <- weighted_networks[[2]] %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to")) %>% drop_na()
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network

# Convert the NicheNet network gene symbols from human to mouse based on one-to-one orthology
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

 #Overview of ligand expression in each cell population
ligands <- data.frame(ligand=sort(unique(colnames(ligand_target_matrix))))
exp_ligands <- ligands[which(ligands$ligand%in%rownames(aggr)),]
receptors <- data.frame(receptors=sort(unique(lr_network$to)))
Lrg1_network <- sig_network %>% filter(from=="LRG1" | to=="LRG1")
ligand_exp_matrix <- aggr@assays[["RNA"]][match(exp_ligands,rownames(aggr@assays[["RNA"]])),]

# STEP1 : Define expressed genes in sender and receiver population 
## Define receiver cells
receiver = "Lrg1+ gCap"
expressed_genes_receiver = get_expressed_genes(receiver, aggr, pct = 0.25,assay_oi ="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cells
sender_celltypes <- c( "Prolif. AM","AM1","AM2","AM3","IM","LCM","Fibroblasts","Pericytes","sCap","Prolif. EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Arterial EC","SV EC","PV EC","Lymphatic EC","Mesothelial cells","AT1","AT2","Multiciliated cells","Mast cells","Neutrophils","Lymphocytes","Myeloid cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, aggr, 0.25,assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# STEP2 : Define the gene set of interest and a background of gene
# geneset= pro angiogenic features expressed in Lrg1+ gCap

SUB <-subset(x = aggr, idents = c("aCap","gCap","Lrg1+ aCap","Lrg1+ gCap","sCap"))
EC_order <- c("Lrg1+ aCap","aCap","Lrg1+ gCap","gCap","sCap")
SUB$cellstate <- factor(SUB$cellstate,levels = rev(EC_order))
SUB <- SUB %>% NormalizeData() %>% ScaleData()
# Expression of Pro-angiogenic markers
angiogenic.features <- c("Hspg2","Plxnd1","Rgcc","Aplnr","Smad1","Sparc","Col4a1","Col4a2","Flt1","Xbp1","Btg1","Rhoj","Cxcl12","Sox17","Adamts1","Adam15","Vegfa","Ackr3","Lrg1","Nrp1","Kdr","Hspb1","Emp2","Cd34","Serpine1","Unc5b","Hif1a","Tgfbr2")
plot <- DotPlot(SUB, features = angiogenic.features,group.by ="cellstate",dot.scale = 5,cols = "RdBu") + RotatedAxis()+FontSize(x.text = 8,y.text = 8)
plot

gCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ gCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
aCap.angiogenic.features <- plot[["data"]] %>% filter(id=="Lrg1+ aCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()
SVEC.angiogenic.features <- plot[["data"]] %>% filter(id=="sCap"& avg.exp >1 & pct.exp >33) %>% pull(features.plot) %>% as.character()

# STEP3 : Define a set of potential ligands 
# ligands that are expressed by the "sender/niche" cell population and bind a receptor expressed by the "receiver/target" population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expr_rec_gCap <- expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# STEP4 : Perform NicheNet's ligand activity analysis on the gene set of interest 
# INPUT :
# geneset = gene set of interest
# background_expressed_genes = all expressed genes in receiver cells
# ligand_target_matrix = ligand target matrix denoting regulatory potential scores
# potential_ligands = expressed ligands in sender cells

# OUTPUT : statistical metrics showing how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes
# pearson corr coeff > 0.1 = good candidate
ligand_activities = predict_ligand_activities(geneset = gCap.angiogenic.features, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
write.table(ligand_activities,file = "ligand_activity_Lrg1gCap.txt",row.names = F,sep = "\t")

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(14, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(14, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
Idents(aggr) <- "celltype"
Sub <- subset(aggr,idents = sender_celltypes)
Sub <- Sub %>% NormalizeData()
DotPlot(Sub, features = best_upstream_ligands %>% rev(), cols = c("#FFFFB2","#B10026")) + RotatedAxis()


# Same for Lrg1+ aCap
# STEP1 : Define expressed genes in sender and receiver population 
## Define receiver cells
receiver = "Lrg1+ aCap"
expressed_genes_receiver = get_expressed_genes(receiver, aggr, pct = 0.25,assay_oi ="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cells
sender_celltypes <- c( "Prolif. AM","AM1","AM2","AM3","IM","LCM","Fibroblasts","Pericytes","sCap","Prolif. EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","Arterial EC","SV EC","PV EC","Lymphatic EC","Mesothelial cells","AT1","AT2","Multiciliated cells","Mast cells","Neutrophils","Lymphocytes","Myeloid cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, aggr, 0.25,assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# STEP2 : Define the gene set of interest and a background of gene
# geneset= pro angiogenic features expressed in Lrg1+ gCap

aCap.angiogenic.features

# STEP3 : Define a set of potential ligands 
# ligands that are expressed by the "sender/niche" cell population and bind a receptor expressed by the "receiver/target" population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expr_rec_aCap <- expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# STEP4 : Perform NicheNet's ligand activity analysis on the gene set of interest 
# INPUT :
# geneset = gene set of interest
# background_expressed_genes = all expressed genes in receiver cells
# ligand_target_matrix = ligand target matrix denoting regulatory potential scores
# potential_ligands = expressed ligands in sender cells

# OUTPUT : statistical metrics showing how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes
# pearson corr coeff > 0.1 = good candidate
ligand_activities = predict_ligand_activities(geneset = aCap.angiogenic.features, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
write.table(ligand_activities,file = "ligand_activity_Lrg1aCap.txt",row.names = F,sep = "\t")

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(14, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(14, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
Idents(aggr) <- "celltype"
Sub <- subset(aggr,idents = sender_celltypes)
Sub <- Sub %>% NormalizeData()
DotPlot(Sub, features = best_upstream_ligands %>% rev(), cols = c("#FFFFB2","#B10026")) + RotatedAxis()


# Same for sCap
# STEP1 : Define expressed genes in sender and receiver population 
## Define receiver cells
receiver = "sCap"
expressed_genes_receiver = get_expressed_genes(receiver, aggr, pct = 0.25,assay_oi ="RNA")
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## Define sender cells
sender_celltypes <- c( "Prolif. AM","AM1","AM2","AM3","IM","LCM","Fibroblasts","Pericytes","sCap","Prolif. EC","Lrg1+ gCap","gCap","Lrg1+ aCap","aCap","SV EC","PV EC","Arterial EC","Lymphatic EC","Mesothelial cells","AT1","AT2","Multiciliated cells","Mast cells","Neutrophils","Lymphocytes","Myeloid cells")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, aggr, 0.25,assay_oi ="RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

# STEP2 : Define the gene set of interest and a background of gene
# geneset= pro angiogenic features expressed in sCap

SVEC.angiogenic.features

# STEP3 : Define a set of potential ligands 
# ligands that are expressed by the "sender/niche" cell population and bind a receptor expressed by the "receiver/target" population
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expr_rec_SVEC <- expressed_receptors
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

# STEP4 : Perform NicheNet's ligand activity analysis on the gene set of interest 
# INPUT :
# geneset = gene set of interest
# background_expressed_genes = all expressed genes in receiver cells
# ligand_target_matrix = ligand target matrix denoting regulatory potential scores
# potential_ligands = expressed ligands in sender cells

# OUTPUT : statistical metrics showing how well a ligand can predict the observed differentially expressed genes compared to the background of expressed genes
# pearson corr coeff > 0.1 = good candidate
ligand_activities = predict_ligand_activities(geneset = SVEC.angiogenic.features, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))
ligand_activities
write.table(ligand_activities,file = "ligand_activity_sCap.txt",row.names = F,sep = "\t")

# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(14, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity

best_upstream_ligands = ligand_activities %>% top_n(14, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
Idents(aggr) <- "celltype"
Sub <- subset(aggr,idents = sender_celltypes)
Sub <- Sub %>% NormalizeData()
DotPlot(Sub, features = best_upstream_ligands %>% rev(), cols = c("#FFFFB2","#B10026")) + RotatedAxis()


#Add best ligand activities for Lrg1+gCap receiver cells
ligand_activities_Lrg1gCap <- read.table( "ligand_activity_Lrg1gCap.txt",header = T) 
ligand_activities_Lrg1aCap <- read.table( "ligand_activity_Lrg1aCap.txt",header = T) 
ligand_activities_sCap <- read.table( "ligand_activity_sCap.txt",header = T) 

best_upstream_ligands <- unique(c(ligand_activities_sCap %>% top_n(14,pearson) %>% pull(test_ligand),ligand_activities_Lrg1aCap %>% top_n(14,pearson) %>% pull(test_ligand),ligand_activities_Lrg1gCap %>% top_n(14,pearson) %>% pull(test_ligand)))

# STEP5 : Infer target genes and receptors of top-ranked ligands and visualize in a heatmap
# ACTIVE TARGET GENE INFERENCE
# Which genes of the geneset are in the top 200 predicted targets of the top potential ligands
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = angiogenic.features, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
# matrix processing for the heatmap
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.1)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% sort() %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network


# RECEPTORS OF TOP-RANKED LIGANDS
# get the ligand-receptor network of the top-ranked ligands
expressed_receptors <- unique(c(expr_rec_aCap,expr_rec_gCap,expr_rec_sCap))
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
# convert to a matrix
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")

order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
order_receptors_large <- order_receptors

# RECEPTORS OF TOP-RANKED LIGANDS, but after considering only bona fide ligand-receptor interactions documented in literature and publicly available databases
# get only bona fide ligand-receptor interactions
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

# get the weights of the bona fide ligand-receptor interactions as used in the NicheNet model
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))
# convert to a matrix
lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_position="right",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict


# SUMMARY VISUALIZATIONS OF THE NICHENET ANALYSIS
# Figure to combine : 1) ligand activity; 2) ligand expression dotplot; 3) ligand logFC; 4) inference of ligand - target gene activity 
# ligand activity heatmap
vis_ligand_pearson = Reduce(function(x, y, ...) merge(x, y,all = TRUE,by="test_ligand"), 
                            list(ligand_activities_Lrg1gCap  %>% mutate(Lrg1.gCap=pearson )%>% dplyr::select(test_ligand,Lrg1.gCap),
                                 ligand_activities_Lrg1aCap  %>% mutate(Lrg1.aCap=pearson )%>% dplyr::select(test_ligand,Lrg1.aCap),
                                 ligand_activities_sCap  %>% mutate(sCap=pearson )%>% dplyr::select(test_ligand,sCap))) %>% filter(test_ligand %in% best_upstream_ligands) %>% column_to_rownames(.,"test_ligand") %>% as.matrix()
vis_ligand_pearson[is.na(vis_ligand_pearson)] <- min(vis_ligand_pearson[!is.na(vis_ligand_pearson)])
vis_ligand_pearson <- vis_ligand_pearson[match(rev(best_upstream_ligands),rownames(vis_ligand_pearson)),]

p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))
p_ligand_pearson


# ligand expression Seurat dotplot
order_ligands_adapted = order_ligands
order_ligands_adapted <- rev(c("Fgf2","Vegfa","Il1b","Tgfb1","Apoe","Col18a1","Hgf","Adam17","Bmp4","Cxcl12","Sema3e","Tnf","Pgf","Gpi1","Tnfsf10","Tgfb3","Il1a"))
aggr$celltype <- factor(aggr$celltype,levels = sender_celltypes)
aggr <- aggr %>% NormalizeData() %>% ScaleData()
rotated_dotplot = DotPlot(aggr %>% subset(celltype %in% sender_celltypes),group.by = "celltype", features = order_ligands_adapted , cols = "RdBu",assay = "RNA") + coord_flip() + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12)) # flip of coordinates necessary because we want to show ligands in the rows when combining all plots
figures_without_legend = cowplot::plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
  rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab("")+xlab(""),
  align = "hv",
  nrow = 1,
  rel_widths = c(ncol(vis_ligand_pearson)+6,  ncol(vis_ligand_target)))
legends = cowplot::plot_grid(
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_pearson)),
  ggpubr::as_ggplot(ggpubr::get_legend(rotated_dotplot)),
  ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)),
  nrow = 1,
  align = "h", rel_widths = c(1.5, 1, 1))

combined_plot = cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")
combined_plot

# ---- SupplFigS5A Ligand activity ----
p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank(),axis.title.x = element_text()) 

# ---- SupplFigS5A Ligand expression ----
rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("") + xlab("") + scale_y_discrete(position = "right")

# ---- SupplFigS5A Predicted targets ----
p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab("")+xlab("")

# Checking for best upstream regulators + bona fide receptors + top predicted genes expression value, fraction of expressing cells and differential expression
# ligand expression Seurat dotplot
sub <- subset(aggr,idents = c("Lrg1+ gCap","Lrg1+ aCap","sCap")) 
sub <- sub %>% NormalizeData() %>% ScaleData()

vis_ligand_receptor_network_final <- vis_ligand_receptor_network[,c("Apoe","Hgf")]%>% as.data.frame()
vis_ligand_receptor_network_final <- t(apply(vis_ligand_receptor_network_final, 1, function(x) replace(x, x != max(x, na.rm = TRUE), 0))) %>% as.data.frame()
vis_ligand_receptor_network_final <- vis_ligand_receptor_network_final[c(rownames(vis_ligand_receptor_network_strict),"Scarb1"),]
a <- rbind(vis_ligand_receptor_network_strict,0)
rownames(a) <- c(rownames(vis_ligand_receptor_network_strict),"Scarb1")
vis_ligand_receptor_network_final <- cbind(vis_ligand_receptor_network_final,a)
vis_ligand_receptor_network_final<-vis_ligand_receptor_network_final[match(c("Fgfr1","Fgfr3","Sdc4","Flt1","Kdr","Nrp1","Nrp2","Il1r1","Acvrl1","Tgfbr2","Tgfbr3","Scarb1","Itgb1","Itga5","Bmpr2","Ackr3","Plxnd1","Tnfrsf1a","Ltbr","Amfr","Tnfrsf11b"),rownames(vis_ligand_receptor_network_final)),
                                                                     match(rev(best_upstream_ligands[which(best_upstream_ligands%in%colnames(vis_ligand_receptor_network_final))]),colnames(vis_ligand_receptor_network_final))]


p_ligand_receptor_network_final = vis_ligand_receptor_network_final %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_position="right",legend_title = "Prior interaction potential")

# ---- Fig3F Ligand-receptor affinity ----
p_ligand_receptor_network_final+ ylab("")+xlab("")+NoLegend()

# ---- Fig3F receptor expression in Lrg1+ PCEC ----
DotPlot(sub, features = rownames(vis_ligand_receptor_network_final),scale=F, cols = c("white","#730d2b"),dot.scale = 4)  + theme(legend.text = element_text(size = 10), legend.title = element_text(size = 12))+ theme(legend.position = "right", axis.ticks = element_blank(), axis.title.x = element_text(size = 12), axis.text.y = element_text(face = "italic", size = 9), axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) + ylab("") + xlab("") + scale_y_discrete(position = "left")+ scale_x_discrete(position = "top")+NoLegend()


#-.- FOLLOW-UP ANALYSIS : Circos plot visualization to show active ligand-target links between interacting cells -.-#

best_upstream_ligands <- colnames(vis_ligand_receptor_network_final)
best_upstream_receptors <- rownames(vis_ligand_receptor_network_final)

Mast_cell_specific_ligands = c("Hgf")
Epithelial_specific_ligands = c("Sema3e","Bmp4")
Neutrophil_specific_ligands = c("Il1b")
EC_specific_ligands = c("Cxcl12", "Vegfa","Tnfsf10")
Macrophages_specific_ligands =c("Apoe","Tnf","Adam17","Il1a")
Mesenchymal_specific_ligands =c("Fgf2","Pgf","Col18a1","Tgfb3")
general_ligands = setdiff(best_upstream_ligands,c(EC_specific_ligands,Epithelial_specific_ligands,Macrophages_specific_ligands,Mesenchymal_specific_ligands,Mast_cell_specific_ligands,Neutrophil_specific_ligands)) %>% print()

ligand_type_indication_df = tibble(
  ligand_type = c(rep("EC-specific", times = EC_specific_ligands %>% length()),
                  rep("Mesenchymal-specific", times = Mesenchymal_specific_ligands %>% length()),
                  rep("Epithelial-specific", times = Epithelial_specific_ligands  %>% length()),
                  rep("General", times = general_ligands %>% length()),
                  rep("Macrophages-specific", times = Macrophages_specific_ligands  %>% length()),
                  rep("Neutrophil-specific", times = Neutrophil_specific_ligands %>% length()),
                  rep("Mast_cell-specific", times = Mast_cell_specific_ligands  %>% length())),
  ligand = c(EC_specific_ligands,Mesenchymal_specific_ligands,Epithelial_specific_ligands, general_ligands,Macrophages_specific_ligands, Neutrophil_specific_ligands,Mast_cell_specific_ligands))

#-.- Visualize ligand-receptor interactions of the prioritized ligands in a circos plot -.-#

# RECEPTORS OF TOP-RANKED LIGANDS

grid_col_ligand =c("General" = "#109e52","Mesenchymal-specific" = "#7d7a78","Macrophages-specific" = "#1D91C0","EC-specific" = "#eda43e","Neutrophil-specific" = "#4300de","Epithelial-specific" = "#2e2f5c","Mast_cell-specific" = "#9315ed")
grid_col_target =c("geneset" = "tomato")

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% dplyr::rename(ligand = from, receptor = to)

lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "receptor") %>% inner_join(ligand_type_indication_df)
grid_col_ligand 
grid_col_receptor =c( "receptor" = "#a37b7a")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
#show only links with a weight higher than a predefined cutoff: links belonging to the 66% of lowest scores were removed.
cutoff_include_all_ligands = circos_links$weight %>% quantile(0.33)
circos_links = circos_links %>% filter(weight>cutoff_include_all_ligands)
links_circle = circos_links %>% dplyr::select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# Prepare the circos visualization: order ligands and receptors
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(EC_specific_ligands,Mesenchymal_specific_ligands,Epithelial_specific_ligands,general_ligands,Macrophages_specific_ligands,Neutrophil_specific_ligands,Mast_cell_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)

order = c(ligand_order,receptor_order)

# Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "EC-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mesenchymal-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Epithelial-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "General") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Macrophages-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Neutrophil-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "Mast_cell-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "receptor") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)

# ---- Fig3G Circosplot ----
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,link.sort = TRUE, order = order, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()


# Differentially expressed receptors and key genes of Angiogenic pathways
library(DESeq2)
genelist <- c("Fgfr1","Fgfr3","Sdc4","Flt1","Kdr","Nrp1","Nrp2","Il1r1","Acvrl1","Tgfbr2","Tgfbr3","Scarb1","Itgb1","Itga5","Bmpr2","Ackr3","Plxnd1","Tnfrsf1a","Ltbr","Amfr","Tnfrsf11b")
deg_receptors <- read.xlsx("~/Supplementary_table_S3.xlsx") %>% filter(celltype %in% c("gCap","aCap","Venous EC")  & gene%in%genelist & baseMean > 3)

aggr <- readRDS(file = "/data/truchi_data/Rdata_objects/Truchi_et_al_seuratobj.rds")
aggr <-subset(x = aggr, subset = cellstate %in% c("Lrg1+ aCap", "aCap","Lrg1+ gCap", "gCap","PV EC") & group %ni% c("D60"))
Idents(aggr) <- "RNA"
aggr$samples <- aggr$mice

cts <- AggregateExpression(aggr,group.by = c("celltype","samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- t(cts$RNA) %>% as.data.frame()
splitRows <- gsub('_.*', '', rownames(cts))
cts.split <- split.data.frame(cts, f = factor(splitRows))
cts.split <- lapply(cts.split, function(x){
  t(x)})

cts <- cts.split[["gCap"]]
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds_gCap <- DESeq(dds,test="Wald")

cts <- cts.split[["aCap"]]
colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('PBS', samples), 'PBS', 'Bleo')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('PBS', 'Bleo'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds_aCap <- DESeq(dds,test="Wald")


DEG_gCap <- read.xlsx("~/DEG.pseudobulks/gCap_pseudobulk.xlsx")
DEG_aCap <- read.xlsx("~/DEG.pseudobulks/aCap_pseudobulk.xlsx")

# ---- SupplFigS5B Differentially expressed receptors ----
for (gene in unique(deg_receptors$gene)) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('PBS', sample), 'PBS', 'Bleo'))
  cts$condition <- factor(cts$condition,levels = c("PBS","Bleo") )
  p <- ggboxplot(cts, x="condition",y= gene,fill = "condition",ylab = gene, palette = c("PBS"=gg_color_hue(2)[2],"Bleo"=gg_color_hue(2)[1]),outlier.shape = NA)+geom_jitter(shape=16, position=position_jitter(0.2),size=1.5)
  stat.test <- tibble(group1=c("PBS"),group2=c("Bleo"),p.adj=c(round(DEG_gCap$padj[which(DEG_gCap$gene==gene)],digits = 3)),y.position=max(counts(dds_gCap, normalized=TRUE)[gene,])+1)
  p + stat_pvalue_manual(stat.test, size = 2.6,label = "p.adj = {p.adj}")+FontSize(x.text = 10,y.text = 10,x.title = 0,y.title = 11)+NoLegend()
  ggsave(paste0("SupplFigS5_B_",gene,"_boxplot.pdf"), width = 3.5, height = 5.8, units = c("cm"), dpi = 200)
}

# ---- SupplFigS5B Differentially expressed angiogenesis associated genes ----
genes <- c("Mmp14",
           "Eng","Acvrl1","Smad6","Id1",
           "Tgfbr1","Smad1","Smad2","Smad5")
for (gene in genes) {
  cts <- data.frame(counts(dds_gCap, normalized=TRUE))[gene,] %>% t() %>% as.data.frame() %>% mutate(sample=rownames(.)) %>% 
    mutate(condition=ifelse(grepl('PBS', sample), 'PBS', 'Bleo'))
  cts$condition <- factor(cts$condition,levels = c("PBS","Bleo") )
  p <- ggboxplot(cts, x="condition",y= gene,fill = "condition",ylab = gene, palette = c("PBS"=gg_color_hue(2)[2],"Bleo"=gg_color_hue(2)[1]),outlier.shape = NA)+geom_jitter(shape=16, position=position_jitter(0.2),size=1.5)
  stat.test <- tibble(group1=c("PBS"),group2=c("Bleo"),p.adj=c(round(DEG_gCap$padj[which(DEG_gCap$gene==gene)],digits = 3)),y.position=max(counts(dds_gCap, normalized=TRUE)[gene,])+1)
  p + stat_pvalue_manual(stat.test, size = 2.6,label = "p.adj = {p.adj}")+FontSize(x.text = 10,y.text = 10,x.title = 0,y.title = 11)+NoLegend()
  ggsave(paste0("SupplFigS5_C_",gene,"_boxplot.pdf"), width = 3.5, height = 5.8, units = c("cm"), dpi = 200)
}


# Compare PCEC signatures induced by bleomycin or by IPF

SUB <- readRDS("/data/truchi_data/Rdata_objects/Habermann_EC.rds")
SUB$samples <- paste0(SUB$Diagnosis,".",SUB$Sample_Name )
SUB <- subset(SUB,subset=samples %in% c("Control.THD0002","Control.THD0005","Control.VUHD65","Control.VUHD67","Control.VUHD68","Control.VUHD70",
                                        "IPF.TILD006","IPF.TILD015","IPF.TILD028","IPF.VUILD59","IPF.VUILD60","IPF.VUILD61","IPF.VUILD63"))

DimPlot(SUB)

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


#For sCap (merge sCap and gCap from controls to have enough cells for pseudobulks)
SUB$id <- paste0(SUB$celltype_MT,"-",SUB$Diagnosis)
SUB_bis <- subset(SUB,subset= id%in% c("sCap-Control","gCap-Control","sCap-IPF"))
SUB_bis$samples <- paste0(SUB_bis$Diagnosis,".",SUB_bis$orig.ident)
DefaultAssay(SUB_bis) <- "RNA" 

cts <- AggregateExpression(SUB_bis,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
DEG.psdblk.sCap <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")

#For SVEC (merge SVEC and PVEC from controls to have enough cells for pseudobulks)
SUB_bis <- subset(SUB,subset= id%in% c("PV EC-Control","SV EC-Control","SV EC-IPF"))
SUB_bis$samples <- paste0(SUB_bis$Diagnosis,".",SUB_bis$orig.ident)
DefaultAssay(SUB_bis) <- "RNA" 

cts <- AggregateExpression(SUB_bis,group.by = c("samples"),return.seurat = F,assays = "RNA",slot = "counts")
cts <- cts$RNA %>% as.data.frame()

colData <- data.frame(samples = colnames(cts)) %>% mutate(condition = ifelse(grepl('Control', samples), 'Control', 'IPF')) %>% column_to_rownames(var = 'samples')
colData$condition <- factor(colData$condition,levels = c('Control', 'IPF'))
dds <- DESeqDataSetFromMatrix(countData = cts,colData = colData,design = ~ condition)
dds <- dds[rowSums(counts(dds)) >=10,]
dds <- DESeq(dds,test="Wald")
DEG.psdblk.SVEC <- results(dds, name = "condition_IPF_vs_Control")%>% as.data.frame()%>% mutate(gene=rownames(.))%>% mutate(test="condition_IPF_vs_Control")



# ---- SupplTableS4 ---- 
# Create Supplementary_table_S4 : Differential expression analysis (DESeq2) between IPF patients and healthy control for aCap, gCap, SV EC, PV EC and Arterial EC from the Habermann et al. Dataset
write.xlsx(list("aCap"=DEG.psdblk.aCap,"gCap"=DEG.psdblk.gCap,"sCap"=DEG.psdblk.sCap,"SVEC"=DEG.psdblk.SVEC,"PVEC"=DEG.psdblk.PVEC,"ArterialEC"=DEG.psdblk.ArterialEC),file ="Supplementary_table_S4.xlsx")



library(nichenetr)
# To compare PCEC signatures induced by bleomycin or by IPF, load differential expression tables obtained for Lrg1+PCEC
DEG.psdblk.sCap <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/Supplementary_table_S4.xlsx",sheet = "sCap")
sCap_IPF_1 <- DEG.psdblk.sCap %>% filter(padj < 0.05) %>% pull(gene) 
sCap_Bleo_1 <- read.xlsx("sCap_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
sCap_Bleo_1 <- sCap_Bleo_1[complete.cases(sCap_Bleo_1)]

intersect_sCap <- intersect(sCap_Bleo_1,sCap_IPF_1) %>% unique()

sCap_IPF <- DEG.psdblk.sCap %>% filter(gene %in% intersect_sCap) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
sCap_Bleo <- read.xlsx("sCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
sCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(sCap_Bleo$gene),Bleo=sCap_Bleo$Bleo) %>% filter(gene %in% intersect_sCap)
sCap <- merge(sCap_Bleo,sCap_IPF)
sCap$comparison <- sign(sCap$Bleo) == sign(sCap$IPF)
sCap <- sCap %>% filter(comparison == TRUE)


#Numbers for Venn Diagramm of figure 5E
nrow(sCap)
length(sCap_IPF_1)-nrow(sCap)
length(sCap_Bleo_1)-nrow(sCap)
total <- sum(nrow(sCap),c(length(sCap_IPF_1)-nrow(sCap)),c(length(sCap_Bleo_1)-nrow(sCap)))

nrow(sCap)/total*100
c(length(sCap_IPF_1)-nrow(sCap))/total*100
c(length(sCap_Bleo_1)-nrow(sCap))/total*100

df <- sCap
df$Bleo[df$Bleo > 4 ] <- 4
df$Bleo[df$Bleo < -4 ] <- -4
df$IPF[df$IPF > 4 ] <- 4
df$IPF[df$IPF < -4 ] <- -4
df <- df %>% mutate(Avg=rowMeans(cbind(IPF,Bleo))) %>% mutate(sd=rowSds(cbind(IPF,Bleo)))
top_up <- df %>% filter(sd < 2 & Avg > 0) %>% slice_max(Avg,n=75)%>% filter(gene %ni% c("ABCB1","RPLP0","RPL12","RPL11", "RPS2")) %>% pull(gene)
top_down <- df %>% filter(sd < 2 & Avg < 0) %>% slice_min(Avg,n=45) %>% pull(gene)
top_genes <-  c(top_down,top_up) 

top_genes_2 <- c("APLN","APLNR","CCND1","CD34","COL15A1","COL4A1","COL4A2","CXCL12","EBF1","ENPP2","GJA1","HIF1A","IGFBP7","INHBB","KDR","LAMA4","LAMB1","LAMC1","NREP","PDGFB","PDLIM1","RGCC","SMAD1","SOX17","SOX4","SPARC","SPRY1","THBS1","TNFRSF10B","VWA1",
                 "ACVRL1","BMPR2","CEBPD","CRYAB","DUSP1","ELN","FOXF1","FGFR3","GATA2","GATA6","HES1","KLF2","KLF4","KLF9","PTGS1","SMAD6","SMAD7","SORBS1","SRGN","TEK","TGFB2","THBD","TIMP3","TM6SF1","VEGFC","VWF",
                 "STXBP6","ABI3BP","TMEM100","PTGIS","HDAC9","ARL4D","INMT","ADGRG6","LTBP1","SAMD5","MGLL","MEOX1","FAM13C","PRDM1","FLT4","RASSF2","ANXA3","ANKRD44","EMP3")



sCap_IPF <- DEG.psdblk.sCap  %>% mutate(sCap_IPF=log2FoldChange)%>% dplyr::select(gene,sCap_IPF)
sCap_Bleo <- read.xlsx("sCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
sCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(sCap_Bleo$gene),sCap_Bleo=sCap_Bleo$Bleo) 
sCap <- merge(sCap_Bleo,sCap_IPF)%>% filter(gene%in%top_genes)


DEG.psdblk.SVEC <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/Supplementary_table_S4.xlsx",sheet = "SVEC")
SVEC_IPF_1 <- DEG.psdblk.SVEC %>% filter(padj < 0.05) %>% pull(gene) 
SVEC_Bleo_1 <- read.xlsx("SVEC_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
SVEC_Bleo_1 <- SVEC_Bleo_1[complete.cases(SVEC_Bleo_1)]

intersect_SVEC <- intersect(SVEC_Bleo_1,SVEC_IPF_1)

SVEC_IPF <- DEG.psdblk.SVEC %>% filter(gene %in% intersect_SVEC) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
SVEC_Bleo <- read.xlsx("SVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
SVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(SVEC_Bleo$gene),Bleo=SVEC_Bleo$Bleo) %>% filter(gene %in% intersect_SVEC)
SVEC <- merge(SVEC_Bleo,SVEC_IPF)
SVEC$comparison <- sign(SVEC$Bleo) == sign(SVEC$IPF)
SVEC <- SVEC %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5E
nrow(SVEC)
length(SVEC_IPF_1)-nrow(SVEC)
length(SVEC_Bleo_1)-nrow(SVEC)
total <- sum(nrow(SVEC),c(length(SVEC_IPF_1)-nrow(SVEC)),c(length(SVEC_Bleo_1)-nrow(SVEC)))

nrow(SVEC)/total*100
c(length(SVEC_IPF_1)-nrow(SVEC))/total*100
c(length(SVEC_Bleo_1)-nrow(SVEC))/total*100

SVEC_IPF <- DEG.psdblk.SVEC  %>% mutate(SVEC_IPF=log2FoldChange)%>% dplyr::select(gene,SVEC_IPF)
SVEC_Bleo <- read.xlsx("SVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
SVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(SVEC_Bleo$gene),SVEC_Bleo=SVEC_Bleo$Bleo) 
SVEC <- merge(SVEC_Bleo,SVEC_IPF)%>% filter(gene%in%top_genes)


DEG.psdblk.PVEC <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/Supplementary_table_S4.xlsx",sheet = "PVEC")
PVEC_IPF_1 <- DEG.psdblk.PVEC %>% filter(padj < 0.05) %>% pull(gene) 
PVEC_Bleo_1 <- read.xlsx("PVEC_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
PVEC_Bleo_1 <- PVEC_Bleo_1[complete.cases(PVEC_Bleo_1)]

intersect_PVEC <- intersect(PVEC_Bleo_1,PVEC_IPF_1)

PVEC_IPF <- DEG.psdblk.PVEC %>% filter(gene %in% intersect_PVEC) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
PVEC_Bleo <- read.xlsx("PVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
PVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(PVEC_Bleo$gene),Bleo=PVEC_Bleo$Bleo) %>% filter(gene %in% intersect_PVEC)
PVEC <- merge(PVEC_Bleo,PVEC_IPF)
PVEC$comparison <- sign(PVEC$Bleo) == sign(PVEC$IPF)
PVEC <- PVEC %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5E
nrow(PVEC)
length(PVEC_IPF_1)-nrow(PVEC)
length(PVEC_Bleo_1)-nrow(PVEC)
total <- sum(nrow(PVEC),c(length(PVEC_IPF_1)-nrow(PVEC)),c(length(PVEC_Bleo_1)-nrow(PVEC)))

nrow(PVEC)/total*100
c(length(PVEC_IPF_1)-nrow(PVEC))/total*100
c(length(PVEC_Bleo_1)-nrow(PVEC))/total*100

PVEC_IPF <- DEG.psdblk.PVEC  %>% mutate(PVEC_IPF=log2FoldChange)%>% dplyr::select(gene,PVEC_IPF)
PVEC_Bleo <- read.xlsx("PVEC_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
PVEC_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(PVEC_Bleo$gene),PVEC_Bleo=PVEC_Bleo$Bleo) 
PVEC <- merge(PVEC_Bleo,PVEC_IPF)%>% filter(gene%in%top_genes)


DEG.psdblk.gCap <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/Supplementary_table_S4.xlsx",sheet = "gCap")
gCap_IPF_1 <- DEG.psdblk.gCap %>% filter(padj < 0.05) %>% pull(gene) 
gCap_Bleo_1 <- read.xlsx("Lrg1pos_gCap_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
gCap_Bleo_1 <- gCap_Bleo_1[complete.cases(gCap_Bleo_1)]

intersect_gCap <- intersect(gCap_Bleo_1,gCap_IPF_1)

gCap_IPF <- DEG.psdblk.gCap %>% filter(gene %in% intersect_gCap) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
gCap_Bleo <- read.xlsx("Lrg1pos_gCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
gCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(gCap_Bleo$gene),Bleo=gCap_Bleo$Bleo) %>% filter(gene %in% intersect_gCap)
gCap <- merge(gCap_Bleo,gCap_IPF)
gCap$comparison <- sign(gCap$Bleo) == sign(gCap$IPF)
gCap <- gCap %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5E
nrow(gCap)
length(gCap_IPF_1)-nrow(gCap)
length(gCap_Bleo_1)-nrow(gCap)
total <- sum(nrow(gCap),c(length(gCap_IPF_1)-nrow(gCap)),c(length(gCap_Bleo_1)-nrow(gCap)))

nrow(gCap)/total*100
c(length(gCap_IPF_1)-nrow(gCap))/total*100
c(length(gCap_Bleo_1)-nrow(gCap))/total*100

gCap_IPF <- DEG.psdblk.gCap  %>% mutate(gCap_IPF=log2FoldChange)%>% dplyr::select(gene,gCap_IPF)
gCap_Bleo <- read.xlsx("Lrg1pos_gCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
gCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(gCap_Bleo$gene),gCap_Bleo=gCap_Bleo$Bleo) 
gCap <- merge(gCap_Bleo,gCap_IPF)%>% filter(gene%in%top_genes)

DEG.psdblk.aCap <- read.xlsx("/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/Supplementary_table_S4.xlsx",sheet = "aCap")
aCap_IPF_1 <- DEG.psdblk.aCap %>% filter(padj < 0.05) %>% pull(gene) 
aCap_Bleo_1 <- read.xlsx("Lrg1pos_aCap_pseudobulk.xlsx") %>% filter(padj < 0.05) %>% pull(gene) %>% convert_mouse_to_human_symbols() 
aCap_Bleo_1 <- aCap_Bleo_1[complete.cases(aCap_Bleo_1)]

intersect_aCap <- intersect(aCap_Bleo_1,aCap_IPF_1)

aCap_IPF <- DEG.psdblk.aCap %>% filter(gene %in% intersect_aCap) %>% mutate(IPF=log2FoldChange)%>% dplyr::select(gene,IPF)
aCap_Bleo <- read.xlsx("Lrg1pos_aCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
aCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(aCap_Bleo$gene),Bleo=aCap_Bleo$Bleo) %>% filter(gene %in% intersect_aCap)
aCap <- merge(aCap_Bleo,aCap_IPF)
aCap$comparison <- sign(aCap$Bleo) == sign(aCap$IPF)
aCap <- aCap %>% filter(comparison == TRUE)

#Numbers for Venn Diagramm of figure 5E
nrow(aCap)
length(aCap_IPF_1)-nrow(aCap)
length(aCap_Bleo_1)-nrow(aCap)
total <- sum(nrow(aCap),c(length(aCap_IPF_1)-nrow(aCap)),c(length(aCap_Bleo_1)-nrow(aCap)))

nrow(aCap)/total*100
c(length(aCap_IPF_1)-nrow(aCap))/total*100
c(length(aCap_Bleo_1)-nrow(aCap))/total*100

aCap_IPF <- DEG.psdblk.aCap  %>% mutate(aCap_IPF=log2FoldChange)%>% dplyr::select(gene,aCap_IPF)
aCap_Bleo <- read.xlsx("Lrg1pos_aCap_pseudobulk.xlsx") %>% mutate(Bleo=log2FoldChange)%>% dplyr::select(gene,Bleo) 
aCap_Bleo <- data.frame(gene=convert_mouse_to_human_symbols(aCap_Bleo$gene),aCap_Bleo=aCap_Bleo$Bleo) 
aCap <- merge(aCap_Bleo,aCap_IPF)%>% filter(gene%in%top_genes)



data <- Reduce(f = inner_join,list(sCap,gCap,aCap,SVEC,PVEC)) %>% column_to_rownames(.,"gene")
data_full <- Reduce(f = full_join,list(sCap,gCap,aCap,SVEC,PVEC)) %>% column_to_rownames(.,"gene")
data <-data_full[c(rownames(data),"COL15A1"),]
cellstate_order <- c("sCap","gCap","aCap","SVEC","PVEC")
sample.order <- c("IPF","Bleo")
annotation <- data.frame(group=colnames(data),sample=gsub(".*_","",colnames(data)),celltype=gsub("_.*","",colnames(data)))
annotation$sample <- factor(annotation$sample,levels = sample.order)
annotation$celltype <- factor(annotation$celltype,levels = cellstate_order)
annotation <-arrange(annotation,celltype,sample)
data <- data[,match(annotation$group,colnames(data))]
annotation <- annotation %>% column_to_rownames(.,"group")
data[data < -4] <- -4
data[data > 4] <- 4
data[is.na(data)] <- 0
subcolors <- c("sCap"="#A6D854","gCap"="#911414","aCap"="#5e1e33","SVEC"="#ccc357","PVEC"="#f2ae5a")
paletteLength = 100
myColor = colorRampPalette(rev(brewer.pal(n = 9, name ="RdBu")))(paletteLength)
Breaks = c(seq(min(data), 0, length.out=ceiling(paletteLength/2) + 1),
           seq(max(data)/paletteLength, max(data), length.out=floor(paletteLength/2)))

# ---- Fig3H ---- 
pheatmap(t(data),annotation_row  = annotation,show_rownames = F,cluster_cols = T,cluster_rows = F,color = myColor,angle_col =90,fontsize_col = 8,treeheight_col=0,
         annotation_colors = list(celltype= subcolors,sample=c("IPF"="#5e5e5e","Bleo"="#b6b6b8")))

#For Table 2
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- sCap$gene
df <- listAttributes(mart)
G_list <- getBM(filters= "hgnc_symbol", attributes= c("description","external_gene_name"),values=genes,mart= mart) %>% mutate(gene=external_gene_name)

sCap <- merge(sCap,G_list,by.x="gene")
sCap$description <- gsub("\\[Source.*","",sCap$description)
sCap$expression_in_fibrosis <- ifelse(sCap$sCap_Bleo>0,"upregulated","downregulated")
sCap <- sCap %>%  dplyr::select(gene,description,expression_in_fibrosis)
write.xlsx(unique(sCap),file = "/home/truchi/LUNG_FIBROSIS/BLEO_Integrated_Analysis_FINAL/Reviews_NatComm/table2.xlsx")


