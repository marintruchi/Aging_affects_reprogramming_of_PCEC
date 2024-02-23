# Paper : Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury
# Code : R script to preprocess the scRNA-seq data
# Author : Marin Truchi

#Load required packages and functions
source("prerequisites.R")

#-------------------------------------------------------------------------------------------------#
#--------------- Load RNA & HTO count tables and integrate them for each time point --------------#
#-------------------------------------------------------------------------------------------------#

## ---------------------------- YoungD14 ---------------------------- ##
#LOAD HTO COUNTS
matrix_dir = "~/GEO_data_submission_Truchi_et_al/YD14/"
barcode.path <- paste0(matrix_dir, "YoungD14_HTOs_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "YoungD14_HTOs_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "YoungD14_HTOs_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
colnames(mat) <- paste0("YoungD14_",colnames(mat))
rownames(mat) = feature.names$V1
hto = mat[1:6,]

#LOAD RNA COUNTS
barcode.path <- paste0(matrix_dir, "YoungD14_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "YoungD14_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "YoungD14_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- gsub("-.*","",barcode.names$V1)
rownames(mat) = feature.names$V2
mat <-  as(mat, "dgCMatrix")
colnames(mat) <- paste0("YoungD14_",colnames(mat))

hashing <- CreateSeuratObject(mat,min.cells = 1,project = "YoungD14")

ncolhto=colnames(hto)
hashing<-subset(hashing, cells=ncolhto) #if hto cell nb != cell nb  

mito.genes <- grep(pattern = "^mt-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts

DefaultAssay(hashing) <- "RNA"
genes <- rownames(hashing)[-which(rownames(hashing)%in%c(mito.genes,ribo.genes))]
hashing$nCount_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nCount_RNA
hashing$nFeature_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nFeature_RNA

hashing[["HTO"]] <- CreateAssayObject(counts = hto)
hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="YoungD14")

VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)
dim(hashing@assays$RNA)
hashing <- subset(hashing, subset = percent.mito < 0.25 & nCount_RNA_noMR < 50000 & nFeature_RNA_noMR > 200)
dim(hashing2@assays$RNA)
VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)

hashing <- SetIdent(hashing, value = "HTO_maxID")
table(hashing@active.ident)
HTOHeatmap(hashing)

hashing[['condition']] <- "PBS"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0304-AAAGCATTCTTCACG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0305-CTTTGTCTTTGTGAG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0306-TATGCTGCCACGGTA"),] <- "Bleo"
hashing$sample <- paste(substr(hashing$hash.ID,1,5),hashing$condition, sep = "_")

print(table(hashing@active.ident, hashing@meta.data$condition))

hashing<-subset(hashing, subset=HTO_classification.global == c("Singlet"))
hashing$barcode <- paste0("YoungD14_",colnames(hashing))
hashing$orig.ident <- "Young_D14"
list.barcodes <- list("YoungD14"=hashing)


## ---------------------------- YoungD28 ---------------------------- ##
#LOAD HTO COUNTS
matrix_dir = "~/GEO_data_submission_Truchi_et_al/YD28/"
barcode.path <- paste0(matrix_dir, "YoungD28_HTOs_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "YoungD28_HTOs_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "YoungD28_HTOs_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
colnames(mat) <- paste0("YoungD28_",colnames(mat))
rownames(mat) = feature.names$V1
hto = mat[1:6,]

#LOAD RNA COUNTS
barcode.path <- paste0(matrix_dir, "YoungD28_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "YoungD28_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "YoungD28_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- gsub("-.*","",barcode.names$V1)
colnames(mat) <- paste0("YoungD28_",colnames(mat))
rownames(mat) = feature.names$V2
mat <-  as(mat, "dgCMatrix")

hashing <- CreateSeuratObject(mat,min.cells = 1,project = "YoungD28")

ncolhto=colnames(hto)
hashing<-subset(hashing, cells=ncolhto) #if hto cell nb != cell nb  

mito.genes <- grep(pattern = "^mt-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts

DefaultAssay(hashing) <- "RNA"
genes <- rownames(hashing)[-which(rownames(hashing)%in%c(mito.genes,ribo.genes))]
hashing$nCount_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nCount_RNA
hashing$nFeature_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nFeature_RNA

hashing[["HTO"]] <- CreateAssayObject(counts = hto)
hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="YoungD28")

VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)
dim(hashing@assays$RNA)
hashing <- subset(hashing, subset = percent.mito < 0.25 & nCount_RNA_noMR < 30000 & nFeature_RNA_noMR > 200)
dim(hashing@assays$RNA)
VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)

hashing <- SetIdent(hashing, value = "HTO_maxID")
table(hashing@active.ident)
HTOHeatmap(hashing)

hashing[['condition']] <- "PBS"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0304-AAAGCATTCTTCACG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0305-CTTTGTCTTTGTGAG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0306-TATGCTGCCACGGTA"),] <- "Bleo"
hashing$sample <- paste(substr(hashing$hash.ID,1,5),hashing$condition, sep = "_")

print(table(hashing@active.ident, hashing@meta.data$condition))

hashing<-subset(hashing, subset=HTO_classification.global == c("Singlet"))
hashing$barcode <- paste0("YoungD28_",colnames(hashing))
hashing$orig.ident <- "Young_D28"
list.barcodes <- list.append(list.barcodes,"YoungD28"=hashing) 


## ---------------------------- YoungD60 ---------------------------- ##
#LOAD HTO COUNTS
matrix_dir = "~/GEO_data_submission_Truchi_et_al/YD60/"
barcode.path <- paste0(matrix_dir, "YoungD60_HTOs_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "YoungD60_HTOs_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "YoungD60_HTOs_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
colnames(mat) <- paste0("YoungD60_",colnames(mat))
rownames(mat) = feature.names$V1
hto = mat[1:6,]

#LOAD RNA COUNTS
barcode.path <- paste0(matrix_dir, "YoungD60_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "YoungD60_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "YoungD60_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- gsub("-.*","",barcode.names$V1)
colnames(mat) <- paste0("YoungD60_",colnames(mat))
rownames(mat) = feature.names$V2
mat <-  as(mat, "dgCMatrix")

hashing <- CreateSeuratObject(mat,min.cells = 1,project = "YoungD60")

ncolhto=colnames(hto)
hashing<-subset(hashing, cells=ncolhto) #if hto cell nb != cell nb  

mito.genes <- grep(pattern = "^mt-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts

DefaultAssay(hashing) <- "RNA"
genes <- rownames(hashing)[-which(rownames(hashing)%in%c(mito.genes,ribo.genes))]
hashing$nCount_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nCount_RNA
hashing$nFeature_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nFeature_RNA

hashing[["HTO"]] <- CreateAssayObject(counts = hto)
hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="YoungD60")

VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)
dim(hashing@assays$RNA)
hashing <- subset(hashing, subset = percent.mito < 0.25 & nCount_RNA_noMR < 40000 & nFeature_RNA_noMR > 200)
dim(hashing@assays$RNA)
VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)

hashing <- SetIdent(hashing, value = "HTO_maxID")
table(hashing@active.ident)
HTOHeatmap(hashing)

hashing[['condition']] <- "PBS"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0304-AAAGCATTCTTCACG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0305-CTTTGTCTTTGTGAG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0306-TATGCTGCCACGGTA"),] <- "Bleo"
hashing$sample <- paste(substr(hashing$hash.ID,1,5),hashing$condition, sep = "_")

print(table(hashing@active.ident, hashing@meta.data$condition))

hashing<-subset(hashing, subset=HTO_classification.global == c("Singlet"))
hashing$barcode <- paste0("YoungD60_",colnames(hashing))
hashing$orig.ident <- "Young_D60"
list.barcodes <- list.append(list.barcodes,"YoungD60"=hashing) 


## ---------------------------- OldD14 ---------------------------- ##
#LOAD HTO COUNTS
matrix_dir = "~/GEO_data_submission_Truchi_et_al/OD14/"
barcode.path <- paste0(matrix_dir, "OldD14_HTOs_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "OldD14_HTOs_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "OldD14_HTOs_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
colnames(mat) <- paste0("OldD14_",colnames(mat))
rownames(mat) = feature.names$V1
hto = mat[1:6,]

#LOAD RNA COUNTS
barcode.path <- paste0(matrix_dir, "OldD14_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "OldD14_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "OldD14_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- gsub("-.*","",barcode.names$V1)
colnames(mat) <- paste0("OldD14_",colnames(mat))
rownames(mat) = feature.names$V2
mat <-  as(mat, "dgCMatrix")

hashing <- CreateSeuratObject(mat,min.cells = 1,project = "OldD14")

ncolhto=colnames(hto)
hashing<-subset(hashing, cells=ncolhto) #if hto cell nb != cell nb  

mito.genes <- grep(pattern = "^mt-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts

DefaultAssay(hashing) <- "RNA"
genes <- rownames(hashing)[-which(rownames(hashing)%in%c(mito.genes,ribo.genes))]
hashing$nCount_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nCount_RNA
hashing$nFeature_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nFeature_RNA


hashing[["HTO"]] <- CreateAssayObject(counts = hto)
hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="OldD14")

VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)
dim(hashing@assays$RNA)
hashing <- subset(hashing, subset = percent.mito < 0.25 & nCount_RNA_noMR < 50000 & nFeature_RNA_noMR > 200)
dim(hashing@assays$RNA)
VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)

hashing <- SetIdent(hashing, value = "HTO_maxID")
table(hashing@active.ident)
HTOHeatmap(hashing)

hashing[['condition']] <- "PBS"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0304-AAAGCATTCTTCACG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0305-CTTTGTCTTTGTGAG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0306-TATGCTGCCACGGTA"),] <- "Bleo"
hashing$sample <- paste(substr(hashing$hash.ID,1,5),hashing$condition, sep = "_")

print(table(hashing@active.ident, hashing@meta.data$condition))

hashing<-subset(hashing, subset=HTO_classification.global == c("Singlet"))
hashing$barcode <- paste0("OldD14_",colnames(hashing))
hashing$orig.ident <- "Old_D14"
list.barcodes <- list.append(list.barcodes,"OldD14"=hashing)

## ---------------------------- OldD28 ---------------------------- ##
#LOAD HTO COUNTS
matrix_dir = "~/GEO_data_submission_Truchi_et_al/OD28/"
barcode.path <- paste0(matrix_dir, "OldD28_HTOs_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "OldD28_HTOs_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "OldD28_HTOs_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
colnames(mat) <- paste0("OldD28_",colnames(mat))
rownames(mat) = feature.names$V1
hto = mat[1:6,]

#LOAD RNA COUNTS
barcode.path <- paste0(matrix_dir, "OldD28_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "OldD28_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "OldD28_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- gsub("-.*","",barcode.names$V1)
colnames(mat) <- paste0("OldD28_",colnames(mat))
rownames(mat) = feature.names$V2
mat <-  as(mat, "dgCMatrix")

hashing <- CreateSeuratObject(mat,min.cells = 1,project = "OldD28")

ncolhto=colnames(hto)
hashing<-subset(hashing, cells=ncolhto) #if hto cell nb != cell nb  

mito.genes <- grep(pattern = "^mt-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts

DefaultAssay(hashing) <- "RNA"
genes <- rownames(hashing)[-which(rownames(hashing)%in%c(mito.genes,ribo.genes))]
hashing$nCount_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nCount_RNA
hashing$nFeature_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nFeature_RNA

hashing[["HTO"]] <- CreateAssayObject(counts = hto)
hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="OldD28")

VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)
dim(hashing@assays$RNA)
hashing <- subset(hashing, subset = percent.mito < 0.25 & nCount_RNA_noMR < 30000 & nFeature_RNA_noMR > 200)
dim(hashing@assays$RNA)
VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)

hashing <- SetIdent(hashing, value = "HTO_maxID")
table(hashing@active.ident)
HTOHeatmap(hashing)

hashing[['condition']] <- "PBS"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0304-AAAGCATTCTTCACG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0305-CTTTGTCTTTGTGAG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0306-TATGCTGCCACGGTA"),] <- "Bleo"
hashing$sample <- paste(substr(hashing$hash.ID,1,5),hashing$condition, sep = "_")

print(table(hashing@active.ident, hashing@meta.data$condition))

hashing<-subset(hashing, subset=HTO_classification.global == c("Singlet"))
hashing$barcode <- paste0("OldD28_",colnames(hashing))
hashing$orig.ident <- "Old_D28"
list.barcodes <- list.append(list.barcodes,"OldD28"=hashing) 


## ---------------------------- OldD60 ---------------------------- ##
#LOAD HTO COUNTS
matrix_dir = "~/GEO_data_submission_Truchi_et_al/OD60/"
barcode.path <- paste0(matrix_dir, "OldD60_HTOs_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "OldD60_HTOs_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "OldD60_HTOs_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) = barcode.names$V1
colnames(mat) <- paste0("OldD60_",colnames(mat))
rownames(mat) = feature.names$V1
hto = mat[1:6,]

#LOAD RNA COUNTS
barcode.path <- paste0(matrix_dir, "OldD60_barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "OldD60_features.tsv.gz")
matrix.path <- paste0(matrix_dir, "OldD60_matrix.mtx.gz")
mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
colnames(mat) <- gsub("-.*","",barcode.names$V1)
colnames(mat) <- paste0("OldD60_",colnames(mat))
rownames(mat) = feature.names$V2
mat <-  as(mat, "dgCMatrix")

hashing <- CreateSeuratObject(mat,min.cells = 1,project = "OldD60")

ncolhto=colnames(hto)
hashing<-subset(hashing, cells=ncolhto) #if hto cell nb != cell nb  

mito.genes <- grep(pattern = "^mt-", x = rownames(hashing@assays$RNA), value = TRUE)
ribo.genes <- grep(pattern = "^Rp[sl]", x = rownames(hashing@assays$RNA), value = TRUE)
dropouts <- Matrix::colSums(hashing@assays$RNA@data == 0)/nrow(hashing@assays$RNA)
percent.mito <- Matrix::colSums(hashing@assays$RNA[mito.genes, ])/Matrix::colSums(hashing@assays$RNA)
percent.ribo <- Matrix::colSums(hashing@assays$RNA[ribo.genes, ])/Matrix::colSums(hashing@assays$RNA)
hashing[['percent.mito']] <- percent.mito
hashing[['percent.ribo']] <- percent.ribo
hashing[['dropouts']] <- dropouts

DefaultAssay(hashing) <- "RNA"
genes <- rownames(hashing)[-which(rownames(hashing)%in%c(mito.genes,ribo.genes))]
hashing$nCount_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nCount_RNA
hashing$nFeature_RNA_noMR <- hashing %>% subset(features = genes) %>% .$nFeature_RNA

hashing[["HTO"]] <- CreateAssayObject(counts = hto)
hashing <- NormalizeData(object = hashing, assay = "HTO", normalization.method = "CLR")
hashing <- HTODemux(hashing)
hashing <- SetIdent(hashing, value="OldD60")

VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)
dim(hashing@assays$RNA)
hashing <- subset(hashing, subset = percent.mito < 0.25 & nCount_RNA_noMR < 50000 & nFeature_RNA_noMR > 200)
dim(hashing@assays$RNA)
VlnPlot(hashing, features = c("nFeature_RNA_noMR", "nCount_RNA_noMR","dropouts","percent.ribo","percent.mito"), ncol=5, cols = "lightsteelblue3",pt.size = 0.25)

hashing <- SetIdent(hashing, value = "HTO_maxID")
table(hashing@active.ident)
HTOHeatmap(hashing)

hashing[['condition']] <- "PBS"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0304-AAAGCATTCTTCACG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0305-CTTTGTCTTTGTGAG"),] <- "Bleo"
hashing[['condition']][which(hashing@meta.data$HTO_maxID == "A0306-TATGCTGCCACGGTA"),] <- "Bleo"
hashing$sample <- paste(substr(hashing$hash.ID,1,5),hashing$condition, sep = "_")

print(table(hashing@active.ident, hashing@meta.data$condition))

hashing<-subset(hashing, subset=HTO_classification.global == c("Singlet"))
hashing$barcode <- paste0("OldD60_",colnames(hashing))
hashing$orig.ident <- "Old_D60"
list.barcodes <- list.append(list.barcodes,"OldD60"=hashing)


for (i in 1:length(list.barcodes)) {
  list.barcodes[[i]]$mice <- paste0(list.barcodes[[i]]$orig.ident, ".",list.barcodes[[i]]$sample)
  list.barcodes[[i]]@meta.data <-list.barcodes[[i]]@meta.data %>%  dplyr::select(barcode,orig.ident,HTO_classification.global,condition,mice,
                                                                             nCount_RNA,nFeature_RNA,percent.mito,percent.ribo,dropouts,nCount_HTO,nFeature_HTO,
                                                                             HTO_maxID, HTO_secondID,HTO_margin) 
  
}

#Aggregate all samples in a single object
aggr <- Reduce(merge,list.barcodes)
condition_order <- c( 'Bleo','PBS')
aggr$condition <- factor(aggr$condition, levels = condition_order)
aggr$sample <- paste0(aggr$orig.ident,".",aggr$condition)
aggr$orig.ident <- gsub("_"," ",aggr$orig.ident)
sampletype_order <- c('Young D14', 'Old D14','Young D28', 'Old D28','Young D60', 'Old D60')
aggr$orig.ident <- factor(aggr$orig.ident, levels = sampletype_order)
sample.colors=c("Old D14"="#B2ABD2","Old D60"="#542788","Old D28"="#8073AC","Young D14"="#FDB863","Young D28"="#E08214","Young D60"="#B35806")
rm(list.barcodes)

#save unnormalized and non integrated data
saveRDS(aggr,"~/before_integration.rds")

