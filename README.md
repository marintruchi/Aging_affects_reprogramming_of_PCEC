# Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury


This is the code from : 

Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury

Marin Truchi*, Grégoire Savary*, Marine Gautier-Isola, Hugo Cadis, Alberto Baeri, Arun Lingampally, Célia Scribe, Virginie Magnone, Cédric Girard-Riboulleau,  Marie-Jeanne Arguel, Clémentine de Schutter, Julien Fassy, Nihad Boukrout, Romain Larrue, Nathalie Martin, Roger Rezzonico, Olivier Pluquet, Michael Perrais, Véronique Hofman, Charles-Hugo Marquette, Paul Hofman, Andreas Günther, Nicolas Ricard, Pascal Barbry, Sylvie Leroy, Kevin Lebrigand, Saverio Bellusci, Christelle Cauffiez, Georges Vassaux, Nicolas Pottier* and Bernard Mari*#

[bioRxiv 2023.07.11.548494; doi: https://doi.org/10.1101/2023.07.11.548494](https://www.biorxiv.org/content/10.1101/2023.07.11.548522v1)


The goal of this repository is to gather materials used for i) the preprocessing of scRNA-seq and spatial transcriptomics data, ii) the integration of scRNA-seq samples, iii) the figures of the manuscript. 



## **Repository Contents**
|Folder | Description |
|:----------|:----------|
|`.rds files`|Contains the .rds files used to produce figures, corresponding to seurat objects of this study (scRNA-seq data + spatial data) as well as public datasets (Strunz et al. and Habermann et al.)|
|`scripts`|Contains the scripts for the preprocessing of scRNA-seq data, their integration and the production of the manuscript's figures|
|`additional files`|Contains additional files like Ingenuity Pathways Analysis outputs and supplementary tables|

 ### **"scripts" repository Contents**   
|Script| Description |
|:----------|:----------|
|`prerequisites.R`|R script to load required packages and functions |
|`scRNAseq_preprocessing.R`|R script to preprocess the scRNA-seq data of each time point |
|`scRNAseq_RPCA_integration.R`|R script to integrate all samples in a single seurat object |
|`Fig1_SupplFigS1.R`|R script to produce the plots of Figure 1 and Suppl. Figure S1 |
|`Fig2_SupplFigS2&S3.R`|R script to produce the plots of Figure 2 and Suppl. Figure S2&3 |
|`Fig3_SupplFigS4&S5.R`|R script to produce the plots of Figure 3 and Suppl. Figure S4&5 |
|`Fig4&5_SupplFigS6.R`|R script to produce the plots of Figure 4&5 and Suppl. Figure S6 |
|`Fig6_SupplFigS7.R`|R script to produce the plots of Figure 6 and Suppl. Figure S7 |
|`scVelo_gCap_Yd14.ipynb`|Jupyter python notebook to perform the scVelo analysis on the gCap subset from young D14 sample |
|`scVelo_gCap_Od28.ipynb`|Jupyter python notebook to perform the scVelo analysis on the gCap subset from old D28 sample |












