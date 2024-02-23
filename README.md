# Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury


This is the code from : 

Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury

Marin Truchi*, Grégoire Savary*, Marine Gautier-Isola, Hugo Cadis, Alberto Baeri, Arun Lingampally, Célia Scribe, Virginie Magnone, Cédric Girard-Riboulleau,  Marie-Jeanne Arguel, Clémentine de Schutter, Julien Fassy, Nihad Boukrout, Romain Larrue, Nathalie Martin, Roger Rezzonico, Olivier Pluquet, Michael Perrais, Véronique Hofman, Charles-Hugo Marquette, Paul Hofman, Andreas Günther, Nicolas Ricard, Pascal Barbry, Sylvie Leroy, Kevin Lebrigand, Saverio Bellusci, Christelle Cauffiez, Georges Vassaux, Nicolas Pottier* and Bernard Mari*#

[bioRxiv 2023.07.11.548494; doi: https://doi.org/10.1101/2023.07.11.548494](https://www.biorxiv.org/content/10.1101/2023.07.11.548522v1)


The goal of this repository is to gather materials used for i) the preprocessing of scRNA-seq and spatial transcriptomics data, ii) the integration of scRNA-seq samples, iii) the figures of the manuscript. 

---

<figure>
  <img src="https://github.com/marintruchi/Aging_affects_reprogramming_of_PCEC/blob/main/graphical_abstract.PNG" alt="SSAE_overview"/>
</figcaption>
</figure>

---

## **Repository Contents**
|Folder | Description |
|:----------|:----------|
|`.rds files`|Contains the .rds files used to produce figures, corresponding to seurat objects of this study (scRNA-seq data + spatial data) as well as public datasets (Strunz et al. and Habermann et al.)|
|`scripts`|Contains the scripts for the preprocessing of scRNA-seq data, their integration and the production of the manuscript's figures|
|`additional files`|Contains additional files like Ingenuity Pathways Analysis outputs and supplementary tables|

 ### **"scripts" repository Contents**   
|Script| Description |
|:----------|:----------|
|`Preprocess_CROPseq_lib.R`|R script to load, manipulate and prepare the count matrices for the SSAE |
|`Run_SSAE_script.py`|Main python script to run the SSAE|
|`Produce_figures.R`|R script to produce the figures of the manuscript|



