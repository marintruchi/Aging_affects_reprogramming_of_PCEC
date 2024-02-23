# Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury


This is the code from : 

Aging affects reprogramming of murine pulmonary capillary endothelial cells after lung injury

Marin Truchi*, Grégoire Savary*, Marine Gautier-Isola, Hugo Cadis, Alberto Baeri, Arun Lingampally, Célia Scribe, Virginie Magnone, Cédric Girard-Riboulleau,  Marie-Jeanne Arguel, Clémentine de Schutter, Julien Fassy, Nihad Boukrout, Romain Larrue, Nathalie Martin, Roger Rezzonico, Olivier Pluquet, Michael Perrais, Véronique Hofman, Charles-Hugo Marquette, Paul Hofman, Andreas Günther, Nicolas Ricard, Pascal Barbry, Sylvie Leroy, Kevin Lebrigand, Saverio Bellusci, Christelle Cauffiez, Georges Vassaux, Nicolas Pottier* and Bernard Mari*#

[bioRxiv 2023.07.11.548494; doi: https://doi.org/10.1101/2023.07.11.548494](https://www.biorxiv.org/content/10.1101/2023.07.11.548522v1)


The goal of this repository is to gather materials used for i) the preprocessing of scRNA-seq and spatial transcriptomics data, ii) the integration of scRNA-seq samples, iii) the figures of the manuscript. 

---

<figure>
  <img src="https://github.com/marintruchi/Aging_affects_reprogramming_of_PCEC/blob/main/graphical_abstract.PNG" alt="SSAE_overview"/>
cells
</figcaption>
</figure>

---

## **Repository Contents**
|Folder | Description |
|:----------|:----------|
|`data`|Contains the input matrices of the SSAE and the link to access the raw and preprocessed sequencing files of the CROP-seq library|
|`scripts`|Contains the scripts for the preprocessing of CROP-seq data, the preparation of SSAE input matrices, the execution of SSAE and the production of the manuscript's figures|
|`functions`|Contains dedicated functions for the execution of SSAE|

 ### **"scripts" repository Contents**   
|Script| Description |
|:----------|:----------|
|`Preprocess_CROPseq_lib.R`|R script to load, manipulate and prepare the count matrices for the SSAE |
|`Run_SSAE_script.py`|Main python script to run the SSAE|
|`Produce_figures.R`|R script to produce the figures of the manuscript|



## **How to run the SSAE** 

We recommend installing and running the SSAE in an anaconda environment.

### Requirements
- python >= 3.8 .
- [Pytorch](https://pytorch.org/get-started/locally/).
- The following packages : [numpy](https://numpy.org/install/), [matplotlib](https://matplotlib.org/stable/users/installing/index.html), [scikit-learn](https://scikit-learn.org/stable/install.html), [pandas](https://pandas.pydata.org/getting_started.html), [shap](https://pypi.org/project/shap/), [captum](https://captum.ai/#quickstart). 



### How to use it

For now, the SSAE can only be run as a raw script (`Run_SSAE_script.py`).

In the working directory of the `Run_SSAE_script.py` script, you must have the `functions` folder, and a `data` folder containing the input count matrix (.csv file) you want to work on.

For each run, you have to edit the `Run_SSAE_script.py` script to specify the name of the input matrix you want to work on in the `data` folder. The string before the first "." of the file name will be used for the names of the SSAE outputs. 

You can run the script in your favorite IDE, or with the command line `python Run_SSAE_script.py`.

ETA is the only parameter which may be ajusted, according to the number of cells/samples. In the manuscript, ETA has been set at 25 for all tests.


### SSAE outputs

For each run, the SSAE produce 2 major outputs :
1) A folder of files in a  `results_stat` folder, containing the following .csv files :

|file| Description |
|:----------|:----------|
|`accuracy_test.csv`|The global accuracy and the accuracy for each class (for each Fold, with the mean and the standard deviation)|
|`metrics_test.csv`|Other metrics for each Fold, with the mean and the standard deviation) |
|`Labelspred_softmax.csv`|The classification scores for each cell ("Proba class 0" = control score, "Proba class 1" = perturbation score)|
|`SSAE_selection.csv`|The results of the SSAE classification (**Ntarget** = nb of gRNA-targeted cells, **Np** = nb of "perturbed" cells, **NNp** = nb of "non-perturbed" cells, **NNegative** = nb of control cells, **NNegative_R** = nb of randomly selected control cells)|
|`proj_l11ball_topGenes_Captum_dl_300.csv`|A list of the most discriminant features between the two compared classes, ranked by their associated weights (for each Fold, with the mean and the standard deviation)|
---
---

2) An output count matrix (`[filename]_after_selection.csv` file within the `data` folder, containing only the raw expression profiles of gRNA-targeted cells classified as "perturbed" and an equivalent number of control cells, randomly selected. If the number of "perturbed" cells is sufficient (~100), this matrix of selected cells can be used as input for a second-turn SSAE run to extract a final list of selected features corresponding to the refined perturbation signature.

