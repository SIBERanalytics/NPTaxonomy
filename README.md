# Taxonomical classification
![alt text](https://github.com/SIBERanalytics/NPTaxonomy/blob/main/overview_figure.png?raw=true)

## Overview

This repository contains source data and python scripts for natural products (NPs) taxonomical classification and structural insights. From the LOTUS database (https://lotus.naturalproducts.net/download), the multiclassification of NPs within five kingdoms (animal, bacteria, chromista, fungi, and plant) are performed for 133,092 NPs.

The generated results and figures are presented in the manuscript titled “Composite Machine Learning Strategy for Natural Products Taxonomical Classification and Structural Insights".

## Environment
The python packages in the conda environment for performing the taxonomical classification and structural insights are provided the yaml file (`environment.yml`).

## Directory structure
- `0_data/`: Contains source data used in the analysis.
- `1_dataset_preparation/`: Contains code for data curation and reproduction of figures.
- `2_models/`: Contains retrieved and organized datasets for analysis and modeling.
  + `2.1_GCNN_model/`
  + `2.2_SVM_models/`
- `3_structural_analysis/`: Contains Python scripts for data preprocessing, visualization, and model development.
  + `3.1_tsne_analysis/`
  + `3.2_tsne_comparison/`
  + `3.3_critical_substructures/`

## Instruction
In the folder `1_dataset_preparation/`
- Run the script `1_data_curation.py` to process the `lotusUniqueNaturalProduct.bson` (obtained from LOTUS database) and curate the dataframe `df_curated.csv` for data exploration.
- Run the script `2_data_exploration_preparation.py` to process the `df_curated.csv` and generate the machine learning dataset `np_5classes.csv` for taxonomical classification. The machine learning dataset `np_5classes.csv` consists of the non-isomeric SMILES of 133092 natural products (NPs) from five different kingdoms ([animal, bacteria, chromista, fungi, plant] = [0,1,2,3,4]). The exploratory figures will also be plotted.

In the folder `2_models/2.1_GCNN_model/`
- Run the script `1_train_GCNN_model.py` to train a GCNN model for the machine learning dataset.
- Run the script `2_LOTUS_MPN_lastFFN_fingerprint.py` to generate the MPN (`lotus_MPN.csv`) and last_FFN (`lotus_last_FFN.csv`) fingerprints of 133092 NPs from LOTUS database.
- Run the script `3_NPAtlas_screening.py` to screen 13136 NPs (bacteria and fungi) from NPAtlas dataset (`npatlas_dataset.csv`) using the trained GCNN model; the prediction is `npatlas_preds.csv`. The MPN (`npatlas _MPN.csv`) and last_FFN (`npatlas _last_FFN.csv`) fingerprints of 13136 NPs of are also generated.

In the folder `2_models/2.2_SVM_models/`
- Run the script `1_train_model_SVM_MPN.py` to construct a SVM model (`model_SVM_MPN`) based the MPN fingerprints of the machine learning dataset (`lotus_MPN.csv`).
- Run the script `2_train_model_SVM_last_FFN.py` to construct a SVM model (`model_SVM_last_FFN`) based the MPN fingerprints of the machine learning dataset (`lotus_ last_FFN.csv`).
- Run the script `3_SVM_prediction_NPAtlas.py` to screen 13136 NPs from NPAtlas dataset (`npatlas_dataset.csv`) using the constructed SVM models (`model_SVM_MPN`and `model_SVM_last_FFN`).

In the folder `3_structural_analysis/3.1_tsne_analysis/`
- Run the script `tsne_analysis.py` to plot and compare the Kullback-Leibler (KL) divergence and Davies-Bouldin (DB) score of three different fingerprints (MAP4, MPN and last_FFN) of the 133092 NPs from LOTUS database. The KL divergence and DB score are plotted as function of the perplexity parameter in t-Distributed Stochastic Neighbor Embedding (t-SNE) analysis.

In the folder `3_structural_analysis/3.2_tsne_comparison/`
- Run the script `tsne_map.py` to t-sne mapping of three different fingerprints (MAP4, MPN and last_FFN) of the machine learning dataset 133092 NPs.

In the folder `3_structural_analysis/3.3_critical_substructures/`
- Run the script `critical_substructures.py` to plot the top ten critical substructures of each kingdom from the machine learning dataset 133092 NPs.
