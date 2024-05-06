# Composite Machine Learning Strategy for Natural Products Taxonomical Classification and Structural Insights
![alt text](https://github.com/SIBERanalytics/NPTaxonomy/blob/main/overview_figure.png?raw=true)

## Overview

This repository contains the source data and python scripts for taxonomical multi-class classification of natural products (NPs) into five kingdoms (animal, bacteria, chromista, fungi, and plant). Training data is obtained from the LOTUS database (https://lotus.naturalproducts.net/download).

The generated results and figures are presented in the manuscript titled “Composite Machine Learning Strategy for Natural Products Taxonomical Classification and Structural Insights".

## Environment
The python packages in the conda environment for performing the taxonomical classification and structural insights are provided the yaml file (`environment.yml`).

## How to use the pre-trained model for taxonomical classification of natural products from SMILES
1. Setup environment (`environment.yml`).
```
conda env create -f environment.yml
conda activate chemenv   
```
2. Go to `4_quickstart/`.
3. Edit (`test_data.csv`) to input the SMILES to be predicted.
4. Run the code under the notebook (`Prediction.ipynb`).
5. Access predictions from (`test_data_preds.csv`).

## Directory structure
- `0_data/`: Contains part of the source data used in the analysis. The rest has been uploaded to Figshare. (https://figshare.com/projects/Composite_Machine_Learning_Strategy_for_Natural_Products_Taxonomical_Classification_and_Structural_Insights/203637)
- `1_dataset_preparation/`: Contains code for data curation and reproduction of figures.
- `2_models/`: Contains retrieved and organized datasets for analysis and modeling.
  + `2.1_GCNN_model/`: Graph Convolutional Neural Network Model
  + `2.2_SVM_models/`: Composite GCNN-SVM (Support Vector Machine) Model
- `3_structural_analysis/`: Contains Python scripts for data processing, visualization, and structural analysis.
  + `3.1_tsne_analysis/`: t-distributed Stochastic Neighbor Embedding (t-SNE) analysis of 3 molecular fingerprints (MAP4, MPN, Last_FFN)
  + `3.2_tsne_comparison/`: Python script to reproduce comparison graphs of the 3 molecular fingerprints
  + `3.3_critical_substructures/`: Python script to reproduce critical substructures of the 5 kingdoms
- `4_quickstart/`: Contains examples for use of pre-trained composite GCNN-SVM model.

## Instructions
In the folder `0_data/`
- machine learning dataset ('np_5classes.csv')
- test dataset ('npatlas_dataset.csv')
  
In the folder `1_dataset_preparation/`
- Run the script `1_data_curation.py` to process the `lotusUniqueNaturalProduct.bson` (10.6084/m9.figshare.25745637) and curate the dataframe `df_curated.csv` (10.6084/m9.figshare.25745325) for data exploration.
- Run the script `2_data_exploration_preparation.py` to process `df_curated.csv` (10.6084/m9.figshare.25745325) and generate the machine learning dataset `np_5classes.csv` for taxonomical classification. The machine learning dataset `np_5classes.csv` consists of the non-isomeric SMILES of 133,092 natural products (NPs) from five different kingdoms ([animal, bacteria, chromista, fungi, plant] = [0,1,2,3,4]). The exploratory figures will also be plotted.

In the folder `2_models/2.1_GCNN_model/`
- Run the script `1_train_GCNN_model.py` to train a GCNN model on the machine learning dataset.
- Run the script `2_LOTUS_MPN_lastFFN_fingerprint.py` to generate the MPN (`lotus_MPN.csv`) and last_FFN (`lotus_last_FFN.csv`) fingerprints (10.6084/m9.figshare.25745448) of 133,092 NPs from the LOTUS database.
- Run the script `3_NPAtlas_screening.py` to classify 13,136 NPs (bacteria and fungi) from the NPAtlas dataset (`npatlas_dataset.csv`) using pre-trained GCNN model, predictions are saved as `npatlas_preds.csv`. The MPN (`npatlas _MPN.csv`) and last_FFN (`npatlas _last_FFN.csv`) fingerprints (10.6084/m9.figshare.25745421) of 13,136 NPs of are also generated.

In the folder `2_models/2.2_SVM_models/`
- Run the script `1_train_model_SVM_MPN.py` to construct a SVM model (`model_SVM_MPN`) (10.6084/m9.figshare.25745634) based the MPN fingerprints of the machine learning dataset (`lotus_MPN.csv`) (10.6084/m9.figshare.25745448).
- Run the script `2_train_model_SVM_last_FFN.py` to construct a SVM model (`model_SVM_last_FFN`) (10.6084/m9.figshare.25745598) based the MPN fingerprints of the machine learning dataset (`lotus_ last_FFN.csv`) (10.6084/m9.figshare.25745448).
- Run the script `3_SVM_prediction_NPAtlas.py` to classify 13,136 NPs from NPAtlas dataset (`npatlas_dataset.csv`) using the constructed SVM models (`model_SVM_MPN` and `model_SVM_last_FFN`).

In the folder `3_structural_analysis/3.1_tsne_analysis/`
- Run the script `tsne_analysis.py` to plot and compare the Kullback-Leibler (KL) divergence and Davies-Bouldin (DB) score of 3 different molecular fingerprints (MAP4, MPN and last_FFN) of the 133,092 NPs from LOTUS database. The KL divergence and DB scores are plotted as function of the perplexity parameter in t-Distributed Stochastic Neighbor Embedding (t-SNE) analysis.

In the folder `3_structural_analysis/3.2_tsne_comparison/`
- Run the script `tsne_map.py` to t-sne mapping of 3 different molecular fingerprints (MAP4, MPN and last_FFN) of the machine learning dataset `np_5classes.csv` of 133,092 NPs.

In the folder `3_structural_analysis/3.3_critical_substructures/`
- Run the script `critical_substructures.py` to plot the top ten critical substructures of each kingdom from the machine learning dataset `np_5classes.csv` of 133,092 NPs.

## Notes
GCNN models are based off Chemprop (https://github.com/chemprop/chemprop)

Literature benchmarking was done using MAP4-SVM models from the manuscript - Capecchi, A., Reymond, JL. "Classifying natural products from plants, fungi or bacteria using the COCONUT database and machine learning" _J. Cheminform_ **13**, 82 (2021).(https://github.com/reymond-group/Coconut-TMAP-SVM)
