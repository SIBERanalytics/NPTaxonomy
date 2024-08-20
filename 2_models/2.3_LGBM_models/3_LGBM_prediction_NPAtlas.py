# Import necessary packages
import pandas as pd
import pickle
from sklearn.metrics import confusion_matrix,f1_score

# Read data
npatlas_dataset = pd.read_csv('npatlas_dataset.csv')
npatlas_class = npatlas_dataset['class']
npatlas_MPN = pd.read_csv('npatlas_MPN.csv')
npatlas_last_FFN = pd.read_csv('npatlas_last_FFN.csv')

# Fingerprints
npatlas_MPN_fp = npatlas_MPN.iloc[:,1:]
npatlas_last_FFN_fp = npatlas_last_FFN.iloc[:,1:]

# load LGBM models
model_LGBM_MPN = pickle.load(open('model_LGBM_MPN.pkl', 'rb'))
model_LGBM_last_FFN = pickle.load(open('model_LGBM_last_FFN.pkl', 'rb'))

# LGBM prediction for different fingerprints
pred_LGBM_MPN = model_LGBM_MPN.predict(npatlas_MPN_fp)
pred_LGBM_last_FFN = model_LGBM_last_FFN.predict(npatlas_last_FFN_fp)

# Classification metrics
# LGBM_MPN
matrix_LGBM_MPN = confusion_matrix(npatlas_class, pred_LGBM_MPN)
accuracy_LGBM_MPN  = matrix_LGBM_MPN.diagonal()/matrix_LGBM_MPN.sum(axis=1) 
f1_LGBM_MPN  = f1_score(npatlas_class, pred_LGBM_MPN, average=None)
# LGBM_last_FFN
matrix_LGBM_last_FFN = confusion_matrix(npatlas_class, pred_LGBM_last_FFN)
accuracy_LGBM_last_FFN = matrix_LGBM_last_FFN.diagonal()/matrix_LGBM_last_FFN.sum(axis=1) 
f1_LGBM_last_FFN = f1_score(npatlas_class, pred_LGBM_last_FFN, average=None)