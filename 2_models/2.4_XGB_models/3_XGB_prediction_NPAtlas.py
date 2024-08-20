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

# load XGB models
model_XGB_MPN = pickle.load(open('model_XGB_MPN.pkl', 'rb'))
model_XGB_last_FFN = pickle.load(open('model_XGB_last_FFN.pkl', 'rb'))

# XGB prediction for different fingerprints
pred_XGB_MPN = model_XGB_MPN.predict(npatlas_MPN_fp)
pred_XGB_last_FFN = model_XGB_last_FFN.predict(npatlas_last_FFN_fp)

# Classification metrics
# XGB_MPN
matrix_XGB_MPN = confusion_matrix(npatlas_class, pred_XGB_MPN)
accuracy_XGB_MPN  = matrix_XGB_MPN.diagonal()/matrix_XGB_MPN.sum(axis=1) 
f1_XGB_MPN  = f1_score(npatlas_class, pred_XGB_MPN, average=None)
# XGB_last_FFN
matrix_XGB_last_FFN = confusion_matrix(npatlas_class, pred_XGB_last_FFN)
accuracy_XGB_last_FFN = matrix_XGB_last_FFN.diagonal()/matrix_XGB_last_FFN.sum(axis=1) 
f1_XGB_last_FFN = f1_score(npatlas_class, pred_XGB_last_FFN, average=None)