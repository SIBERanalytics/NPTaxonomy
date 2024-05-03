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

# load SVM models
model_SVM_MPN = pickle.load(open('model_SVM_MPN', 'rb'))
model_SVM_last_FFN = pickle.load(open('model_SVM_last_FFN', 'rb'))

# SVM prediction for different fingerprints
pred_SVM_MPN = model_SVM_MPN.predict(npatlas_MPN_fp)
pred_SVM_last_FFN = model_SVM_last_FFN.predict(npatlas_last_FFN_fp)

# Classification metrics
# SVM_MPN
matrix_SVM_MPN = confusion_matrix(npatlas_class, pred_SVM_MPN)
accuracy_SVM_MPN  = matrix_SVM_MPN.diagonal()/matrix_SVM_MPN.sum(axis=1) 
f1_SVM_MPN  = f1_score(npatlas_class, pred_SVM_MPN, average=None)
# SVM_last_FFN
matrix_SVM_last_FFN = confusion_matrix(npatlas_class, pred_SVM_last_FFN)
accuracy_SVM_last_FFN = matrix_SVM_last_FFN.diagonal()/matrix_SVM_last_FFN.sum(axis=1) 
f1_SVM_last_FFN = f1_score(npatlas_class, pred_SVM_last_FFN, average=None)