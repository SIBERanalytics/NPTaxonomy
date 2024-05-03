# Import necessary packages
import pandas as pd
from sklearn.svm import SVC
import pickle

# Read data
fp = pd.read_csv('lotus_last_FFN.csv')
data = pd.read_csv('np_5classes.csv')
X = fp.iloc[:,1:]
y = data['labels']

# Train model
clf =  SVC(kernel='linear',probability=True, C = 0.01)
clf.fit(X, y)
        
# Save model        
pickle.dump(clf, open('model_SVM_last_FFN', "wb"), protocol=4)