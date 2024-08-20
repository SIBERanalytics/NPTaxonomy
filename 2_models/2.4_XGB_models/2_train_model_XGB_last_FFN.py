# Import necessary packages
import pandas as pd
import xgboost as xgb
import pickle

# Read data
fp = pd.read_csv('lotus_last_FFN.csv')
data = pd.read_csv('np_5classes.csv')
X = fp.iloc[:,1:]
y = data['labels']

# Train model
clf =  xgb.XGBClassifier()
clf.fit(X, y)
        
# Save model        
pickle.dump(clf, open('model_XGB_last_FFN.pkl', "wb"), protocol=4)