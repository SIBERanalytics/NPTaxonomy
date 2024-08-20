# Import necessary packages
import pandas as pd
import lightgbm as lgb
import pickle

# Read data
fp = pd.read_csv('0_data/lotus_last_FFN.csv')
data = pd.read_csv('0_data/np_5classes.csv')
X = fp.iloc[:,1:]
y = data['labels']

# Train model
clf = lgb.LGBMClassifier(verbosity=-1)
clf.fit(X, y)
        
# Save model        
pickle.dump(clf, open('model_LGBM_last_FFN.pkl', "wb"), protocol=4)