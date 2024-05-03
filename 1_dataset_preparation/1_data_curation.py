# Import necessary packages
import pandas as pd
import bson
from rdkit import Chem

# Import data from LOTUS database
file ="lotusUniqueNaturalProduct.bson"
with open(file,'rb') as f:
    list_data_lotus = bson.decode_all(f.read())
    
# Convert list to dataframe
df_data_lotus =  pd.DataFrame(list_data_lotus)

# Create curated dataframe
df_curated = df_data_lotus[['lotus_id']]
df_data_lotus['allTaxa'] = df_data_lotus['allTaxa'].astype(str)
 
# Canonicalize 3D SMILES
print('running MolFromSmiles for CanSMILES3D')
mols3D = df_data_lotus['smiles'].progress_apply(Chem.MolFromSmiles)
print('running MolToSmiles for CanSMILES3D')
df_curated.loc[:,'CanSMILES3D'] = mols3D.progress_apply(Chem.MolToSmiles,isomericSmiles=True)

# Canonicalize 2D SMILES
print('running MolFromSmiles for CanSMILES2D')
mols = df_data_lotus['smiles2D'].progress_apply(Chem.MolFromSmiles)
print('running MolToSmiles for CanSMILES2D')
df_curated.loc[:,'CanSMILES2D'] = mols.progress_apply(Chem.MolToSmiles,isomericSmiles=False)

# Identify kingdom 
kingdom_name = ['Animalia','Archaea','Bacteria', 'Chromista','Fungi','Plantae', 'Protozoa']
for i in range(len(kingdom_name)):
    df_curated.loc[:,kingdom_name[i]] = df_data_lotus['allTaxa'].str.contains(kingdom_name[i]).astype(int)   
# rename kingdom
df_curated.rename(columns={'Animalia':'Animal', 'Plantae':'Plant'}, inplace=True)

# Save curated dataframe
df_curated.to_csv('df_curated.csv', index=False)