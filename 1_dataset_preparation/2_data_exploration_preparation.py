# Import necessary packages
# Generate seven figures (Figure 2a,2b,3,4,5a,5b,S1)
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt 
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Draw
plt.rcParams['font.family'] = ['Arial']
matplotlib.rcParams.update({'font.size':14})

# Read dataset
data_selected = pd.read_csv('df_curated.csv')
data_columns = data_selected.columns
kingdom_name = data_columns[3:]

# Remove duplicates
duplicated_rows = data_selected[data_selected.iloc[:,1:].duplicated()]
df_cleaned = data_selected.drop_duplicates(subset=data_selected.columns.difference(['lotus_id'])).reset_index(drop=True)

# Check for unique SMILES
# SMILES 2D
CanSMILES2D_unqiue = data_selected['CanSMILES2D'].value_counts().reset_index()
# SMILES 3D
CanSMILES3D_unqiue = data_selected['CanSMILES3D'].value_counts().reset_index()

# Number of kingdoms for each isomeric SMILES
kingdom_sum =  df_cleaned[kingdom_name].sum(axis=1)
counts = kingdom_sum.value_counts()
plt.figure()
ax1 = sns.barplot(x = counts.index, y = counts, palette = 'magma_r')
for bar in ax1.patches:
    ax1.annotate(format(bar.get_height(), '.0f'),
                (bar.get_x() + bar.get_width() / 2,
                 bar.get_height()), ha='center', va='center',
                size=11, xytext=(0, 8),textcoords='offset points', color='saddlebrown')
plt.xticks(fontsize=10)
plt.ylabel('Occurrences'),plt.ylim(0,300000)
plt.xlabel('Number of different kingdoms'),plt.ylabel('Occurrences')
plt.title(str(len(df_cleaned)) + ' isomeric SMILES')
plt.tight_layout()
plt.savefig('figure_isomeric_SMILES_kingdom_number.jpg', dpi=300)

# Structures of isomeric SMILES belonging to 6 or more kingdoms
ind6 = kingdom_sum[kingdom_sum == 6].index.to_series()
ind7 = kingdom_sum[kingdom_sum == 7].index.to_series()
ind6_7 = pd.concat([ind6,ind7],axis=0).reset_index(drop=True)
df_cleaned_ind6_7 = df_cleaned.loc[ind6_7].reset_index(drop=True)
molecule_list = df_cleaned_ind6_7['lotus_id'].tolist()
molecule_names = df_cleaned_ind6_7['CanSMILES3D'].tolist()

plt.figure()
opts = Draw.MolDrawOptions()
opts.legendFraction = 0.20
opts.legendFontSize = 20
mols = [Chem.MolFromSmiles(smi) for smi in molecule_names]
grid_image=Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200),legends=molecule_list, drawOptions=opts, returnPNG=False)
grid_image.save('kingdom_6_7_classes.jpg', dpi=(200, 200))

# Number of isomeric SMILES for each kingdom
df_cleaned_kingdom_single = df_cleaned[(kingdom_sum == 1)]
kingdom_single = df_cleaned_kingdom_single[kingdom_name]
kingdom_single_counts = kingdom_single.sum(axis=0).sort_values(ascending=False)
plt.figure()
ax2=sns.barplot(x = kingdom_single_counts.index, y = kingdom_single_counts, palette = 'muted',)
for bar in ax2.patches:
    ax2.annotate(format(bar.get_height(), '.0f'),
                (bar.get_x() + bar.get_width() / 2,
                 bar.get_height()), ha='center', va='center',
                size=11, xytext=(0, 8),textcoords='offset points', color='saddlebrown')
plt.xticks(fontsize=10)
plt.ylabel('Occurrences'),plt.ylim(0,200000)
plt.title(str(len(df_cleaned_kingdom_single)) + ' isomeric SMILES from single kingdom', fontsize = 14)
plt.tight_layout()
plt.savefig('figure_isomeric_SMILES_kingdom.jpg', dpi=300)


# Analysis of canonicalized 2D SIMLES 
CanSMILES2D_unqiue_kingdom_single = df_cleaned_kingdom_single['CanSMILES2D'].value_counts(sort=False).reset_index()
# Group by CanSMILES2D and calculate the sum directly
CanSMILES2D_unqiue_kingdom_single_sum = df_cleaned_kingdom_single.groupby('CanSMILES2D', sort=False)[kingdom_name].sum()
CanSMILES2D_unqiue_kingdom_single_sum = CanSMILES2D_unqiue_kingdom_single_sum.reset_index()

# Count of canonicalized 2D SIMLES  
single_label_raw = CanSMILES2D_unqiue_kingdom_single_sum
CanSMILES2D_unqiue_kingdom_total = np.count_nonzero(CanSMILES2D_unqiue_kingdom_single_sum.iloc[:,1:], axis=1)
CanSMILES2D_unqiue_kingdom_total = pd.DataFrame(CanSMILES2D_unqiue_kingdom_total)
CanSMILES2D_unqiue_kingdom_total_count = CanSMILES2D_unqiue_kingdom_total[0].value_counts()

plt.figure()
ax3=sns.barplot(x = CanSMILES2D_unqiue_kingdom_total_count.index, y = CanSMILES2D_unqiue_kingdom_total_count, palette = 'magma_r')
for bar in ax3.patches:
    ax3.annotate(format(bar.get_height(), '.0f'),
                (bar.get_x() + bar.get_width() / 2,
                 bar.get_height()), ha='center', va='center',
                size=11, xytext=(0, 8),textcoords='offset points', color='saddlebrown')
plt.xticks(fontsize=10)
plt.ylabel('Occurrences'),plt.ylim(0,150000)
plt.xlabel('Number of different kingdoms'),plt.ylabel('Occurrences')
plt.title(str(len(CanSMILES2D_unqiue_kingdom_total)) + ' non-isomeric SMILES', fontsize = 14)
plt.tight_layout()
plt.savefig('figure_non-isomeric_SMILES_kingdom_number.jpg', dpi=300)


# Counts of canonicalized 2D SIMLES from single kingdom
CanSMILES2D_all_kingdom = CanSMILES2D_unqiue_kingdom_single_sum.iloc[np.where(CanSMILES2D_unqiue_kingdom_total == 1)[0],:]
def find_kingdom_name(row):
    return row.idxmax()
CanSMILES2D_all_kingdom['kingdom'] = CanSMILES2D_all_kingdom.iloc[:,1:].apply(find_kingdom_name, axis=1)
kingdom_counts = CanSMILES2D_all_kingdom['kingdom'].value_counts()

# Figure for the distribution of non-isomeric SMILES for each kingdom
plt.figure()
ax4 = sns.barplot(x = kingdom_counts.index, y = kingdom_counts, palette = 'muted',)
for bar in ax4.patches:
    ax4.annotate(format(bar.get_height(), '.0f'),
                (bar.get_x() + bar.get_width() / 2,
                 bar.get_height()), ha='center', va='center',
                size=11, xytext=(0, 8),textcoords='offset points', color='saddlebrown')
plt.xticks(fontsize=10)
plt.ylabel('Occurrences'),plt.ylim(0,100000)
plt.title(str(len(CanSMILES2D_all_kingdom)) + ' non-isomeric SMILES from single kingdom', fontsize = 14)
plt.tight_layout()
plt.savefig('figure_non-isomeric_SMILES_kingdom.jpg', dpi=300)

# 2D SMILES from same kingdom
selected_smiles2D ='COC(C=C(C)C)CC(C)C1CCC2(C)C3C=CC45OCC3(CCC12C)C4CCC(OC1OC(CO)C(O)C(O)C1O)C5(C)C'
kingdom_same = df_cleaned_kingdom_single[df_cleaned_kingdom_single['CanSMILES2D'] == selected_smiles2D]

mols = [Chem.MolFromSmiles(smi) for smi in kingdom_same['CanSMILES3D']]
molecule_list = kingdom_same['lotus_id'].tolist() 
molecule_list = [x +  '\n(Plant)' for x in molecule_list]

opts = Draw.MolDrawOptions()
opts.legendFraction = 0.20
opts.legendFontSize = 20
grid_image=Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200),legends=molecule_list, drawOptions=opts, returnPNG=False)
grid_image.save('isomeric_SMILES_same_kingdom.jpg', dpi=(200, 200))

# 2D SMILES from different kingdoms
selected_smiles2D ='OCC1OC(OC2C(CO)OC(OC3C(CO)OC(O)C(O)C3O)C(O)C2O)C(O)C(O)C1O'
kingdom_different = df_cleaned_kingdom_single[df_cleaned_kingdom_single['CanSMILES2D'] == selected_smiles2D]

mols = [Chem.MolFromSmiles(smi) for smi in kingdom_different['CanSMILES3D']]
molecule_list = kingdom_different['lotus_id'].tolist() 

from operator import add
kingdom_diff = ['\n(Plant)','\n(Fungi)','\n(Animal)','\n(Bacteria)','\n(Plant)','\n(Plant)','\n(Animal)','\n(Plant)']
molecule_list = list( map(add, molecule_list, kingdom_diff) )

opts = Draw.MolDrawOptions()
opts.legendFraction = 0.20
opts.legendFontSize = 20
grid_image=Draw.MolsToGridImage(mols, molsPerRow=4, subImgSize=(200, 200),legends=molecule_list, drawOptions=opts, returnPNG=False)
grid_image.save('isomeric_SMILES_different_kingdom.jpg', dpi=(200, 200))


# Curate dataset for taxonomical classification
kingdom_selected = ['Animal','Bacteria','Chromista','Fungi','Plant']
CanSMILES2D_5_kingdom = CanSMILES2D_all_kingdom[CanSMILES2D_all_kingdom["kingdom"].isin(kingdom_selected)]
CanSMILES2D_5_kingdom = CanSMILES2D_5_kingdom.reset_index(drop=True)

kingdom_labels = pd.DataFrame({'kingdom': kingdom_selected,
                   'labels': [0, 1, 2, 3,4]})
kingdom_labels.set_index('kingdom', inplace=True)
CanSMILES2D_5_kingdom = CanSMILES2D_5_kingdom.join(kingdom_labels, on="kingdom", how="left", rsuffix="_index")

# Save curated dataset
np_5classes = CanSMILES2D_5_kingdom[['CanSMILES2D','labels']]
np_5classes.columns = ['smiles','labels']
np_5classes.to_csv('np_5classes.csv',index=False)