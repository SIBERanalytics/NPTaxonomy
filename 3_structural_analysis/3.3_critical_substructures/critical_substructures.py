# Import necessary packages
import pandas as pd 
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
from IPython.display import display

# Read data
substructures = pd.read_csv('critical_substructures.csv')
labels = ['A','B','C','D','E']

# Plot critical substructures for each kingdom
for (indx, kingdom) in enumerate(substructures):
    # list of SMILES
    smiList = substructures[kingdom]
    # list of labels
    label_list = [labels[indx] + str(l+1) for l in range(len(substructures))]
    # Create RDKit molecular objects
    mols = [Chem.MolFromSmiles(m) for m in smiList]
    
    # display critical substructures
    plt.figure()
    display(Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200),legends=label_list))
    grid_image = Draw.MolsToGridImage(mols, molsPerRow=5, subImgSize=(200, 200),legends=label_list)
    grid_image.save(f'figure_critical_structures_{kingdom}.jpg', dpi=(200, 200))