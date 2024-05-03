# Import necessary packages
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
plt.rcParams['font.family'] = ['Arial']
plt.rcParams.update({'font.size':14})

# Load data
classes = pd.read_csv('np_5classes.csv')
labels_num = classes.labels

# Kingdom names
kingdom_name = ['Animal','Bacteria', 'Chromista','Fungi','Plant']
# Define colourmap
muted_palette = sns.color_palette('muted')[:5]
reordered_palette = [muted_palette[i] for i in [3, 2, 4, 1, 0]]

# Plots of tsne for different fingerprints
fp = ["MAP4", "MPN", "last_FFN"]
for i in range(len(fp)):
    # Prepare data
    trans_data = pd.read_csv(f'tsne_{fp[i]}.csv')
    data = pd.concat([trans_data,labels_num], axis=1)
    # Plot figure
    plt.figure()
    ax=sns.scatterplot(data=data,x='0',y='1',hue='labels',
                       s=1,palette=reordered_palette,linewidth=0,alpha = 0.5)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, kingdom_name)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(0.95, 1), title='Kingdom', frameon=False )
    plt.gca().set_xticklabels([]),plt.gca().set_yticklabels([])
    plt.tick_params(left = False),plt.tick_params(bottom  = False)
    plt.xlabel('Dimension 1' ),plt.ylabel('Dimension 2' )
    plt.title("$\\mathit{t}$-SNE " + f'({fp[i]})' )
    plt.tight_layout()
    plt.savefig(f'tsne_{fp[i]}.jpg', dpi=300)
    