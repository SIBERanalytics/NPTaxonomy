# Import necessary packages
import pandas as pd
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib
plt.rcParams['font.family'] = ['Arial']
matplotlib.rcParams.update({'font.size':14})

# Import data
metrics_mean = pd.read_csv('results_class_mean.csv') 
metrics_std = pd.read_csv('results_class_std.csv')

# Rename columns
column_names = ['algor','Balanced accuracy','MCC']
metrics_mean.columns = column_names
metrics_std.columns = column_names


# Melt the DataFrame 
melted_metrics_mean = pd.melt(metrics_mean, id_vars='algor', value_vars=['Balanced accuracy','MCC'], var_name='metrics', value_name='score')
melted_metrics_std = pd.melt(metrics_std, id_vars='algor', value_vars=['Balanced accuracy','MCC'], var_name='metrics', value_name='errors')

# Figures for performance comparison
fig = plt.figure()
plt.figure()
ax = sns.barplot(x='metrics', y='score', hue='algor', data=melted_metrics_mean, palette='muted')
x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
y_coords = [p.get_height() for p in ax.patches]
ax.errorbar(x=x_coords, y=y_coords, yerr=melted_metrics_std["errors"], fmt="none", c="k")
for bar in ax.patches:
    ax.annotate(format(bar.get_height(), '.1f'),
                (bar.get_x() + bar.get_width() / 2,
                 bar.get_height()), ha='center', va='center',
                size=11, xytext=(0, 12),textcoords='offset points', color='purple')
plt.xlabel(None),plt.ylabel('Score (%)'),plt.ylim(0,110)
plt.legend(loc='center',bbox_to_anchor=(0.5, 1.075),ncol=3, fontsize=12)
plt.tight_layout()
plt.savefig('figure_performance_comparison.jpg', dpi=300)



