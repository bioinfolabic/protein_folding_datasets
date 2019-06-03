# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# configuring seaborn
sns.set_context("notebook", rc={"axes.labelsize":20}, font_scale=1.5)



# HEATMAP  BEGIN
# loading the initial file
a = np.load('/home/bruna/DM/generate_dataset/heatmap_5NAZ_begin.npy') ########################################################################## CHANGE HERE
# creating the heatmap
ax = sns.heatmap(a/6, cmap="YlGnBu", cbar_kws={'label': 'Normalized Kabsch RMSD'}, vmin=0, vmax=1) # tirei vmax=1 e a/alguma coisa
ax.figure.axes[-1].yaxis.label.set_size(22)
# writing the title and name of the axes
plt.title("Initial folding states", size = 25)
plt.axis('off')
# saving the figure
plt.savefig('heatmap_5NAZ_begin.eps', orientation='landscape', dpi=100, bbox_inches='tight', metadata='eps')
plt.show()


# HEATMAO END
# loading the final file
b = np.load('/home/bruna/DM/generate_dataset/heatmap_5NAZ_end.npy') ################################################################################ CHANGE HERE
# creating the heatmap
bx = sns.heatmap(b/6, cmap="YlGnBu", cbar_kws={'label': 'Normalized Kabsch RMSD'}, vmin=0, vmax=1) # tirei vmax=1 e b/alguma coisa
bx.figure.axes[-1].yaxis.label.set_size(22)
# writing the title and name of the axes
plt.title("Final folding states", size = 25)
plt.axis('off')
# saving the figure
plt.savefig('heatmap_5NAZ_end.eps', orientation='landscape', dpi=100, bbox_inches='tight', metadata='eps')
plt.show()