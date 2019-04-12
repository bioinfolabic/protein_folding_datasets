# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# configurando o seaborn
sns.set_context("notebook", rc={"axes.labelsize":20}, font_scale=1.5)


'''
# HEATMAP  BEGIN
# carregando o arquivo inicial
a = np.load('/home/bruna/heatmap/generate_dataset/heatmap_13_fibonacci_begin.npy') ########################################################################## CHANGE HERE
# criando os heatmap
ax = sns.heatmap(a/2.5, cmap="YlGnBu", cbar_kws={'label': 'Normalized Kabsch RMSD'}, vmin=0, vmax=1)
# escrevendo o título e o nome dos eixos
plt.title("Initial folding states", size = 20)
plt.axis('off')
# salvando a figura
plt.savefig('heatmap_13_fibonacci_begin.png', orientation='landscape', dpi=100, bbox_inches='tight')
plt.show()


'''
# HEATMAO END
# carregando o arquivo final
b = np.load('/home/bruna/heatmap/generate_dataset/heatmap_13_fibonacci.npy') ################################################################################ CHANGE HERE
# criando os heatmap
bx = sns.heatmap(b/2.5, cmap="YlGnBu", cbar_kws={'label': 'Kabsch RMSD'}, vmin=0, vmax=1)
# escrevendo o título e o nome dos eixos em inglês
plt.title("Final folding states", size = 20)
plt.axis('off')
# salvando a figura
plt.savefig('heatmap_13_fibonacci_end.png', orientation='landscape', dpi=100, bbox_inches='tight')
plt.show()