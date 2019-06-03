# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# number of files to be read
n_arquivos1 = 500
n_dados1 = 1000
max_step1 = 8000000
step1 = 8000
id = 160

# creating the vector of potential energy
vetor_Ep1 = [] # a dynamic vector for the potential energy of a file

vetor_all_Ep1 = []

vetor_all_Ep1_media = [] # Vector with the average of all vectors

vetor_all_Ep1_media_normalizado = []

vetor_step1 = []

data = []

# taking the values of each step1
vetor_step1.append(1)
for i in xrange(2, max_step1):
        if (i%step1==0):
                vetor_step1.append(i)
# taking the energies of all the files
for i in xrange(1, n_arquivos1+1):
    with open("/home/bruna/DM/Proteinas/5NAZ/5NAZ/5NAZ_229_pathway_"+str(i)+".txt", "r") as input: ################################################ CHANGE HERE
        vetor_Ep1 = [] # cleaning the vector
        for line in input: # going from line to line searching the data in the file
                if 'Potential' in line:
                        vetor_Ep1.append(float(line.replace("Potential Energy = ","")))
        vetor_all_Ep1.append(vetor_Ep1)
# averaging all vectors
for k in xrange(0, n_dados1):
        soma1 = 0
        for j in xrange(0, n_arquivos1):
                soma1 = soma1 + vetor_all_Ep1[j][k] # averaging the first data of each file, then we have to go through all the files quickly
        soma1 = soma1/n_arquivos1 # so we have the average data of the n files
        vetor_all_Ep1_media.append(soma1) # vector with the mean of all vectors
max = vetor_all_Ep1_media[990]
print 'max = ' + str(max)
for k in xrange(0, n_dados1):
        vetor_all_Ep1_media_normalizado.append(vetor_all_Ep1_media[k]/(-1*max)) # divide por -1*max

# calculating the standard deviation of the last potential energy value
'''for i in xrange (0, n_arquivos1):
        data.append(vetor_all_Ep1[i][999])
v = np.var(data)
d = np.sqrt(v)
print "Media da enerdia potencial final = " + str(vetor_all_Ep1_media[999])
print "Desvio padrão = " + str(d)
print "fim! :D"'''
# now graphing the data collected
plt.plot(vetor_step1, vetor_all_Ep1_media_normalizado)
# leaving in scientific notation
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
# writing the title and names of the axes
plt.xlabel('Pathway step (x1e6)', labelpad=10, fontsize=25)
plt.ylabel('Norm. average Ep', labelpad=10, fontsize=25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.ylim(ymin=-1.05,ymax=0.05)
plt.grid(True)
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize = 15)
plt.rcParams['figure.figsize'] = (11,7)
plt.savefig('/home/bruna/DM/graficos_ep/graphic_average_ep_5NAZ.png', orientation='landscape', dpi=100, bbox_inches='tight')
plt.show()
