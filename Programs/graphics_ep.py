# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# numero de arquivos a serem lidos
n_arquivos1 = 1000
n_dados1 = 1000
max_step1 = 3000000
step1 = 3000
id = 160

# criando o vetor da energia exponencial
vetor_Ep1 = [] # um vetor dinâmico para a energia exponencial de um arquivo

vetor_all_Ep1 = []

vetor_all_Ep1_media = [] # vetor com a média de todos os vetores

vetor_all_Ep1_media_normalizado = []

vetor_step1 = []

data = []

# pegando os valores de cada step1
vetor_step1.append(1)
for i in xrange(2, max_step1):
        if (i%step1==0):
                vetor_step1.append(i)
# pegando as enegias de todos os arquivos
for i in xrange(1, n_arquivos1+1):
    with open("/home/bruna/DM/Proteinas/13FIBO/pathways/ABBABBABABBAB_pathwaystep3000_"+str(i)+"_pathways.txt", "r") as input: ################################################ CHANGE HERE
        vetor_Ep1 = [] # limpando o vetor
        for line in input: # passando de linha em linha procurando os dados no arquivo
                if 'Potential' in line:
                        vetor_Ep1.append(float(line.replace("Potential Energy = ","")))
        vetor_all_Ep1.append(vetor_Ep1)
# fazendo a média de todos os vetores
for k in xrange(0, n_dados1):
        soma1 = 0
        for j in xrange(0, n_arquivos1):
                soma1 = soma1 + vetor_all_Ep1[j][k] # fazendo a média do primeiro dado de cada arquivo, então temos que passar por todos os arquivos rapidamente
        soma1 = soma1/n_arquivos1 # assim temos a média dos dados dos n arquivos
        vetor_all_Ep1_media.append(soma1) # vetor com a média de todos os vetores
max = vetor_all_Ep1_media[990]
print 'max = ' + str(max)
for k in xrange(0, n_dados1):
        vetor_all_Ep1_media_normalizado.append(vetor_all_Ep1_media[k]/(-1*max)) # divide por -1*max

# calculando o desvio padrão do último valor da energia potencial
'''for i in xrange (0, n_arquivos1):
        data.append(vetor_all_Ep1[i][999])
v = np.var(data)
d = np.sqrt(v)
print "Media da enerdia potencial final = " + str(vetor_all_Ep1_media[999])
print "Desvio padrão = " + str(d)
print "fim! :D"'''
# agora fazendo o gráfico dos dados recolhidos
plt.plot(vetor_step1, vetor_all_Ep1_media_normalizado) # tirei , label="5NAZ"
# deixando em notação científica
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
# escrevendo o título e o nome dos eixos em inglês
plt.xlabel('Pathway step (x1e6)', labelpad=10, fontsize=25)
plt.ylabel('Average Ep (x24.92)', labelpad=10, fontsize=25)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.ylim(ymin=-1.05,ymax=0.05)
plt.grid(True)
#plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize = 15)
plt.rcParams['figure.figsize'] = (11,7)
plt.savefig('/home/bruna/DM/graficos_ep/graphic_average_ep_13FIBO.png', orientation='landscape', dpi=100, bbox_inches='tight')
plt.show()
