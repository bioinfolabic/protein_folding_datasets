# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# numero de arquivos a serem lidos
n_arquivos = 1
n_dados = 1000
max_step = 3000000
step = 3000
# numero de arquivos a serem lidos
n_arquivos4 = 99
n_dados4 = 1000
max_step4 = 8000000
step4 = 8000

# criando o vetor da energia exponencial
vetor_Ep1 = [] # um vetor dinâmico para a energia exponencial de um arquivo
vetor_Ep2 = []
vetor_Ep3 = []
vetor_Ep4 = []

vetor_all_Ep1 = []
vetor_all_Ep2 = []
vetor_all_Ep3 = []
vetor_all_Ep4 = []

vetor_all_Ep1_media = [] # vetor com a média de todos os vetores
vetor_all_Ep2_media = []
vetor_all_Ep3_media = []
vetor_all_Ep4_media = []

vetor_step = []
vetor_step4 = []

data1 = []
data2 = []
data3 = []
data4 = []

# pegando os valores de cada Step
vetor_step.append(1)
for i in xrange(1, max_step):
        if (i%step==0):
                vetor_step.append(i)
# pegando os valores de cada step1
vetor_step4.append(1)
for i in xrange(2, max_step4):
        if (i%step4==0):
                vetor_step4.append(i)
# pegando as enegias de todos os arquivos
for i in xrange(1, n_arquivos+1):
    print i
    with open("/home/bruna/DM/Proteinas/13FIBO/pathways/ABBABBABABBAB_pathwaystep3000_"+str(i)+"_pathways.txt", "r") as input: ################################################ CHANGE HERE 13FIBO
        vetor_Ep1 = [] # limpando o vetor
        for line in input: # passando de linha em linha procurando os dados no arquivo
            if 'Potential' in line:
                vetor_Ep1.append(float(line.replace("Potential Energy = ","")))
    with open("/home/bruna/DM/Proteinas/2GB1/2GB1/pathways56_"+str(i)+".txt", "r") as input: # ################################################################################# CHANGE HERE 2GB1
        vetor_Ep2 = [] # limpando o vetor
        for line in input: # passando de linha em linha procurando os dados no arquivo
            if 'Potential' in line:
                vetor_Ep2.append(float(line.replace("Potential Energy = ","")))
    with open("/home/bruna/DM/Proteinas/1PLC/pathways/1PCY/pathways99_"+str(i)+".txt", "r") as input: # ################################################################################# CHANGE HERE 1PLC
        vetor_Ep3 = [] # limpando o vetor
        for line in input: # passando de linha em linha procurando os dados no arquivo
            if 'Potential' in line:
                vetor_Ep3.append(float(line.replace("Potential Energy = ","")))
        vetor_all_Ep1.append(vetor_Ep1)
        vetor_all_Ep2.append(vetor_Ep2)
        vetor_all_Ep3.append(vetor_Ep3)
# pegando as enegias de todos os arquivos
for i in xrange(1, n_arquivos4+1):
    with open("/home/bruna/Área de Trabalho/5NAZ/5NAZ_229_pathway_"+str(i)+".txt", "r") as input: ################################################ CHANGE HERE
        vetor_Ep4 = [] # limpando o vetor
        for line in input: # passando de linha em linha procurando os dados no arquivo
                if 'Potential' in line:
                        vetor_Ep4.append(float(line.replace("Potential Energy = ","")))
        vetor_all_Ep4.append(vetor_Ep4)
# fazendo a média de todos os vetores
for k in xrange(0, n_dados):
        soma1 = 0
        soma2 = 0
        soma3 = 0
        for j in xrange(0, n_arquivos):
                soma1 = soma1 + vetor_all_Ep1[j][k] # fazendo a média do primeiro dado de cada arquivo, então temos que passar por todos os arquivos rapidamente
                soma2 = soma2 + vetor_all_Ep2[j][k]
                soma3 = soma3 + vetor_all_Ep3[j][k]
        soma1 = soma1/n_arquivos # assim temos a média dos dados dos n arquivos
        soma2 = soma2/n_arquivos
        soma3 = soma3/n_arquivos
        vetor_all_Ep1_media.append(soma1) # vetor com a média de todos os vetores
        vetor_all_Ep2_media.append(soma2)
        vetor_all_Ep3_media.append(soma3)
for k in xrange(0, n_dados4):
        soma4 = 0
        for j in xrange(0, n_arquivos4):
                print 'j: '+str(j)
                soma4 = soma4 + vetor_all_Ep4[j][k] # fazendo a média do primeiro dado de cada arquivo, então temos que passar por todos os arquivos rapidamente
        soma4 = soma4/n_arquivos4 # assim temos a média dos dados dos n arquivos
        vetor_all_Ep4_media.append(soma4) # vetor com a média de todos os vetores
# print vetor_all_media
'''print vetor_all_Ep1_media[999]
print vetor_all_Ep2_media[999]
print vetor_all_Ep3_media[999]'''
# calculando o desvio padrão do último valor da energia potencial
'''for i in xrange (0, n_arquivos):
        data1.append(vetor_all_Ep1[i][999])
        data2.append(vetor_all_Ep2[i][999])
        data3.append(vetor_all_Ep3[i][999])
v1 = np.var(data1)
v2 = np.var(data2)
v3 = np.var(data3)
d1 = np.sqrt(v1)
d2 = np.sqrt(v2)
d3 = np.sqrt(v3)
print "Desvio padrão = " + str(d)'''
print "fim! :D"
# agora fazendo o gráfico dos dados recolhidos
plt.plot(vetor_step, vetor_all_Ep1_media, label="13 fibonacci")
plt.plot(vetor_step, vetor_all_Ep2_media, label="2GB1")
plt.plot(vetor_step, vetor_all_Ep3_media, label="1PLC")
plt.plot(vetor_step4, vetor_all_Ep4_media, label="5NAZ")
# deixando em notação científica
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
# escrevendo o título e o nome dos eixos em inglês
plt.xlabel('Pathway step (x1e6)', labelpad=10, fontsize=15)
plt.ylabel('Average Ep (x1e2)', labelpad=10, fontsize=15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.grid(True)
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0., fontsize = 15)
plt.rcParams['figure.figsize'] = (11,7)
plt.savefig('/home/bruna/DM/graficos_ep/average_Ep.png', orientation='landscape', dpi=100, bbox_inches='tight')
plt.show()
