# -*- coding: utf-8 -*-
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


# numero de arquivos a serem lidos
n_arquivos = 1000
n_dados = 1000
max_step = 1000000
step = 1000

# criando o vetor da energia exponencial
vetor_RgAll = [] # um vetor dinâmico para a energia exponencial de um arquivo
vetor_RgH = []
vetor_RgP = []

vetor_all_RgAll = []
vetor_all_RgH = []
vetor_all_RgP = []

vetor_all_RgAll_media = [] # vetor com a média de todos os vetores
vetor_all_RgH_media = []
vetor_all_RgP_media = []

vetor_all_RgAll_media_normalizado = []
vetor_all_RgH_media_normalizado = []
vetor_all_RgP_media_normalizado = []

vetor_step = []

vetor_step.append(1)
# pegando os valores de cada Step
for i in xrange(2, max_step):
        if i%step==0 :
                vetor_step.append(i)
# pegando as enegias de todos os arquivos
for i in xrange(1, n_arquivos+1):
    with open("/home/bruna/DM/Proteinas/1PLC/pathways/1PCY/pathways99_"+str(i)+".txt", "r") as input: ################################################################################ CHANGE HERE
        vetor_RgAll = [] # limpando o vetor
        vetor_RgH = []
        vetor_RgP = []
        for line in input: # passando de linha em linha procurando os dados no arquivo
            if 'rGA' in line:
                vetor_RgAll.append(float(line.replace("rGAll = ","")))
            if 'rGH' in line:
                vetor_RgH.append(float(line.replace("rGH = ","")))
            if 'rGP' in line:
               vetor_RgP.append(float(line.replace("rGP = ","")))
        vetor_all_RgAll.append(vetor_RgAll)
        vetor_all_RgH.append(vetor_RgH)
        vetor_all_RgP.append(vetor_RgP)
# fazendo a média de todos os vetores
for k in xrange(0, n_dados):
        soma1 = 0
        soma2 = 0
        soma3 = 0
        for j in xrange(0, n_arquivos):
                soma1 = soma1 + vetor_all_RgAll[j][k] # fazendo a média do primeiro dado de cada arquivo
                soma2 = soma2 + vetor_all_RgH[j][k]
                soma3 = soma3 + vetor_all_RgP[j][k]
        soma1 = soma1/n_arquivos # assim temos a média dos dados dos n arquivos
        soma2 = soma2/n_arquivos
        soma3 = soma3/n_arquivos
        vetor_all_RgAll_media.append(soma1) # vetor com a média de todos os vetores
        vetor_all_RgH_media.append(soma2)
        vetor_all_RgP_media.append(soma3)
max_rgall = vetor_all_RgAll_media[0]
max_rgh = vetor_all_RgH_media[0]
max_rgp = vetor_all_RgP_media[0]
print 'max_rgall = ' + str(max_rgall)
print 'max_rgh = ' + str(max_rgh)
print 'max_rgp = ' + str(max_rgp)
for k in xrange(0, n_dados):
        vetor_all_RgAll_media_normalizado.append(vetor_all_RgAll_media[k]/(max_rgall)) # divide por -1*max
        vetor_all_RgH_media_normalizado.append(vetor_all_RgH_media[k]/(max_rgall)) # divide por -1*max
        vetor_all_RgP_media_normalizado.append(vetor_all_RgP_media[k]/(max_rgall)) # divide por -1*max


dataRgAll = []
dataRgH = []
dataRgP = []
# calculando o desvio padrão do último valor da energia potencial
'''
for i in xrange (0, n_arquivos):
        dataRgAll.append(vetor_all_RgAll[i][999])
        dataRgH.append(vetor_all_RgH[i][999])
        dataRgP.append(vetor_all_RgP[i][999])
vRgAll = np.var(dataRgAll)
vRgH = np.var(dataRgH)
vRgP = np.var(dataRgP)

dRgAll = np.sqrt(vRgAll)
dRgH = np.sqrt(vRgH)
dRgP = np.sqrt(vRgP)

print "Desvio padrão RgAll = " + str(dRgAll)
print "Desvio padrão RgH = " + str(dRgH)
print "Desvio padrão RgP = " + str(dRgP)
print "fim! :D"
# print vetor_all_media
print "fim! :)"
print "Media RgAll = " +str(vetor_all_RgAll_media[999]) # vetor com a média de todos os vetores
print "Media rgH = " + str(vetor_all_RgH_media[999])
print "Media RgP = " + str(vetor_all_RgP_media[999])'''
# agora fazendo o gráfico dos dados recolhidos
plt.plot(vetor_step, vetor_all_RgAll_media_normalizado, label="Average RgAll")
plt.plot(vetor_step, vetor_all_RgH_media_normalizado, label="Average RgH")
plt.plot(vetor_step, vetor_all_RgP_media_normalizado, label="Average RgP")
# colocando limites nos eixos e especificando a notação
plt.axis(ymin = -0.1, ymax = 1.1)
plt.ticklabel_format(axis='both', style='sci', scilimits=(0,0), useMathText=True)
# escrevendo o título e o nome dos eixos em inglês
text = 'Average Rg (x' + str(round(max_rgall, 2)) + ')'
plt.xlabel('Pathway step (x1e6)', labelpad=15, fontsize=20)
plt.ylabel(text, labelpad=15, fontsize=20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.grid(True)
plt.legend(loc="upper right", fontsize=20)
plt.rcParams['figure.figsize'] = (11,7)
plt.savefig('/home/bruna/DM/graficos_rg/rg_1PLC.png', orientation='landscape', dpi=100, bbox_inches='tight')
plt.show()
