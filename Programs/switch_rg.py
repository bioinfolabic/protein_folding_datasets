# -*- coding: utf-8 -*-
import numpy as np

# numero de arquivos a serem lidos
n_arquivos = 1000

# criando o vetor para armazenas os rG
vetor_RgAll = []
vetor_RgH = []
vetor_RgP = []

# rodando por todos os arquivos
for i in xrange(1, n_arquivos+1):

    # armazenando as linhas corretas do arquivo
    with open("/home/bruna/heatmap/13_fibonacci_rG/rg"+str(i)+".txt", "r") as a: # pegando o arquivo ############################################################################### CHANGE HERE
        vetor_RgAll = []
        vetor_RgH = []
        vetor_RgP = []
        for line in a: # percorrendo todas as linhas de um arquivo
            if 'rGA' in line:
                vetor_RgAll.append(line.replace("",""))
            if 'rGH' in line:
                vetor_RgH.append(line.replace("",""))
            if 'rGP' in line:
                vetor_RgP.append(line.replace("",""))

    # substituindo os valores errados pelos certos
    with open("/home/bruna/heatmap/pathways/13_fibonacci/pathways13_"+str(i)+".txt", "a") as b:
        with open("/home/bruna/heatmap/pathways_rG_errado/13_fibonacci/ABBABBABABBAB_pathwaystep3000_"+str(i)+"_pathways.txt", "r+w") as c: # pegando o arquivo ##################### CHANGE HERE
                for line in c: # passando de linha em linha procurando os dados no arquivo
                        if 'rGA' in line:
                                b.writelines(vetor_RgAll.pop(0))
                        elif 'rGH' in line:
                                b.writelines(vetor_RgH.pop(0))
                        elif 'rGP' in line:
                                b.writelines(vetor_RgP.pop(0))
                        else:
                                b.writelines(line)
                
