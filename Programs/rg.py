# -*- coding: utf-8 -*-
import numpy as np

# se for mudar o programa, mudar o arquivo, a sequência e o número de arquivos

def RgH(mol, sequence):
    RgH = 0.
    x_avg = 0.; y_avg = 0.; z_avg = 0.;
    n = 0
    for i in xrange(mol.shape[0]):
        if(sequence[i] == 'A'):
            x_avg = x_avg + mol[i,0]
            y_avg = y_avg + mol[i,1]
            z_avg = z_avg + mol[i,2]
            n = n + 1;
    x_avg = x_avg / n
    y_avg = y_avg / n
    z_avg = z_avg / n
    
    for i in xrange(mol.shape[0]):
        if(sequence[i] == 'A'):
            RgH = RgH + (mol[i,0] - x_avg) * (mol[i,0] - x_avg) + (mol[i,1] - y_avg) * (mol[i,1] - y_avg) + (mol[i,2] - z_avg) * (mol[i,2] - z_avg)

    RgH = np.sqrt((RgH/n))
    return RgH



def RgP(mol, sequence):
    RgP = 0.
    x_avg = 0.; y_avg = 0.; z_avg = 0.;
    n = 0
    for i in xrange(mol.shape[0]):
        if(sequence[i] == 'B'):
            x_avg = x_avg + mol[i,0]
            y_avg = y_avg + mol[i,1]
            z_avg = z_avg + mol[i,2]
            n = n + 1;
    x_avg = x_avg / n
    y_avg = y_avg / n
    z_avg = z_avg / n
    
    for i in xrange(mol.shape[0]):
        if(sequence[i] == 'B'):
            RgP = RgP + (mol[i,0] - x_avg) * (mol[i,0] - x_avg) + (mol[i,1] - y_avg) * (mol[i,1] - y_avg) + (mol[i,2] - z_avg) * (mol[i,2] - z_avg)

    RgP = np.sqrt((RgP/n))
    return RgP




def RgAll(mol, sequence):
    RgAll = 0.
    x_avg = 0.; y_avg = 0.; z_avg = 0.;
    n = 0
    for i in xrange(mol.shape[0]):
        x_avg = x_avg + mol[i,0]
        y_avg = y_avg + mol[i,1]
        z_avg = z_avg + mol[i,2]
        n = n + 1;
    x_avg = x_avg / n
    y_avg = y_avg / n
    z_avg = z_avg / n
    
    for i in xrange(mol.shape[0]):
        RgAll = RgAll + (mol[i,0] - x_avg) * (mol[i,0] - x_avg) + (mol[i,1] - y_avg) * (mol[i,1] - y_avg) + (mol[i,2] - z_avg) * (mol[i,2] - z_avg)

    RgAll = np.sqrt((RgAll/n))
    return RgAll
'''
rgall = RgAll(mol, sequence)
rgp = RgP(mol, sequence)
rgh = RgH(mol, sequence)

print "rGH:   ", rgh
print "rGP:   ",rgp
print "rGALL: ", rgall'''

def rg(mol, sequece):
    rgh   = RgH(mol, sequece)
    rgp   = RgP(mol, sequece)
    rgall = RgAll(mol, sequece)

    #print "rGH:   ", rgh
    #print "rGP:   ",rgp
    #print "rGALL: ", rgall
    return rgall, rgp, rgh


##########################################################################################################################################################################################
########################################################################################## CHANGE HERE ###################################################################################
##########################################################################################################################################################################################
n_arquivos = 1000 ########################################################################################################################################################################
s = "ABAAAAABBABAAAAABBABABAABBAAABBBAAAABBAAABBBBAABAABABBABABBBBAABABABBABAAABBBABBBABABABBAAAAAABABAB"#################################################################################
sequence = list(s)
print len(sequence)
for i in xrange(1, n_arquivos+1):
    # abrindo o arquivo para salvar os rg nele
    with open("/home/bruna/heatmap/1PCY_rG/rg"+str(i)+".txt", "w") as saida: # pegando o arquivo para escrever ##########################################################################
        flag = 1
        mol = [] # limpando meu arquivo
        # extraindo os dados do rg
        with open("/home/bruna/heatmap/1PCY_99/pathways99_"+str(i)+".txt", "r") as f: # pegando o arquivo para ler ######################################################################
            mol = [] # limpando o vetor da sequência
            for line in f.readlines():
                if 'N	x	y	z' in line:
                    flag = 1
                elif len(line) == 1 and len(mol) != 0:
                    flag = 0
                    rgall, rgp, rgh = rg( np.array(mol), sequence )
                    saida.writelines("rGAll = " + str(rgall) + "\n")
                    saida.writelines("rGH = " + str(rgh) + "\n")
                    saida.writelines("rGP = " + str(rgp) + "\n")
                    mol = []
                elif flag == 1:
                    aux = []
                    for entry in line.split():
                        aux.append(float(entry))
                    aux.pop(0)
                    mol.append(aux)
                
        
