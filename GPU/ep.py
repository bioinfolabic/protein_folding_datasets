# -*- coding: utf-8 -*-
import numpy as np
from decimal import Decimal



def computeAngleForces(mol):
    
    c12 = 0. 
    uChainAngles = 0.   
    dr1 = np.ndarray(shape=[3])
    dr2 = np.ndarray(shape=[3])

    for i in xrange(mol.shape[0]-2):
        dr1[0]  = mol[i+1,0] - mol[i,0]
        dr1[1]  = mol[i+1,1] - mol[i,1]
        dr1[2]  = mol[i+1,2] - mol[i,2]

        dr2[0]  = mol[i+2,0] - mol[i+1,0]
        dr2[1]  = mol[i+2,1] - mol[i+1,1] 
        dr2[2]  = mol[i+2,2] - mol[i+1,2] 
        
        #c11 = dr1x * dr1x + dr1y * dr1y + dr1z * dr1z    
        c12 = (dr1[0] * dr2[0]) + (dr1[1] * dr2[1]) + (dr1[2] * dr2[2])
        #c22 = dr2x * dr2x + dr2y * dr2y + dr2z * dr2z    
         
        uChainAngles = uChainAngles + c12
    return uChainAngles
#uChainAngles = computeAngleForces(mol)

def computeTorsionForces(mol): 
    dr1 = np.ndarray(shape=[3])
    dr2 = np.ndarray(shape=[3])
    dr3 = np.ndarray(shape=[3])
    uTorsion = 0.
    c11=0.; c12=0.; c13=0.; c22=0.; c23=0.; c33=0.; pi=0.; qia=0.; qib=0.;

    for i in xrange(0,mol.shape[0]-3):

        dr1[0] = mol[i + 1,0] - mol[i,0]
        dr1[1] = mol[i + 1,1] - mol[i,1]
        dr1[2] = mol[i + 1,2] - mol[i,2]

        dr2[0] = mol[i + 2,0] - mol[i + 1,0]
        dr2[1] = mol[i + 2,1] - mol[i + 1,1]
        dr2[2] = mol[i + 2,2] - mol[i + 1,2]

        dr3[0] = mol[i + 3,0] - mol[i + 2,0]
        dr3[1] = mol[i + 3,1] - mol[i + 2,1]
        dr3[2] = mol[i + 3,2] - mol[i + 2,2]

        #c11 = dr1[0] * dr1[0] + dr1[1] * dr1[1] + dr1[2] * dr1[2]
        #c12 = dr1[0] * dr2[0] + dr1[1] * dr2[1] + dr1[2] * dr2[2]    
        #c22 = dr2[0] * dr2[0] + dr2[1] * dr2[1] + dr2[2] * dr2[2]
        c13 = dr1[0] * dr3[0] + dr1[1] * dr3[1] + dr1[2] * dr3[2]    
        #c23 = dr2[0] * dr3[0] + dr2[1] * dr3[1] + dr2[2] * dr3[2]
        #c33 = dr3[0] * dr3[0] + dr3[1] * dr3[1] + dr3[2] * dr3[2]      

        uTorsion = uTorsion + (-0.5) * c13
    return uTorsion


#uTorsion = computeTorsionForces(mol)

def computeLennard_Jones(mol, sequence):
    uLJ      = 0.
    r_ij     = np.ndarray(shape=[3],dtype=float)
    r2       = 0.
    forceLJ  = 0. 
    U_LJ_pair= 0.
    soma = 0
    
    #print mol
    for n in xrange(mol.shape[0]-2):
        soma = 0
        for j in xrange(n+2,mol.shape[0]):
            
            #print("1) %f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2] ))
            r_ij[0] = round(Decimal(mol[n,0]),6) -1. * round(Decimal(mol[j,0]),6)
            r_ij[1] = round(Decimal(mol[n,1]),6) -1. * round(Decimal(mol[j,1]),6)
            r_ij[2] = round(Decimal(mol[n,2]),6) -1. * round(Decimal(mol[j,2]),6)
            
            #print("2) %f\t%f\t%f\t%f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2], r_ij[0], r_ij[1], r_ij[2] ))
            r2        = (r_ij[0] * r_ij[0]) + (r_ij[1] * r_ij[1]) + (r_ij[2] * r_ij[2])
            
            forceLJ   = 24. * ((2. * pow(r2, -7)) -1. * pow(r2, -4))
            U_LJ_pair = 4. * ((pow(r2, -6)) -1. * pow(r2, -3))
            
            #print("3) %f\t%f\t%f\t%f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2], forceLJ, U_LJ_pair, r2 ))

            if((sequence[n] == 'A' and sequence[j] == 'B') or (sequence[n] == 'B' and sequence[j] == 'A') or (sequence[n] == 'B' and sequence[j] == 'B')):
                forceLJ   = (forceLJ * 0.5)
                U_LJ_pair = (U_LJ_pair * 0.5)

            #print("4) %f\t%f\t%f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2], forceLJ, U_LJ_pair ))

            #mol[n,0] = np.float64(mol[n,0]) + ((r_ij[0] * forceLJ))
            #mol[n,1] = np.float64(mol[n,1]) + ((r_ij[1] * forceLJ))
            #mol[n,2] = np.float64(mol[n,2]) + ((r_ij[2] * forceLJ))

            #print("5) %f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2] ))
            #mol[j,0] = np.float64(mol[j,0]) - ((r_ij[0] * forceLJ))
            #mol[j,1] = np.float64(mol[j,1]) - ((r_ij[1] * forceLJ))
            #mol[j,2] = np.float64(mol[j,2]) - ((r_ij[2] * forceLJ))

            #print("\n\nn:%d\tj:%d\tuLJ:%.3f\tforceLJ:%.3f\tr2:%.3f\tU_LJ_pair:%.3f" % (n, j, uLJ, forceLJ, r2, U_LJ_pair))
            #print("%f\t%f\t%f" % (r_ij[0], r_ij[1], r_ij[2] ))
            #print("6) %f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2] ))
            #print("\n\n\n")
            #r_ij = [0.,0.,0.]
            #r2   = 0.
            soma = soma + U_LJ_pair;
            uLJ = uLJ + U_LJ_pair;
            #if(n == 70):
            #    print "i: "+str(n)+"    j: "+str(j)+"      d_uLJComp["+str(j)+"] = "+str(U_LJ_pair)
        #print "i:"+str(n)+" soma = "+str(soma)        
    #print "uLJ: "+str(uLJ)
    return uLJ

#uLJ = computeLennard_Jones(mol,sequence)


#print "torsion: ", uTorsion
#print "angle: ", uChainAngles
#print "lennard jones: ", uLJ
#potEnergy = uTorsion + uChainAngles + uLJ
#print "pEnergy: ", potEnergy

def penergy(mol, sequence):
    uTorsion     = computeTorsionForces(mol)
    uChainAngles = computeAngleForces(mol)
    uLJ          = computeLennard_Jones(mol,sequence)
    
    potEnergy = uTorsion + uChainAngles + uLJ
    #print "torsion:       ", uTorsion
    #print "angle:         ", uChainAngles
    #print "lennard jones: ", uLJ
    #print "pEnergy:       ", potEnergy
    return potEnergy , uTorsion, uChainAngles, uLJ
    #return potEnergy

def open_file(filename):
	f = open(filename, 'r')
	return f 
def close_file(f):
	f.close()



##########################################################################################################################################################################################
########################################################################################## CHANGE HERE ###################################################################################
##########################################################################################################################################################################################
n_arquivos = 1 ########################################################################################################################################################################
s = "BAABAAAABBBBBBBBAABAABABABABBAABAABABABBBABABBAABAABAABBABBAAAAAABABBAABAABBBBBBBAABBABAAAABABAABABBABAAABBAAAABAAAAAAAABBBBABAABAABAABBAAAABBAAABBBAAABABABAAABAABAABBABBAAAABABABABABBBABBBBAAAABABABBAABBABBBBABAABABBBABBABAAABBB"#################################################################################
#s = "ABBABBABABBAB"
#s = "ABBBAAABABBABABBBBBAABAABABBAABBBABBBAABABABBBBABBBABABB"
sequence = list(s)
print len(sequence)

for i in xrange(0, n_arquivos):
    # abrindo o arquivo para salvar os rg nele
    with open("/home/users/brunap/DM_lucas/Proteinas/5NAZ/5NAZ_229_ep_7.txt", "w") as saida: # pegando o arquivo para escrever ##########################################################################
        flag = 1
        mol = [] # limpando meu arquivo
        # extraindo os dados do rg
        with open("/home/users/brunap/DM_lucas/Proteinas/5NAZ/5NAZ_229_pathway_7.txt", "r") as f: # pegando o arquivo para ler ######################################################################
            mol = [] # limpando o vetor da sequência
            for line in f.readlines():
                if 'N	x	y	z' in line:
                    flag = 1
                elif len(line) == 1 and len(mol) != 0:
                    flag = 0
                    #print type(mol)
                    mol = np.array(mol)
                    #print type(mol)
                    #print mol.shape
                    #print type(mol[0,0])
                    #mol = np.array(mol, dtype = np.float32)
                    #print type(mol[0,0])
                    #print mol.shape
                    potEnergy , uTorsion, uChainAngles, uLJ = penergy(mol, sequence)
                    saida.writelines("Potential Energy = " + str(potEnergy) + "\n")
                    saida.writelines("uLJ = " + str(uLJ) + "\n")
                    saida.writelines("uTorsion = " + str(uTorsion) + "\n")
                    saida.writelines("uChainAngles = " + str(uChainAngles) + "\n\n\n")
                    mol = []
                elif flag == 1:
                    aux = []
                    for entry in line.split():
                        aux.append(float(entry))
                    aux.pop(0)
                    mol.append(aux)

