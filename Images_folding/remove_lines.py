# -*- coding: utf-8 -*- 

for s in range(1000, 1001):
    i = 0 # i é o número de linhas
    with open("/home/bruna/DM_leandro/pathways/1PCY_errado/pathways99_" + str(s) + ".txt","r") as input:
        with open("/home/bruna/DM_leandro/pathways/1PCY/pathways99_" + str(s) + ".txt","wb") as output: 
            for line in input:
                i = i + 1
                if i < 109109:
                    output.write(line)
                else:
                    break
    print s
