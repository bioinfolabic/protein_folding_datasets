# Lucas DF - 13/11/2018 - LABIC
# Cálculo de média e desvio padrão de arquivos de saída da simulação de agregação.
# Versão 2 - Cálculo para BestEnergy, uInter, rGH e rGAll

from os import listdir
from os.path import isfile, join
from sys import argv
from math import pow, sqrt
import numpy as np
import csv

def get_data(mypath, keyword, particles):
    p_energy = 0
    p_energy_inter = 0
    rgh = 0
    rg = 0
    data = []

    for item in particles:
        with open(mypath+'/'+item, 'r') as f:
            for line in f:
                if(' = -nan' in line):
                    print(line)
                    None
                elif(' = 0.000000' in line):
                    print(line)
                    None
                elif('BestEnergy = ' in line):
                    p_energy = float(line[line.find('BestEnergy = ')+len('BestEnergy = '):])
                    print('\np_energy = ', p_energy)
                elif('BestInterEnergy = ' in line):
                    p_energy_inter = float(line[line.find('BestInterEnergy = ')+len('BestInterEnergy = '):])
                    print('p_energy_inter = ', p_energy_inter)
                elif('rGH = ' in line):
                    rgh = float(line[line.find('rGH = ')+len('rGH = '):])
                    print('rgh = ', rgh)
                elif('rGAll = ' in line):
                    rg = float(line[line.find('rGAll = ')+len('rGAll = '):])
                    print('rg = ', rg)
        f.closed
        #if('b' in keyword and p_energy and rgh and rg and len(data) < 10):
        #    print("Energias: ",p_energy, p_energy_inter, rgh, rg)
        #    data.append([p_energy, p_energy_inter, rgh, rg])
        #    print("Data: ", data)
        #elif(p_energy and p_energy_inter and rgh and rg and len(data) < 10):
        #    print("Energias: ",p_energy, p_energy_inter, rgh, rg)
        #    data.append([p_energy, p_energy_inter, rgh, rg])
        #    print("Data: ", data)
        if('b' in keyword and p_energy and rgh and rg and len(data) < 10):
            print("Energias: ",p_energy, p_energy_inter, rgh, rg)
            data.append([p_energy, p_energy_inter, rgh, rg])
            print("Data: ", data)
        elif(p_energy and p_energy_inter and rgh and rg and len(data) < 10):
            print("Energias: ",p_energy, p_energy_inter, rgh, rg)
            data.append([p_energy, p_energy_inter, rgh, rg])
            print("Data: ", data)

        p_energy = 0
        p_energy_inter = 0
        rgh = 0
        rg = 0

    if(len(data) == 10):
        return data
    else:
        print("Erro: quantidade de amostras insuficiente!")
        return

def calc_media(data):

    print('len(data[0]) = ', len(data[0]))

    sum = np.zeros(len(data[0]))
    media = np.zeros(len(data[0]))

    for item in data:
        for x in range(0,len(item)):
            sum[x] += item[x]

    for x in range(0,len(sum)):
        media[x] = sum[x]/len(data)

    return media


def calc_dp(data, media):

    var = np.zeros(len(data[0]))
    dp = np.zeros(len(data[0]))

    for item in data:
        for x in range(0,len(item)):
            var[x] += pow((item[x] - media[x]),2)

    for x in range(0,len(var)):
        dp[x] = sqrt(var[x]/len(data))

    return dp


def main():

    particles = []
    data = []
    lista = []

    for i in range(0,len(argv)):
        print('argv['+str(i)+'] = ', argv[i])

    mypath = argv[1]
    keyword = argv[2]
    fileOut = argv[3]

    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

    for item in files:
        if(keyword in item):
            particles.append(item)

    print('Particles:')
    print('\n'.join(particles))

    data = get_data(mypath, keyword, particles)

    print('Data_main:', data)
    print('len(data_main): ', len(data))

    media = calc_media(data)
    dp = calc_dp(data, media)

    for i in range(0, len(media)):
        lista.append([media[i], dp[i]])

    print('Lista: ', lista)

    print('\nMedia\tDesvio Padrao')
    print('BestEnergy\tBestEnergyInter\trGH\trGAll')
    print(keyword+'\tLista:', lista)

    with open(fileOut+'.csv', 'a', newline='') as f:
        writer = csv.writer(f)
        var = [y for x in lista for y in x]
        var.insert(0, mypath+'_'+keyword)
        writer.writerow(var)
    f.closed

    with open(fileOut+'.txt', 'a') as f:
        f.write(mypath+'_'+keyword+' ')
        # var = str(media[0] + "\t + str(media[1]
        f.write(str(media))
        f.write(str(dp))
        f.write('\n')
    f.closed


if __name__ == "__main__":
    main()
