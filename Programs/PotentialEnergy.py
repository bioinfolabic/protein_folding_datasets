import numpy as np
from decimal import Decimal

'''
mol = np.asarray(
[[34.969292,35.623087,34.133969],
[34.230047,35.066515,33.754833],
[34.524095,35.852426,33.210885],
[33.849863,35.250414,32.783106],
[33.888530,36.200478,32.473458],
[33.633721,35.587619,31.725478],
[33.754531,36.524880,31.398453],
[32.952941,36.395653,31.982196],
[32.986393,37.312947,31.585392],
[32.822086,36.664611,30.841977],
[32.187283,37.404157,30.618150],
[32.131615,36.715672,29.895038],
[32.288040,37.619708,29.497230]])'''
'''
mol = np.asarray(
[[35.000004,35.000014,35.000004],
[34.792181,35.187607,34.039994],
[35.164437,36.005223,33.600754],
[34.217742,35.954461,33.282649],
[34.118338,35.842953,32.293868],
[33.619593,35.259641,31.652773],
[32.812354,35.449049,31.093767],
[32.832539,36.060080,31.885118],
[33.061227,37.016705,31.704643],
[32.803457,37.288132,30.777343],
[32.494538,37.876455,30.030051],
[32.132192,37.611031,29.136599],
[31.169835,37.666814,29.402603]])'''



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
    
    #print mol
    for n in xrange(mol.shape[0]-2):
        for j in xrange(n+2,mol.shape[0]):
            #print("1) %f\t%f\t%f" % (mol[n,0], mol[n,1], mol[n,2] ))
            r_ij[0] = round(Decimal(mol[n,0]),4) -1. * round(Decimal(mol[j,0]),4)
            r_ij[1] = round(Decimal(mol[n,1]),4) -1. * round(Decimal(mol[j,1]),4)
            r_ij[2] = round(Decimal(mol[n,2]),4) -1. * round(Decimal(mol[j,2]),4)
            
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
            
            uLJ = uLJ + U_LJ_pair;
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
fs = open_file("2GB1_56.txt")
ab_seq = []
sequencia = fs.read()
sequencia = sequencia.replace('\n','')
for  i in sequencia:
	if i == 'n':
		break
	ab_seq.append(i)
ab_seq = np.asarray(ab_seq)
close_file(fs)

print ab_seq

mol = np.asarray([[54.322090,55.691379,50.621181],
[54.289876,55.117862,51.439742],
[55.163920,55.024808,50.962887],
[54.794490,54.180934,51.351989],
[54.974776,54.075277,50.374066],
[54.101161,54.538240,50.523940],
[54.189015,54.120995,49.619403],
[53.873542,55.069929,49.617897],
[54.838093,55.056897,49.881465],
[54.902744,55.931454,49.400872],
[55.590140,55.781887,50.111587],
[55.827933,55.262759,49.290636],
[56.633690,55.171507,49.875811],
[55.837648,54.647706,50.179042],
[56.694873,54.271169,50.530304],
[55.884999,54.056246,51.076116],
[56.282522,53.219318,50.699916],
[55.297815,53.182156,50.870125],
[55.587450,52.552130,50.149582],
[54.764516,53.080607,49.941038],
[55.187172,52.616243,49.162753],
[54.661939,53.439234,48.946385],
[55.643849,53.500227,48.767123],
[55.239560,54.227549,49.321697],
[55.787972,53.580583,49.851490],
[56.364403,54.234352,49.361273],
[57.029989,53.557575,49.675866],
[57.347896,54.289133,49.072740],
[56.783270,55.073223,48.815053],
[55.992115,54.586483,48.444706],
[56.061344,55.557579,48.216274],
[55.347585,55.108289,47.678980],
[54.973803,55.913429,48.139458],
[54.347422,55.178921,47.878408],
[54.096359,55.594049,48.752843],
[54.918915,55.025368,48.751376],
[54.130383,54.417705,48.656701],
[54.971496,54.198302,48.162339],
[54.181203,53.648708,47.891446],
[54.348647,54.308514,47.158904],
[54.902764,53.480840,47.069949],
[55.403762,54.336219,46.938323],
[55.827439,53.759346,47.636688],
[56.289152,54.602881,47.362334],
[56.921433,53.829957,47.309307],
[57.092572,54.513440,48.018926],
[57.589141,53.669384,48.221385],
[56.619213,53.676885,48.464662],
[57.087875,52.805550,48.610023],
[56.317420,52.839767,49.246600],
[56.223491,51.921982,48.860795],
[56.103303,52.666665,48.204287],
[55.371914,52.004537,48.041021],
[55.056037,52.953329,48.037189],
[54.412585,52.293211,47.649617],
[54.082091,52.737100,48.482523]])




potEnergy , uTorsion, uChainAngles, uLJ = penergy(mol, ab_seq)


print potEnergy

'''h = 0
if(h == 0):
    mol = np.asarray([[ 0.,          0.,          0.  ],
 [ 0.60839294, -0.65171681 , 0.45290533],
 [ 0.88366831,  0.27113066 , 0.72230471],
 [ 1.4771602 , -0.44713656 , 1.08542945],
 [ 2.32307844,  0.08412439 , 1.13216418],
 [ 2.56696412, -0.68616842 , 1.72137447],
 [ 2.2477314 ,  0.22413247 , 1.45785454],
 [ 2.30927724, -0.54063186 , 2.09921842],
 [ 1.48995124, -0.05099926 , 1.8009487 ],
 [ 1.82853764, -0.5152943  , 2.6193553 ],
 [ 1.17785722,  0.22870433 , 2.46743091],
 [ 1.87483957,  0.47488078 , 3.14093871],
 [ 1.19819347,  1.18277314  ,2.9383592 ]])
else:
    mol = np.asarray([[  33.964075  , 34.688672  , 30.918651],
[   34.673607  , 33.994754  , 30.796008],
[   34.475531  , 34.226219  , 31.748474],
[   34.689401  , 33.261159  , 31.597078],
[   34.841340  , 33.574751  , 32.534401],
[   34.374818  , 32.694979  , 32.442975],
[   33.804219  , 33.496577  , 32.264482],
[   33.665839   ,32.743351  , 31.621441],
[   33.741849   ,33.631235  , 31.167699],
[   32.860284  , 33.486961  , 31.617173],
[  32.887681  , 34.264960  , 30.989506],
[  33.378591 ,  34.434608  , 31.844040],
[  33.109316,   35.269941  , 31.364759]])

mol = mol + 100
sequence = ["A", "B", "B", "A", "B", "B", "A", "B", "A", "B", "B", "A", "B"]




potEnergy , uTorsion, uChainAngles, uLJ = penergy(mol, sequence)
print "\n\n\n\n"
print "torsion:       ", uTorsion
print "angle:         ", uChainAngles
print "lennard jones: ", uLJ
print "pEnergy:       ", potEnergy

def verify(xyz):
    r       =  np.sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2])
    print "r: ", r


for i in xrange(mol.shape[0]-1):
    mol = mol - mol[i]
    verify(mol[i+1].reshape(3))
'''