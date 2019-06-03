import numpy as np
import re
import argparse
import sys

def kabsch_rmsd(P, Q):
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    P = kabsch_rotate(P, Q)
    #print "Kabsch Algoritm Values "
    #print "Pc", Pc
    #print "Qc", Qc
    #print "P", P
    #print "Q", Q
    #print "P rot: ", P
    return rmsd(P, Q)

def kabsch_rotate(P, Q):
    U = kabsch(P, Q)
    #print "U ", U
    # Rotate P
    P = np.dot(P, U)
    return P

def kabsch(P, Q):
    # Computation of the covariance matrix
    C = np.dot(np.transpose(P), Q)

    V, S, W = np.linalg.svd(C)
    #print "========================================"

    #print "V", V
    #print "S", S
    #print "W", W
    #print "========================================"

    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
        #print "here....................................................."
    # Create Rotation matrix U
    U = np.dot(V, W)

    return U

def rmsd(V, W):
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)

def centroid(X):
    C = X.mean(axis=0)
    return C


def get_kabsch(p_all, q_all):
    P = p_all
    Q = q_all 
    return kabsch_rmsd(P, Q) 
    
def main():
    '''
    p_all   = np.asarray([[ 35.000004,35.000014,35.000004],
 [ 34.792181,35.187607,34.039994],
 [ 35.164437,36.005223,33.600754],
 [ 34.217742,35.954461,33.282649],
 [ 34.118338,35.842953,32.293868],
 [ 33.619593,35.259641,31.652773],
 [ 32.812354,35.449049,31.093767],
 [ 32.832539,36.06008, 31.885118],
 [ 33.061227,37.016705,31.704643],
 [ 32.803457,37.288132,30.777343],
 [ 32.494538,37.876455,30.030051],
 [ 32.132192,37.611031,29.136599],
 [ 31.169835,37.666814,29.402603]])

    q_all   = np.asarray([[ 34.969292,35.623087,34.133969],
 [ 34.230047,35.066515,33.754833],
 [ 34.524095,35.852426,33.210885],
 [ 33.849863,35.250414,32.783106],
 [ 33.88853, 36.200478,32.473458],
 [ 33.633721,35.587619,31.725478],
 [ 33.754531,36.52488, 31.398453],
 [ 32.952941,36.395653,31.982196],
 [ 32.986393,37.312947,31.585392],
 [ 32.822086,36.664611,30.841977],
 [ 32.187283,37.404157,30.61815 ],
 [ 32.131615,36.715672,29.895038],
 [ 32.28804, 37.619708,29.49723 ]])
    '''
    p_all = np.asarray([[ -0.98353,1.81095,-0.03140],
[  0.12683,1.80418,-0.03242],
[ -2.05051,4.88314,0.47604],
[ -1.35042,1.15351,0.78475],
[ -1.35043,1.42183,-1.00450],
[ -1.68362,5.27226,1.44914],
[ -3.16088,4.88992,0.47705],
[ -1.68362,5.54059,-0.34012]], dtype=np.float32)

    q_all = np.asarray([[-1.18012,1.83558,-0.02389],
[-0.07891,1.97662,-0.00383],
[-1.87442,3.17118,0.18101],
[-1.47136,1.13188,0.78415],
[-1.47377,1.40501,-1.00437],
[-1.58077,3.60175,1.16149],
[-2.97563,3.03014,0.16095],
[-1.58318,3.87488,-0.62704]], dtype=np.float32)

    P = p_all
    Q = q_all
    #print type(P)
    #print type(Q)
    # Create the centroid of P and Q which is the geometric center of a
    # N-dimensional region and translate P and Q onto that center.
    #Pc = centroid(P)
    #Qc = centroid(Q)
    #P -= Pc
    #Q -= Qc

    #print("Kabsch RMSD: {0}".format(kabsch_rmsd(P, Q)))
    #Kabsch RMSD: 0.842936625152

if __name__ == "__main__":
    main()
    






