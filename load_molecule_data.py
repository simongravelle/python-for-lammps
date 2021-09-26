import os
import numpy as np
def loaddata(datapath):
    '''Import all molecule data'''
    P, B, A, D, I = 0, 0, 0, 0, 0
    if os.path.isfile(datapath+'Positions.dat')==True :
        P = np.loadtxt(datapath+'Positions.dat')
    if os.path.isfile(datapath+'Bonds.dat')==True :
        B = np.loadtxt(datapath+'Bonds.dat')
    if os.path.isfile(datapath+'Angles.dat')==True :
        A = np.loadtxt(datapath+'Angles.dat')
    if os.path.isfile(datapath+'Dihedrals.dat')==True :
        D = np.loadtxt(datapath+'Dihedrals.dat')
    if os.path.isfile(datapath+'Impropers.dat')==True :
        I = np.loadtxt(datapath+'Impropers.dat')
    return P, B, A, D, I
