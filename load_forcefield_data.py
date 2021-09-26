import os
import numpy as np
def loadff(datapath):
    '''Import all force field parameters'''
    M, PC, B, A, D, I = 0, 0, 0, 0, 0, 0
    if os.path.isfile(datapath+'Masses.dat')==True :
        M = np.loadtxt(datapath+'Masses.dat')
    if os.path.isfile(datapath+'PairCoeffs.dat')==True :
        PC = np.loadtxt(datapath+'PairCoeffs.dat')
    if os.path.isfile(datapath+'BondCoeffs.dat')==True :
        B = np.loadtxt(datapath+'BondCoeffs.dat')
    if os.path.isfile(datapath+'AngleCoeffs.dat')==True :
        A = np.loadtxt(datapath+'AngleCoeffs.dat')
    if os.path.isfile(datapath+'DihedralCoeffs.dat')==True :
        D = np.loadtxt(datapath+'DihedralCoeffs.dat')
    if os.path.isfile(datapath+'ImproperCoeffs.dat')==True :
        I = np.loadtxt(datapath+'ImproperCoeffs.dat')
    return M, PC, B, A, D, I
