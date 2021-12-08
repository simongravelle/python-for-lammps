import copy
import numpy as np
from numpy.linalg import norm
def neighborsearch(neighbor,molecule,cptatm, x, y, z, Lx, Ly, Lz):
    '''Search all neighbor to a molecule in a box and return the closest distance'''
    box = np.array([Lx, Ly, Lz])
    minr = 10
    for m in molecule.T:
        x0 = m[0] + x
        y0 = m[1] + y
        z0 = m[2] + z
        dxdydz = np.remainder(neighbor[:cptatm].T - np.array([x0,y0,z0]) + box/2., box) - box/2.
        minr = np.min([minr,np.min(norm(dxdydz,axis=1))])
    return minr
