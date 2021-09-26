import copy
import numpy as np
def neighborsearch(neighbor,molecule,cptatm, x, y, z, Lx, Ly, Lz):
    '''Search neighbor in a box and return the closest distance found'''
    mind = 10
    XYZ = copy.deepcopy(neighbor.T[0:cptatm].T)
    for xx in [-Lx,0,Lx]:
        for yy in [-Ly,0,Ly]:
            for zz in [-Lz,0,Lz]:
                XYZ = np.append(XYZ.T,neighbor.T[0:cptatm]+np.array([xx,yy,zz]), axis=0).T
    XYZ = XYZ.T[cptatm:].T
    for m in molecule.T:
        d = np.sqrt((m[0]+x-XYZ[0])**2+(m[1]+y-XYZ[1])**2+(m[2]+z-XYZ[2])**2)
        mind = np.min([np.min(d),mind])
    return mind
