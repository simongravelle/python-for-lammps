from copy import copy, deepcopy
from ReadDataLAMMPS import *
import numpy as np

def main_replicate(datafile, Replicate):
	LAMMPSData = readdatalammps(datafile,'full')
	Init_Atoms, Init_Bonds, Init_Angles, Lx, Ly, Lz = initialize(LAMMPSData)
	seq_replicate = set_replicate(Replicate)
	for jj in range(len(seq_replicate)):
	    ReplicatedAtoms, ReplicatedBonds, ReplicatedAngles, newLx, newLy, newLz = ReplicateSystem(Init_Atoms,Init_Bonds,Init_Angles,seq_replicate[jj],Lx,Ly,Lz)
	    del Init_Atoms,Init_Bonds,Init_Angles
	    Init_Atoms,Init_Bonds, Init_Angles, Lx, Ly, Lx = ReplicatedAtoms, ReplicatedBonds, ReplicatedAngles, newLx, newLy, newLz
	return  Init_Atoms,Init_Bonds, Init_Angles, Lx, Ly, Lx

def initialize(LAMMPSData):
	Init_Atoms = deepcopy(LAMMPSData.AtomsProperties.AllAtomsProperties)
	Init_Bonds = deepcopy(LAMMPSData.Bonds)
	Init_Angles = deepcopy(LAMMPSData.Angles)
	Lx = np.float(LAMMPSData.xloxhi[1])-np.float(LAMMPSData.xloxhi[0])
	Ly = np.float(LAMMPSData.yloyhi[1])-np.float(LAMMPSData.yloyhi[0])
	Lz = np.float(LAMMPSData.zlozhi[1])-np.float(LAMMPSData.zlozhi[0])
	
	return Init_Atoms, Init_Bonds, Init_Angles, Lx, Ly, Lz

def set_replicate(Replicate):
	seq_replicate = np.zeros(np.sum(Replicate)-3, dtype=int)
	ii = 0
	for repX in range(Replicate[0]-1): # replicate along x 
		seq_replicate[ii] = 0
		ii += 1
	for repY in range(Replicate[1]-1): # replicate along x 
		seq_replicate[ii] = 1
		ii += 1
	for repZ in range(Replicate[2]-1): # replicate along x 
		seq_replicate[ii] = 2
		ii += 1
	return seq_replicate

def distance(x0, x1, dimensions):
    '''Measure distance between points including boundary conditions'''
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

def CleanBond(NewBonds, InitialBonds, NewAtoms, InitialAtoms, newLx, newLy, newLz, replicate):

    for cptb, bond_id in enumerate(InitialBonds[0]):

        a1i = InitialBonds[2][cptb]
        a2i = InitialBonds[3][cptb]

        a1r = NewBonds[2][cptb]
        a2r = NewBonds[3][cptb]

        c1i = InitialAtoms.T[InitialAtoms[0,:] == a1i][0][4:7]
        c2i = InitialAtoms.T[InitialAtoms[0,:] == a2i][0][4:7]

        c1r = NewAtoms.T[NewAtoms[0,:] == a1r][0][4:7]
        c2r = NewAtoms.T[NewAtoms[0,:] == a2r][0][4:7]

        if replicate == 0:
            di = distance(c1i[0],c2i[0],newLx)
            if di > 5:
                NewBonds[2][cptb] = a1i
                InitialBonds[2][cptb] = a1r
        elif replicate == 1:
            di = distance(c1i[1],c2i[1],newLy)
            if di > 5:
                NewBonds[2][cptb] = a1i
                InitialBonds[2][cptb] = a1r
        elif replicate == 2:
            di = distance(c1i[2],c2i[2],newLz)
            if di > 5:
                NewBonds[2][cptb] = a1i
                InitialBonds[2][cptb] = a1r
                
    return NewBonds, InitialBonds

def CleanAngle(NewAngles, InitialAngles, NewAtoms, InitialAtoms,  newLx, newLy, newLz, replicate):

    for cpta, angle_id in enumerate(InitialAngles[0]):

        a1i = InitialAngles[2][cpta]
        a2i = InitialAngles[3][cpta]
        a3i = InitialAngles[4][cpta]

        a1r = NewAngles[2][cpta]
        a2r = NewAngles[3][cpta]
        a3r = NewAngles[4][cpta]

        c1i = InitialAtoms.T[InitialAtoms[0,:] == a1i][0][4:7]
        c2i = InitialAtoms.T[InitialAtoms[0,:] == a2i][0][4:7]
        c3i = InitialAtoms.T[InitialAtoms[0,:] == a3i][0][4:7]

        c1r = NewAtoms.T[NewAtoms[0,:] == a1r][0][4:7]
        c2r = NewAtoms.T[NewAtoms[0,:] == a2r][0][4:7]
        c3r = NewAtoms.T[NewAtoms[0,:] == a3r][0][4:7]

        if replicate == 0:
            di1 = distance(c1i[0],c2i[0],newLx);
            di2 = distance(c2i[0],c3i[0],newLx);
            if di1 > 5 and di2 < 5:
                InitialAngles[2][cpta] = a1r
                NewAngles[2][cpta] = a1i
            if di2 > 5 and di1 < 5:
                InitialAngles[4][cpta] = a3r
                NewAngles[4][cpta] = a3i
            if di2 > 5 and di1 > 5:
                InitialAngles[3][cpta] = a2r
                NewAngles[3][cpta] = a2i

        if replicate == 1:
            di1 = distance(c1i[1],c2i[1],newLy);
            di2 = distance(c2i[1],c3i[1],newLy);
            if di1 > 5 and di2 < 5:
                InitialAngles[2][cpta] = a1r
                NewAngles[2][cpta] = a1i
            if di2 > 5 and di1 < 5:
                InitialAngles[4][cpta] = a3r
                NewAngles[4][cpta] = a3i
            if di2 > 5 and di1 > 5:
                InitialAngles[3][cpta] = a2r
                NewAngles[3][cpta] = a2i
                
        if replicate == 2:
            di1 = distance(c1i[2],c2i[2],newLz);
            di2 = distance(c2i[2],c3i[2],newLz);
            if di1 > 5 and di2 < 5:
                InitialAngles[2][cpta] = a1r
                NewAngles[2][cpta] = a1i
            if di2 > 5 and di1 < 5:
                InitialAngles[4][cpta] = a3r
                NewAngles[4][cpta] = a3i
            if di2 > 5 and di1 > 5:
                InitialAngles[3][cpta] = a2r
                NewAngles[3][cpta] = a2i
    
    return NewAngles, InitialAngles

def ReplicateSystem(Init_Atoms,Init_Bonds,Init_Angles,replicate,Lx,Ly,Lz):
    
    Natoms = len(Init_Atoms.T)
    Nbonds = len(Init_Bonds.T)
    Nangles = len(Init_Angles.T)
    
    NewAtoms = deepcopy(Init_Atoms)
    NewBonds = deepcopy(Init_Bonds)
    NewAngles = deepcopy(Init_Angles)
    
    newLx = Lx
    newLy = Ly
    newLz = Lz
    
    if replicate == 0:
        NewAtoms[4] += Lx
        newLx *= 2
    elif replicate == 1:
        NewAtoms[5] += Ly
        newLy *= 2
    elif replicate == 2:
        NewAtoms[6] += Lz
        newLz *= 2
        
    NewAtoms[0] += Natoms

    NewBonds[0] += Nbonds
    NewBonds[2] += Natoms
    NewBonds[3] += Natoms

    NewAngles[0] += Nangles
    NewAngles[2] += Natoms
    NewAngles[3] += Natoms
    NewAngles[4] += Natoms
    
    NewBonds, InitialBonds = CleanBond(NewBonds, Init_Bonds, NewAtoms, Init_Atoms, newLx, newLy, newLz, replicate)
    NewAngles, InitialAngles = CleanAngle(NewAngles, Init_Angles, NewAtoms, Init_Atoms, newLx, newLy, newLz, replicate)

    TempAtoms = np.concatenate((Init_Atoms.T,NewAtoms.T))
    TempBonds = np.concatenate((Init_Bonds.T,NewBonds.T))
    TempAngles = np.concatenate((Init_Angles.T,NewAngles.T))

    ReplicatedAtoms = TempAtoms.T
    ReplicatedBonds = TempBonds.T
    ReplicatedAngles = TempAngles.T
    Natoms = len(ReplicatedAtoms.T)
    Nbonds = len(ReplicatedBonds.T)
    Nangles = len(ReplicatedAngles.T)
        
    return ReplicatedAtoms, ReplicatedBonds, ReplicatedAngles, newLx, newLy, newLz
