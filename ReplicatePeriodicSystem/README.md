## Note

This code is in development 

## How to use

from ReplicateSystem import *

datafile = 'mydata.lammps'

Replicate = 2,2,1 # choose how many times to replicate a system along x, y, and z respectively

Init_Atoms,Init_Bonds, Init_Angles, Lx, Ly, Lz = main_replicate(datafile, Replicate)
