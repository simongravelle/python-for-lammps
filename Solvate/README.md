## Note

This code is in development. Solvate a molecule with a desired solvent.

## How to use

from HydrateLAMMPS import *

lx = 16 # size for the waterbox
ly = 16 # size for the waterbox
lz = 16 # size for the waterbox
d = 2.8 # distance between solute
hydrated = hydrate(SoluteData,SolventData,lx,ly,lz,d)

SoluteData and SolventData can be generated using ReadDataLAMMPS
