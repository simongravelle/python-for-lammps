import random
def randomlocation(Lx,Ly,Lz):
    '''Choose a random location within a given box
    OLD, see utils.py instead'''
    txlo, txhi = -Lx/2, Lx/2
    tylo, tyhi = -Ly/2, Ly/2
    tzlo, tzhi = -Lz/2, Lz/2    
    x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
    y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
    z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo
    return x, y, z
    
import numpy as np
def generate_random_location(Lx,Ly,Lz):
    '''generate a random location within a given box'''    
    x = np.random.rand()*Lx-Lx/2
    y = np.random.rand()*Ly-Ly/2
    z = np.random.rand()*Lz-Lz/2
    return x, y, z
