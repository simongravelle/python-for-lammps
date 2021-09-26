import random
def randomlocation(Lx,Ly,Lz):
    '''Choose a random location within a given box'''
    txlo, txhi = -Lx/2, Lx/2
    tylo, tyhi = -Ly/2, Ly/2
    tzlo, tzhi = -Lz/2, Lz/2    
    x = random.randint(1,1000)/1000*(txhi-txlo)+txlo
    y = random.randint(1,1000)/1000*(tyhi-tylo)+tylo
    z = random.randint(1,1000)/1000*(tzhi-tzlo)+tzlo
    return x, y, z
