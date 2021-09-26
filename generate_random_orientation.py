import random
import numpy as np
from scipy.spatial.transform import Rotation as R
def randomorientation(XYZ):
    '''3D aleatory rotation of molecule/particule coordinate'''
    rotation_degrees = random.randint(0,9000)/100
    rotation_radians = np.radians(rotation_degrees)
    rotation_axis = np.array([random.randint(0,100)/100, random.randint(0,100)/100, random.randint(0,100)/100])
    rotation_axis /= np.linalg.norm(rotation_axis)
    rotation_vector = rotation_radians * rotation_axis
    rotation = R.from_rotvec(rotation_vector)
    mol_rotated = rotation.apply(XYZ)    
    return mol_rotated.T
