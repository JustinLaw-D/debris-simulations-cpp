# creating and saving basic single-celled starlink system

import sys
sys.path.append('./../')

from NCellShell import NCell
import numpy as np
R = 6371 # radius of earth in km
alt = [575, 625] # altitude band of Starlink satellites (km)
V = 4*np.pi*50*(R+600)**2 # volume of band
S_i = [0]
S_di = [0]
D_i = [0]
N_i = 2.5e-8*V
lam = 1000
T = 50
atmosphere = NCell([S_i], [S_di], [D_i], [N_i], [600], alt, [lam], min_lifetime=0)
atmosphere.save('./data/', 'BasicStarlinkTest')