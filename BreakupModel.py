# stripped down implementation of NASA breakup model, used to generate initial values

import numpy as np

def randL(num, L_min, L_max, typ):
    '''
    generates num random characteristic lengths for debris from a collision/explosion
    Parameter(s):
    num : number of random lengths to generate
    L_min : minimum characteristic length to consider (m)
    L_max : maximum characteristic length to consider (m)
    typ : one of 'coll' (collision) or 'expl' (explosion)
    Keyword Parameter(s): None
    Output(s):
    L : array of random characteristic lengths (m)
    Note(s): returns 0 on an invalid type
    '''

    lam_min, lam_max = np.log10(L_min), np.log10(L_max)
    if typ == 'coll' : beta = -1.71
    elif typ == 'expl' : beta = -1.6
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0
    P = np.random.uniform(size=num) # get random P values
    lam = np.log10(10**(beta*lam_min) - P*(10**(beta*lam_min) - 10**(beta*lam_max)))/beta
    return 10**lam
