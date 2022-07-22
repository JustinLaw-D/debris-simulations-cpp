# stripped down implementation of NASA breakup model, used to generate initial values

def L_cdf(L, L_min, L_max, typ):
    '''
    calculates the cumulative distribution function for characteristic lengths
    from a collision/explosion at length L, assuming the distribution is truncated 
    at L_min and L_max

    Parameter(s):
    L : characteristic length (m)
    L_min : minimum characteristic length to consider (m)
    L_max : maximum characteristic length to consider (m)
    typ : one of 'coll' (collision) or 'expl' (explosion)

    Keyword Parameter(s): None

    Output(s):
    P : value of CDF at L

    Note(s): returns 0 on an invalid type
    '''

    if typ == 'coll' : beta = -1.71
    elif typ == 'expl' : beta = -1.6
    else:
        print('WARNING: Invalid Debris Generation Type')
        return 0
    return (L_min**beta - L**beta)/(L_min**beta - L_max**beta)