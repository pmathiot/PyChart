import numpy as np

def get_lvls(bnds):
    """ 
    Compute an array of discrete levels.

    Parameters:
    - bnds (list): List of levels with length 2, 3, or more.

    Returns:
    - np.array: Array of levels.

    Raises:
    - SystemExit: If the levels are not properly defined.
    """
    if len(bnds)==2:
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvl=np.linspace(lvlmin, lvlmax, num=10)
    elif len(bnds)==3 :
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvlint = bnds[2]
        lvlmax = lvlmin+round((lvlmax-lvlmin)/lvlint)*lvlint
        lvl= np.arange(lvlmin,lvlmax+0.000001,lvlint)
    elif len(bnds) > 3:
        lvl=np.array(bnds[:])
    else:
        print(' Need definition of levels (min,max) at least.')
        sys.exit(42)
    return lvl

