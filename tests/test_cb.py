import pytest
import numpy as np
from cb import get_lvl

def test_get_lvl():
    lvls_in=[0,10,3]
    lvls_out=np.array([0.,3.,6.,9.])
    assert np.all(get_lvl(lvls_in) == lvls_out)

    lvls_in=[0,9]
    lvls_out=np.array([0.,1.,2.,3.,4.,5.,6.,7.,8.,9.])
    assert np.all(get_lvl(lvls_in) == lvls_out)

    lvls_in=[0,3,100,200]
    lvls_out=np.array([0.,3.,100.,200.])
    assert np.all(get_lvl(lvls_in) == lvls_out)
