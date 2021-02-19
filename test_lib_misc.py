import pytest
import numpy as np
from lib_misc import get_lvl, get_name

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

def test_get_name():
    cregex='toto|titi|tata'
    cvarlst=['abcd','titi','efgh']
    assert get_name(cregex,cvarlst)=='titi'

    cregex='toto|titi|tata'
    cvarlst=['abcd','titi','toto']
    assert get_name(cregex,cvarlst)=='titi'
