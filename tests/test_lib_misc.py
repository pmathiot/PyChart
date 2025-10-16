import pytest
import numpy as np
from lib_misc import get_name

def test_get_name():
    cregex='toto|titi|tata'
    cvarlst=['abcd','titi','efgh']
    assert get_name(cregex,cvarlst)=='titi'

    cregex='toto|titi|tata'
    cvarlst=['abcd','titi','toto']
    assert get_name(cregex,cvarlst)=='titi'
