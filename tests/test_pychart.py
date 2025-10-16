import pytest
from pychart import get_var_lst, get_file_and_varname

def test_get_var_lst():
    cvar=['toto']
    cfile=['file1.nc','file2.nc','file3.nc']
    assert get_var_lst(cvar,cfile)==['toto','toto','toto']

    cvar=['toto']
    cfile=['file1.nc']
    assert get_var_lst(cvar,cfile)==['toto']

def test_get_file_and_varname():
    cvar=['toto']
    cfile=['file1.nc','file2.nc','file3.nc']
    cdir='datadir/'
    out1,out2=get_file_and_varname(cdir,cfile,cvar)
    assert out1==['datadir/file1.nc','datadir/file2.nc','datadir/file3.nc']
    assert out2==['toto','toto','toto']
