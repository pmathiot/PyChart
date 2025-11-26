"""
io_utils.py

This module provides utility functions for handling NetCDF files and extracting data. It includes functions for retrieving latitude and longitude variables, determining variable shapes, and loading 2D data from NetCDF files.

Functions:
- get_name(regex, varlst): Matches a regex pattern to a list of variable names.
- get_latlon_var(cfile): Retrieves latitude and longitude variable names from a NetCDF file.
- get_latlon(cfile, offsety=None): Retrieves 2D latitude and longitude data from a NetCDF file.
- get_variable_shape(ncvar): Determines the shape of a NetCDF variable.
- get_dim(cfile, cdir): Retrieves the size of a specific dimension in a NetCDF file.
- get_dims(cfile): Retrieves the sizes of all dimensions (x, y, z, t) in a NetCDF file.
- get_2d_data(cfile, cvar, ktime=0, klvl=0, offsety=None, lmask=True): Loads 2D data from a NetCDF file.
"""

import numpy as np
import netCDF4 as nc
import re
import sys
from .log import debug, info

# ================================= NETCDF ===============================
def get_name(regex, varlst):
    """
    Matches a regex pattern to a list of variable names.

    Parameters:
    - regex (str): Regular expression pattern to match.
    - varlst (list): List of variable names to search.

    Returns:
    - str: The first matching variable name.
    """
    revar = re.compile(r'\b%s\b' % regex, re.I)
    cvar = list(filter(revar.match, varlst))
    if len(cvar) > 1:
        print(f'{regex} name list is longer than 1 or 0')
        print(f'{cvar[0]} is selected')
    if len(cvar) == 0:
        print(f'no match between {regex} and ')
        print(varlst)
        sys.exit(42)
    return cvar[0]

def get_latlon_var(cfile):
    """
    Retrieve the latitude and longitude variable names from a NetCDF file.

    Parameters:
    - cfile (str): Path to the NetCDF file.

    Returns:
    - tuple: Names of the latitude and longitude variables.
    """
    ncid = nc.Dataset(cfile)
    clon = get_name("(glamt|nav_lon.*|lon|longitude)", ncid.variables.keys())
    clat = get_name("(gphit|nav_lat.*|lat|latitude)", ncid.variables.keys())
    ncid.close()
    return clat, clon

def get_latlon(cfile, offsety=None):
    """
    Retrieve 2D latitude and longitude data from a NetCDF file.

    Parameters:
    - cfile (str): Path to the NetCDF file.
    - offsety (int, optional): Offset for the y-dimension.

    Returns:
    - tuple: 2D arrays of latitude and longitude data.
    """
    clat, clon = get_latlon_var(cfile)
    lat2d = get_2d_data(cfile, clat, offsety=offsety)
    lon2d = get_2d_data(cfile, clon, offsety=offsety)
    lon2d[lon2d >= 180] = lon2d[lon2d >= 180.] - 360.
    delta_lon = np.abs(np.diff(lon2d))
    for i, start in enumerate(np.argmax(delta_lon > 180, axis=1)):
        lon2d[i, start + 1:] += 360
    return lat2d, lon2d

def get_variable_shape(ncvar):
    """
    Determine the shape of a NetCDF variable.

    Parameters:
    - ncvar (netCDF4.Variable): NetCDF variable object.

    Returns:
    - str: Shape descriptor (e.g., 'XY', 'XYZ', etc.).
    """
    redimt = re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    redimz = re.compile(r"\b(z|dep|depth|deptht|nav_lev)\b", re.I)
    redimy = re.compile(r"\b(j|y|y_grid_.+|latitude|lat|nj|ny)\b", re.I)
    redimx = re.compile(r"\b(i|x|x_grid_.+|lon|longitude|long|ni|nx)\b", re.I)
    dimlst = ncvar.dimensions
    if (len(ncvar.shape) == 1) and redimx.match(dimlst[0]):
        cshape = 'X'
    elif (len(ncvar.shape) == 1) and redimy.match(dimlst[0]):
        cshape = 'Y'
    elif (len(ncvar.shape) == 2) and redimx.match(dimlst[1]) and redimy.match(dimlst[0]):
        cshape = 'XY'
    elif (len(ncvar.shape) == 3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape = 'XYT'
    elif (len(ncvar.shape) == 3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimz.match(dimlst[0]):
        cshape = 'XYZ'
    elif (len(ncvar.shape) == 4) and redimx.match(dimlst[3]) and redimy.match(dimlst[2]) and redimz.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape = 'XYZT'
    else:
        print('cshape undefined, error')
        print(dimlst, len(ncvar.shape), redimx.match(dimlst[2]), redimy.match(dimlst[1]), redimt.match(dimlst[0]))
        sys.exit(42)
    return cshape

def get_dim(cfile, cdir):
    """
    Retrieve the size of a specific dimension in a NetCDF file.

    Parameters:
    - cfile (str): Path to the NetCDF file.
    - cdir (str): Dimension name ('x', 'y', 'z', or 't').

    Returns:
    - int: Size of the specified dimension.
    """
    dncdim = {'x': re.compile(r"\b(x|nx|x_grid_.+|lon|longitude|long)\b", re.I),
              'y': re.compile(r"\b(y|ny|y_grid_.+|latitude|lat)\b", re.I),
              'z': re.compile(r"\b(z|dep|depth|deptht)\b", re.I),
              't': re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
              }

    ncid = nc.Dataset(cfile)
    cdim = list(filter(dncdim[cdir].match, ncid.dimensions.keys()))
    if len(cdim) > 1:
        print(regex + ' name list is longer than 1; error')
        print(cdim)
        sys.exit(42)
    elif len(cdim) == 0:
        print(f'{cdir} dim in {cfile} is 0.')
        ndim = 0
    else:
        cdim = cdim[0]
        ndim = len(ncid.dimensions[cdim])
    ncid.close()
    return ndim

def get_dims(cfile):
    """
    Retrieve the sizes of all dimensions (x, y, z, t) in a NetCDF file.

    Parameters:
    - cfile (str): Path to the NetCDF file.

    Returns:
    - tuple: Sizes of dimensions (nx, ny, nz, nt).
    """
    nx = get_dim(cfile, 'x')
    ny = get_dim(cfile, 'y')
    nz = get_dim(cfile, 'z')
    nt = get_dim(cfile, 't')
    return nx, ny, nz, nt

def get_2d_data(cfile, cvar, ktime=0, klvl=0, offsety=None, lmask=True):
    """
    Load 2D data from a NetCDF file.

    Parameters:
    - cfile (str): Path to the NetCDF file.
    - cvar (str): Variable name to extract.
    - ktime (int, optional): Time index (default: 0).
    - klvl (int, optional): Level index (default: 0).
    - offsety (int, optional): Offset for the y-dimension.
    - lmask (bool, optional): Whether to apply masking (default: True).

    Returns:
    - numpy.ndarray: 2D data array.
    """
    info(' reading ' + cvar + ' in ' + cfile + ' ...')

    nx, ny, _, _ = get_dims(cfile)
    if not offsety:
        nx, ny, _, _ = get_dims(cfile)
        offsety = ny

    ncid = nc.Dataset(cfile)
    lmask = True
    ncid.set_auto_maskandscale(lmask)
    clvar = get_name(cvar, ncid.variables.keys())
    var = ncid.variables[clvar]

    shape = get_variable_shape(var)

    dslice = {
        'XY': (                                                  slice(0, offsety, None), slice(0, None, None)),
        'XYT': (slice(ktime - 1, ktime, None),                        slice(0, offsety, None), slice(0, None, None)),
        'XYZ': (                          slice(klvl - 1, klvl, None), slice(0, offsety, None), slice(0, None, None)),
        'XYZT': (slice(ktime - 1, ktime, None), slice(klvl - 1, klvl, None), slice(0, offsety, None), slice(0, None, None))
    }

    debug(f"Data shape: {var.shape}, type: {shape}")

    if shape == 'X':
        print(' 1d variable X => extend it 2d')
        tmp = np.zeros(shape=(ny,))
        var, _ = np.meshgrid(var[:], tmp)
        dat2d = var[dslice['XY']].squeeze()
    elif shape == 'Y':
        print(' 1d variable Y => extend it 2d')
        tmp = np.zeros(shape=(nx,))
        _, var = np.meshgrid(tmp, var[:])
        dat2d = var[dslice['XY']].squeeze()
    elif (shape == 'XY') or (shape == 'XYT') or (shape == 'XYZ') or (shape == 'XYZT'):
        dat2d = var[dslice[shape]].squeeze()
    else:
        print(f'{cvar} contains {len(var.shape)} dimensions')
        print(f'dimension names are {var.dimensions}')
        print(f' shape {shape} is unknown, exit ')
        sys.exit(42)

    ncid.close()

    debug(f" {cvar} read: shape = {dat2d.shape} shape expected = ({ny},{nx}) real shape = {dat2d.shape} shape, {dslice[shape]}, {ktime}, {klvl}")

    return dat2d
