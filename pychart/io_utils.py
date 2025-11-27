"""
io_utils.py

This module provides utility functions for handling NetCDF files and extracting data. It includes functions for retrieving latitude and longitude variables, determining variable shapes, and loading 2D data from NetCDF files.

Functions:
- get_name(regex, varlst): Matches a regex pattern to a list of variable names.
- get_latlon_var(cfile): Retrieves latitude and longitude variable names from a NetCDF file.
- get_latlon(cfile, offsety=None): Retrieves 2D latitude and longitude data from a NetCDF file.
- get_dims(cfile): Retrieves the sizes of all dimensions (x, y, z, t) in a NetCDF file.
- get_2d_data(cfile, cvar, ktime=0, klvl=0, offsety=None, lmask=True): Loads 2D data from a NetCDF file.
"""

import xarray as xr
import numpy as np
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
    ds = xr.open_dataset(cfile)
    clon = get_name("(glamt|nav_lon.*|lon|longitude)", ds.variables.keys())
    clat = get_name("(gphit|nav_lat.*|lat|latitude)", ds.variables.keys())
    ds.close()
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

    lon2d = lon2d.copy()
    lon2d[lon2d >= 180] -= 360.0

    delta_lon = np.abs(np.diff(lon2d))
    if delta_lon.size:
        for i, start in enumerate(np.argmax(delta_lon > 180, axis=1)):
            if start > 0:
                lon2d[i, start + 1 :] += 360.0
    return lon2d, lat2d

def get_dims(cfile):
    """
    Retrieve the sizes of all dimensions (x, y, z, t) in a NetCDF file.

    Parameters:
    - cfile (str): Path to the NetCDF file.

    Returns:
    - tuple: Sizes of dimensions (nx, ny, nz, nt).
    """
    dncdim = {
        'x': re.compile(r"\b(x|nx|x_grid_.+|lon|longitude|long)\b", re.I),
        'y': re.compile(r"\b(y|ny|y_grid_.+|latitude|lat)\b", re.I),
        'z': re.compile(r"\b(z|dep|depth|deptht)\b", re.I),
        't': re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    }

    ds = xr.open_dataset(cfile)
    sizes = ds.sizes

    nx = next((sizes[dim] for dim in sizes if dncdim['x'].match(dim)), 0)
    ny = next((sizes[dim] for dim in sizes if dncdim['y'].match(dim)), 0)
    nz = next((sizes[dim] for dim in sizes if dncdim['z'].match(dim)), 0)
    nt = next((sizes[dim] for dim in sizes if dncdim['t'].match(dim)), 0)

    ds.close()
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

    ds = xr.open_dataset(cfile, mask_and_scale=lmask)
    var = ds[cvar]

    dncdim = {
        'x': re.compile(r"\b(x|nx|x_grid_.+|lon|longitude|long)\b", re.I),
        'y': re.compile(r"\b(y|ny|y_grid_.+|latitude|lat)\b", re.I),
        'z': re.compile(r"\b(z|dep|depth|deptht)\b", re.I),
        't': re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    }

    if offsety is None:
        offsety = var.sizes.get(next((dim for dim in var.dims if dncdim['y'].match(dim)), 'y'), var.shape[-2])

    if any(dncdim['t'].match(dim) for dim in var.dims):
        var = var.isel({dim: ktime - 1 for dim in var.dims if dncdim['t'].match(dim)})
    if any(dncdim['z'].match(dim) for dim in var.dims):
        var = var.isel({dim: klvl - 1 for dim in var.dims if dncdim['z'].match(dim)})
    if any(dncdim['y'].match(dim) for dim in var.dims):
        var = var.isel({dim: slice(0, offsety) for dim in var.dims if dncdim['y'].match(dim)})

    dat2d = var.values
    ds.close()

    debug(f" {cvar} read: shape = {dat2d.shape}")
    debug(f" offsety = {offsety}")

    return dat2d
