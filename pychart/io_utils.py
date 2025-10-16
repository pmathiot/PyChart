import numpy as np
import netCDF4 as nc
import re
import sys

def load_data(file_path):
    # Placeholder for loading data (e.g., from NetCDF)
    return np.random.rand(10, 10)

# ================================= NETCDF ===============================
def get_name(regex,varlst):
    revar = re.compile(r'\b%s\b'%regex,re.I)
    cvar  = list(filter(revar.match, varlst))
    if len(cvar) > 1 :
        print(regex+' name list is longer than 1 or 0')
        print(cvar[0]+' is selected')
    if len(cvar) == 0 :
        print('no match between '+regex+' and :')
        print(varlst)
        sys.exit(42)
    return cvar[0]

def get_latlon_var(cfile):
    ncid   = nc.Dataset(cfile)
    clon=get_name("(glamt|nav_lon.*|lon|longitude)",ncid.variables.keys())
    clat=get_name("(gphit|nav_lat.*|lat|latitude)",ncid.variables.keys())
    ncid.close()
    return clat,clon

def get_latlon(cfile,offsety=None):
    clat,clon=get_latlon_var(cfile)
    lat2d=get_2d_data(cfile,clat,offsety=offsety)
    lon2d=get_2d_data(cfile,clon,offsety=offsety)
    lon2d[lon2d>=180] = lon2d[lon2d>=180.] - 360.
    delta_lon=np.abs(np.diff(lon2d))
    for i, start in enumerate(np.argmax(delta_lon > 180, axis=1)):
        lon2d[i, start+1:] += 360
    return lat2d,lon2d

def get_variable_shape(ncvar):
    redimt=re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    redimz=re.compile(r"\b(z|dep|depth|deptht|nav_lev)\b", re.I)
    redimy=re.compile(r"\b(j|y|y_grid_.+|latitude|lat|nj|ny)\b", re.I)
    redimx=re.compile(r"\b(i|x|x_grid_.+|lon|longitude|long|ni|nx)\b", re.I)
    dimlst = ncvar.dimensions
    if (len(ncvar.shape)==1) and redimx.match(dimlst[0]):
        cshape='X'
    elif (len(ncvar.shape)==1) and redimy.match(dimlst[0]):
        cshape='Y'
    elif (len(ncvar.shape)==2) and redimx.match(dimlst[1]) and redimy.match(dimlst[0]):
        cshape='XY'
    elif (len(ncvar.shape)==3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape='XYT'
    elif (len(ncvar.shape)==3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimz.match(dimlst[0]):
        cshape='XYZ'
    elif (len(ncvar.shape)==4) and redimx.match(dimlst[3]) and redimy.match(dimlst[2]) and redimz.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape='XYZT'
    else:
        print('cshape undefined, error')
        print(dimlst,len(ncvar.shape),redimx.match(dimlst[2]),redimy.match(dimlst[1]),redimt.match(dimlst[0]))
        sys.exit(42)
    return cshape

def get_dim(cfile,cdir):

    dncdim={'x':re.compile(r"\b(x|nx|x_grid_.+|lon|longitude|long)\b", re.I),
            'y':re.compile(r"\b(y|ny|y_grid_.+|latitude|lat)\b", re.I),
            'z':re.compile(r"\b(z|dep|depth|deptht)\b", re.I),
            't':re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
           }

    ncid   = nc.Dataset(cfile)
    cdim=list(filter(dncdim[cdir].match, ncid.dimensions.keys()))
    if len(cdim) > 1 :
        print(regex+' name list is longer than 1; error')
        print(cdim)
        sys.exit(42)
    elif len(cdim) == 0 :
        print(cdir+' dim in '+cfile+' is 0.')
        ndim=0
    else:
        cdim=cdim[0]
        ndim=len(ncid.dimensions[cdim])
    ncid.close()
    return ndim

def get_dims(cfile):
    nx = get_dim(cfile,'x')
    ny = get_dim(cfile,'y')
    nz = get_dim(cfile,'z')
    nt = get_dim(cfile,'t')
    return nx,ny,nz,nt

# get_2d_data
def get_2d_data(cfile,cvar,ktime=0,klvl=0,offsety=None,lmask=True):
    print(' reading '+cvar+' in '+cfile+' ...')
    #if (klvl > 0) and (ktime > 0) :
    #    print('error klvl or ktime larger than 0 (klvl = '+str(klvl)+', ktime = '+str(ktime)+')')
    #    sys.exit(42)

    nx,ny,_,_=get_dims(cfile)
    if not offsety:
        nx,ny,_,_=get_dims(cfile)
        offsety=ny

    ncid   = nc.Dataset(cfile)
    lmask=True
    ncid.set_auto_maskandscale(lmask)
    clvar   = get_name(cvar,ncid.variables.keys())
    var    = ncid.variables[clvar]
    print(var.shape)
    shape  = get_variable_shape(var)

    dslice={
            'XY'  :(                                                  slice(0,offsety,None),slice(0,None,None) ),
            'XYT' :(slice(ktime-1,ktime,None),                        slice(0,offsety,None),slice(0,None,None) ),
            'XYZ' :(                          slice(klvl-1,klvl,None),slice(0,offsety,None),slice(0,None,None) ),
            'XYZT':(slice(ktime-1,ktime,None),slice(klvl-1,klvl,None),slice(0,offsety,None),slice(0,None,None) )
           }

    if shape=='X' :
        print(' 1d variable X => extend it 2d')
        tmp=np.zeros(shape=(ny,))
        var,_=np.meshgrid(var[:],tmp)
        dat2d=var[dslice['XY']].squeeze()
    elif shape=='Y' :
        print(' 1d variable Y => extend it 2d')
        tmp=np.zeros(shape=(nx,))
        _,var=np.meshgrid(tmp,var[:])
        dat2d=var[dslice['XY']].squeeze()
    elif (shape=='XY') or (shape=='XYT') or (shape=='XYZ') or (shape=='XYZT') :
        dat2d=var[dslice[shape]].squeeze()
    else:
        print(cvar+' contains '+str(len(var.shape))+' dimensions')
        print('dimension names are ',var.dimensions)
        print(' shape '+shape+' is unknown, exit ')
        sys.exit(42)

    ncid.close()

    print(' '+cvar+' read: shape = ',dat2d.shape,' shape expected = (',ny,',',nx,')',' real shape = ',dat2d.shape,' ',shape, dslice[shape], ktime, klvl)

    return dat2d
