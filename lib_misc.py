import re
import cartopy
import subprocess
import sys

import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# ============================ output argument list in txt file ================
def save_output(cfile, fig):
    """
    save command line as string in a text file.
    save figure.

    Parameters
    ----------
    parameter 1: string
        file name used for the figure and command line record (without extension)

    Returns
    -------
    None

    Raises
    ------
    No raise
    """

    # save argument list
    print(cfile)
    fid = open(cfile+'.txt',"w")
    fid.write(' python '+subprocess.list2cmdline(sys.argv[0:])+'\n')
    fid.close()

    # save figure
    fig.savefig(cfile+'.png', format='png', dpi=150)

# ============================ plot utility ========================================
def add_land_features(ax,cfeature_lst):
    dfeature={'isf':cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none'),
              'lakes':cartopy.feature.NaturalEarthFeature('physical', 'lakes'                    , '50m', facecolor='none'),
              'coast':cartopy.feature.NaturalEarthFeature('physical', 'coastline'                , '50m', facecolor='none'),
              'land' :cartopy.feature.NaturalEarthFeature('physical', 'land'                     , '50m', facecolor='none'),
#              'coast':cartopy.feature.NaturalEarthFeature('physical', 'coastline'                , '50m', facecolor='0.75'),
#              'land' :cartopy.feature.NaturalEarthFeature('physical', 'land'                     , '50m', facecolor='0.75'),
              'bathy_z1000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_J_1000'  , '10m', facecolor='none'),
              'bathy_z2000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_I_2000'  , '10m', facecolor='none'),
              'bathy_z3000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_H_3000'  , '10m', facecolor='none')
             }

    for _,cfeat in enumerate(cfeature_lst):
        ax.add_feature(dfeature[cfeat],linewidth=0.5,edgecolor='k')

def def_projection(proj_name):
    dproj={
           'ortho_natl'  :[ ccrs.Orthographic(central_longitude=-60.0, central_latitude=45.0)   , \
                            [(-180, 180, -90, -60),ccrs.PlateCarree()] ],
           'south_stereo':[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)   , \
                            [(-180, 180, -90, -60),ccrs.PlateCarree()] ],
           'south_ocean' :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)   , \
                            [(-180, 180, -90, -45),ccrs.PlateCarree()] ],
           'ant'         :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)   , \
                            [(-180, 180, -90, -65),ccrs.PlateCarree()] ],
           'arctic'      :[ ccrs.Stereographic(central_latitude= 90.0, central_longitude=0.0)   , \
                            [(-180, 180, 60, 90)  ,ccrs.PlateCarree()] ],
           'ross'        :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=-180.0), \
                            [(-6.67e5,8.33e5,1.05e6,2.47e6), 'cproj' ] ],
           'fris'        :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0), \
                            [(-1.718e6,-4.114e5,1.259e5,1.411e6), 'cproj' ] ],
           'pig'         :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)   , \
                            [(  -96, -105, -73.9, -76),ccrs.PlateCarree()] ],
           'ispig'       :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)   , \
                            [(-1.756e6, -1.342e6, -5.439e5, -1.299e5), 'cproj' ] ],
           'amu'         :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)   , \
                            [(  -99, -130, -70, -78),ccrs.PlateCarree()] ], 
           'global'         :[ ccrs.PlateCarree()                    , ['global'] ],
           'global_robinson':[ ccrs.Robinson(central_longitude=0)    , ['global'] ],
           'global_mercator':[ ccrs.Mercator(central_longitude=-90.0), ['global'] ]
          }
#    elif proj_name=='natl' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[-4.039e6,2.192e6,-1.429e6,4.805e6]
#    elif proj_name=='greenland' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[-1.124e6,0.897e6,1.648e6,5.198e6]
#    elif proj_name=='ovf' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[-3.553e5,2.141e6,9.915e5,3.4113e6]
#    elif proj_name=='ovf_larger' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[-2.5e5,2.45e6,0.9e6,3.6e6]
#    elif proj_name=='irminger' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[-2.164e5,1.395e6,1.635e6,3.265e6]
#    elif proj_name=='japan' :
#        proj=ccrs.LambertConformal(150, 30,cutoff=10)
#        XY_lim=[-3.166e6,2.707e6,-1.008e6,4.865e6]
#    elif proj_name=='feroe' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[1.56e6,2.145e6,1.973e6,2.555e6]
#    elif proj_name=='gulf_stream' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[(-4.250e6,1.115e5,-1.546e6,2.8155e6), proj]
#    else:
#        print('projection '+proj_name+' unknown')
#        print('should be ross, gulf_stream, feroe, global_mercator, global_robinson, japan'
#              ', ovf, greenland, natl, global, south_stereo, ant')
#        sys.exit(42)

    proj=dproj[proj_name][0]

    XY_lim=dproj[proj_name][1]
    if XY_lim[-1]=='cproj':
        XY_lim[-1]=proj

    return proj, XY_lim

# ============================ LEGEND ==================================
def get_corner(ax):
    """
    find the coordinate of the axes

    Parameters
    ----------
    parameter 1: axe

    Returns
    -------
    output 1: list
        position box [x left, x right, y bottom, y top]

    Raises
    ------
    No raise
    """

    box=ax.get_position()
    return [box.x0,box.x1,box.y0,box.y1]

# ======================= COLORBAR =======================================
def get_plt_bound(ax_lst,nplt):
    """
    find the box containing of a list of axes

    Parameters
    ----------
    parameter 1: list of axes

    Returns
    -------
    output 1: list
        position box [x left, x right, y bottom, y top]

    Raises
    ------
    No raise
    """
    x0=1.0; x1=0.0; y0=1.0; y1=0.0
    for iplt in range(0,nplt):
        ax_lst[iplt].apply_aspect()
        box = get_corner(ax_lst[iplt])
        x0=np.min([x0,box[0]])
        x1=np.max([x1,box[1]])
        y0=np.min([y0,box[2]])
        y1=np.max([y1,box[3]])
    return [x0, y0, x1, y1]

# ======================= TITLE =======================================
def get_subplt_title(args,nplt):
    # supplot title prefix
    csubplt_title=['']
    if nplt > 1:
        csubplt_title=['a) ','b) ','c) ','d) ','e) ','f) ','g) ','h) ']

    # get title list for each subplot
    if args.spfid:
        crun_title = args.spfid[:]

    # get subtitle extention
    if (args.mapreff or args.cntreff) and args.sprid:
        cref_title=[' '+args.maprefop[0]+' '+cchar for i,cchar in enumerate(args.sprid)]
    else:
        cref_title=['']*nplt

    # build final title
    ctitle=[None]*nplt
    for iplt in range(nplt):
        ctitle[iplt]=csubplt_title[iplt]+crun_title[iplt]+cref_title[iplt]

    return ctitle

def add_title(ctitle,boxxy):
    cax  = plt.axes([boxxy[0], boxxy[3], boxxy[2]-boxxy[0], 1-boxxy[3]])
    cax.text(0.5,0.5,ctitle,horizontalalignment='center',verticalalignment='bottom',fontsize=20)
    cax.axis('off')

# ================================ extra information on plot =============
def plot_section_line(fig,cfile_lst):
    for fsection in cfile_lst:
        lat,lon=get_latlon(fsection)
        fig.plot(lon.squeeze(),lat.squeeze(),'k-',linewidth=2.0,transform=ccrs.PlateCarree())

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
    delta_lon=np.abs(np.diff(lon2d))
    for i, start in enumerate(np.argmax(delta_lon > 180, axis=1)):
        lon2d[i, start+1:] += 360
    return lat2d,lon2d

def get_variable_shape(ncvar):
    redimt=re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    redimz=re.compile(r"\b(z|dep|depth|deptht|nav_lev)\b", re.I)
    redimy=re.compile(r"\b(j|y|y_grid_.+|latitude|lat|nj)\b", re.I)
    redimx=re.compile(r"\b(i|x|x_grid_.+|lon|longitude|long|ni)\b", re.I)
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
        print(dimlst)
        sys.exit(42)
    return cshape

def get_dim(cfile,cdir):

    dncdim={'x':re.compile(r"\b(x|x_grid_.+|lon|longitude|long)\b", re.I),
            'y':re.compile(r"\b(y|y_grid_.+|latitude|lat)\b", re.I),
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
def get_2d_data(cfile,cvar,ktime=0,klvl=0,offsety=None):
    print(' reading '+cvar+' in '+cfile+' ...')

    if (klvl > 0) and (ktime > 0) :
        print('error klvl or ktime larger than 0 (klvl = '+str(klvl)+', ktime = '+str(ktime)+')')
        sys.exit(42)

    if not offsety:
        nx,ny,_,_=get_dims(cfile)
        offsety=ny

    ncid   = nc.Dataset(cfile)
    clvar   = get_name(cvar,ncid.variables.keys())
    var    = ncid.variables[clvar]
    shape = get_variable_shape(var)

    dslice={
            'XY'  :(                                                  slice(0,offsety,None),slice(0,None,None) ),
            'XYT' :(slice(ktime,ktime+1,None),                        slice(0,offsety,None),slice(0,None,None) ),
            'XYZ' :(                          slice(klvl,klvl+1,None),slice(0,offsety,None),slice(0,None,None) ),
            'XYZT':(slice(ktime,ktime+1,None),slice(klvl,klvl+1,None),slice(0,offsety,None),slice(0,None,None) )
           }

    if shape=='X' :
        print(' 1d variable X => extend it 2d')
        tmp=np.zeros(shape=(ny,))
        dat2d,_=np.meshgrid(var[:],tmp)
    elif shape=='Y' :
        print(' 1d variable Y => extend it 2d')
        tmp=np.zeros(shape=(nx,))
        _,dat2d=np.meshgrid(tmp,var[:])
    elif (shape=='XY') or (shape=='XYT') or (shape=='XYZ') or (shape=='XYZT') :
        dat2d=var[dslice[shape]].squeeze()
    else:
        print(cvar+' contains '+str(len(var.shape))+' dimensions')
        print('dimension names are ',var.dimensions)
        print(' shape '+shape+' is unknown, exit ')
        sys.exit(42)

    ncid.close()

    return dat2d

def get_k_level(zlvl, cfile):
    # open netcdf
    ncid = nc.Dataset(cfile)
    cvar = get_name("(z|dep|depth|deptht|depthw)",ncid.variables.keys())
    zdep = ncid.variables[cvar][:].squeeze()

    z0=zlvl[0]
    err0=99999.
    for jk in range(0,len(zdep)):
        if np.abs(zdep[jk]-z0) < err0:
            jk0=jk
            err0=np.abs(zdep[jk]-z0)
    print('the closest level to the requiered depth is: ', jk0)

    ncid.close()

    return jk0, zdep[jk0]
