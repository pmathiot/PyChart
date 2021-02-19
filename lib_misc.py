import re
import seawater
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
    # save argument list
    fid = open(cfile+'.txt',"w")
    fid.write(' python '+subprocess.list2cmdline(sys.argv[0:])+'\n')
    fid.close()

    # save figure
    fig.savefig(cfile+'.png', format='png', dpi=150)

# ============================ file parser =====================================
def parse_dbfile(cfile,key_lst):
    print('open file '+cfile)
    val_lst=[None]*len(key_lst)
    with open(cfile) as fid:
        for ikey,ckey in enumerate(key_lst):
            val_lst[ikey]=find_key(ckey,fid)
    # return value
    return val_lst

def find_key(char,fid):
    for cline in fid:
        lmatch = re.findall(char,cline)
        if lmatch :
            return re.split(' *= *| *',cline.strip().strip('\n'))[-1]
    return 'N/A'
# ============================ file parser end =====================================

# ============================ plot utility ========================================
def add_land_features(ax,cfeature_lst):
    dfeature={'isf':cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none'),
              'lakes':cartopy.feature.NaturalEarthFeature('physical', 'lakes'                    , '50m', facecolor='none'),
              'coast':cartopy.feature.NaturalEarthFeature('physical', 'coastline'                , '50m', facecolor='0.75'),
              'land' :cartopy.feature.NaturalEarthFeature('physical', 'land'                     , '50m', facecolor='0.75'),
              'bathy_z1000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_J_1000'  , '10m', facecolor='none'),
              'bathy_z2000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_I_2000'  , '10m', facecolor='none'),
              'bathy_z3000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_H_3000'  , '10m', facecolor='none')
             }

    for _,cfeat in enumerate(cfeature_lst):
        ax.add_feature(dfeature[cfeat],linewidth=0.5,edgecolor='k')

# ============================ CMAP ====================================
def get_lvl(bnds):
    if len(bnds)==3 :
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvlint = bnds[2]
        lvlmax = lvlmin+round((lvlmax-lvlmin)/lvlint)*lvlint
        lvl= np.arange(lvlmin,lvlmax+0.000001,lvlint)
    elif len(bnds)==2:
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvl=np.linspace(lvlmin, lvlmax, num=10)
    else:
        lvl=np.array(bnds[:])
    return lvl

def get_cmap(cpal, bnds, cnorm, cext='neither'):

    if cext not in ['neither','both','max','min']:
        print('colorbar extension should be neither, both, max or min')
        sys.exit(42)

    cmap = plt.get_cmap(cpal)

    if bnds:
        lvl=get_lvl(bnds)
    else:
        print(' Need definition of levels (min,max) at least.')
        sys.exit(42)

    if cnorm == 'BoundaryNorm':
        norm = colors.BoundaryNorm(lvl, cmap.N, extend=cext)
        #norm = colors.BoundaryNorm(lvl, cmap.N)
    elif cnorm == 'LogNorm':
        norm = colors.LogNorm(vmin=lvl[0],vmax=lvl[-1])
    elif cnorm == 'Normalize':
        norm = colors.Normalize(vmin=lvl[0],vmax=lvl[-1])

    return cmap,norm

# ============================ LEGEND ==================================
def get_corner(ax):
    x0=ax.get_position().x0
    x1=ax.get_position().x1
    y0=ax.get_position().y0
    y1=ax.get_position().y1
    return [x0,x1,y0,y1]

def add_legend(lh,ll,ncol=4,lframe=False,loc='bottom'):
    if loc=='bottom':
        lax=plt.axes([0.0, 0.0, 1.0, 0.05])
    else:
        print('legend location not yet supported, supported location are: bottom,')
    plt.legend(lh,ll,loc='center left',ncol=ncol,frameon=lframe,columnspacing=1)
    lax.set_axis_off()

# ======================= COLORBAR =======================================
def get_plt_bound(ax_lst,nplt):
# get plot corner position
    bc_lst = [None] * nplt
    x0=1.0; x1=0.0; y0=1.0; y1=0.0
    for iplt in range(0,nplt):
        ax_lst[iplt].apply_aspect()
        bc_lst[iplt] = ax_lst[iplt].get_position()
        x0=np.min([x0,bc_lst[iplt].x0])
        x1=np.max([x1,bc_lst[iplt].x1])
        y0=np.min([y0,bc_lst[iplt].y0])
        y1=np.max([y1,bc_lst[iplt].y1])
    return [x0, y0, x1, y1]

def add_colorbar(cb,boxxy,cunit='',cfmt='%5.2f',cext='neither',fontsize=16,cboffset=0.02,cbw=0.02):
    cax  = plt.axes([boxxy[2]+cboffset, boxxy[1], cbw, boxxy[3]-boxxy[1]])
    cbar = plt.colorbar(cb, cax=cax, format=cfmt,extend=cext)
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.set_title(cunit,fontsize=fontsize,y=1.0)

def get_subplt_title(args,nplt):
    # supplot title prefix
    csubplt_title=['']
    if nplt > 1:
        csubplt_title=['a) ','b) ','c) ','d) ','e) ','f) ']

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
    cax  = plt.axes([boxxy[0], boxxy[3], boxxy[1]-boxxy[0], 1-boxxy[3]])
    cax.text(0.5,0.5,ctitle,horizontalalignment='center',verticalalignment='bottom',fontsize=20)
    cax.axis('off')

# ======================= TS diag ====================================
def plot_s0_line(tmin,smin,tmax,smax,siglvl=[27.88, 27.8, 27.68, 27.55]):
# Create empty grid of zeros
    ydim=xdim=100
    dens = np.zeros((ydim,xdim))
# Create temp and salt vectors of appropiate dimensions
    ti = np.linspace(tmin,tmax,ydim)
    si = np.linspace(smin,smax,xdim)
# Loop to fill in grid with densities
    for j in range(0,int(ydim)):
        for i in range(0, int(xdim)):
            dens[j,i]=seawater.dens0(si[i],ti[j])
# Substract 1000 to convert to sigma-t
    dens = dens - 1000
# Plot data ***********************************************
    CS = plt.contour(si, ti, dens, levels=siglvl, linestyles='dashed', colors='k')
    plt.clabel(CS, fontsize=12, inline=1, fmt='%4.2f', inline_spacing=1, use_clabeltext=1) # Label every second level

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
    redimz=re.compile(r"\b(z|dep|depth|deptht)\b", re.I)
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
    ncid   = nc.Dataset(cfile)
    if cdir=='x' :
        redim=re.compile(r"\b(x|x_grid_.+|lon|longitude|long)\b", re.I)
    elif cdir=='y' :
        redim=re.compile(r"\b(y|y_grid_.+|latitude|lat)\b", re.I)
    elif cdir=='z' :
        redim=re.compile(r"\b(z|dep|depth|deptht)\b", re.I)
    elif cdir=='t' :
        redim=re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    else:
        print('dimension direction unknown, need to be x, y, z or k')
        sys.exit(42)

    cdim=list(filter(redim.match, ncid.dimensions.keys()))
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
    if not offsety:
        nx,ny,_,_=get_dims(cfile)
        offsety=ny

    ncid   = nc.Dataset(cfile)
    clvar   = get_name(cvar,ncid.variables.keys())
    var    = ncid.variables[clvar]
    shape = get_variable_shape(var)

    if shape=='X' :
        print(' 1d variable X => extend it 2d')
        tmp=np.zeros(shape=(ny,))
        dat2d,_=np.meshgrid(var[:],tmp)
    elif shape=='Y' :
        print(' 1d variable Y => extend it 2d')
        tmp=np.zeros(shape=(nx,))
        _,dat2d=np.meshgrid(tmp,var[:])
    elif shape=='XY' :
        print(' 2d variable XY')
        if (klvl > 0) and (ktime > 0) :
            print('error klvl or ktime larger than 0 (klvl = '+str(klvl)+', ktime = '+str(ktime)+')')
            sys.exit(42)
        dat2d=var[0:offsety,:]
    elif shape=='XYT' :
        print(' 3d variable XYT')
        if (klvl > 0) and (ktime > 0) :
            print('error klvl or ktime larger than 0 (klvl = '+str(klvl)+', ktime = '+str(ktime)+')')
            sys.exit(42)
        dat2d=var[ktime,0:offsety,:]
    elif shape=='XYZ' :
        print(' 3d variable XYZ')
        if ktime > 0 :
            print('error ktime larger than 0 (ktime = '+str(ktime)+')')
            sys.exit(42)
        dat2d=var[klvl,0:offsety,:]
    elif len(var.shape)==4 :
        print(' 4d variable XYZT')
        dat2d=var[ktime,klvl,0:offsety,:]
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
    return jk0, zdep[jk0]
