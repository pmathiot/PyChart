import re
import cartopy.crs as ccrs
import netCDF4 as nc
import sys
sys.path.insert(0, '/home/h05/pmathiot/PYTHON/MISC/seawater/')
import seawater
import cartopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# ============================ output argument list in txt file ================
def output_argument_lst(cfile, arglst):
    fid = open(cfile,"w")
    fid.write(' python2.7 '+' '.join(arglst))
    fid.close()

# ============================ file parser =====================================
def parse_dbfile(cfile,key_lst):
    print 'open file '+cfile
    val_lst=[None]*len(key_lst)
    with open(cfile) as fid:
       for ikey,ckey in enumerate(key_lst):
           val_lst[ikey]=find_key(ckey,fid)
    # return value
    return val_lst

def find_key(char,fid):
    for cline in fid:
        lmatch = re.findall(char,cline)
        if (lmatch) :
            return re.split(' *= *| *',cline.strip().strip('\n'))[-1]
    return 'N/A'
# ============================ file parser end =====================================

# ============================ plot utility ========================================
def add_land_features(ax,cfeature_lst):
# get isf groiunding line, ice shelf front and coastline
    for ifeat,cfeat in enumerate(cfeature_lst):
        if cfeat=='isf':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', '50m',facecolor='none',edgecolor='k')
        elif cfeat=='lakes':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'lakes'                      , '50m',facecolor='none',edgecolor='k')
        elif cfeat=='coast':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'coastline'                  , '50m',facecolor='0.75',edgecolor='k')
        elif cfeat=='land':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'land'                       , '50m',facecolor='0.75',edgecolor='k')
        elif cfeat=='bathy_z1000':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_J_1000'          , '10m',facecolor='none',edgecolor='k')
        elif cfeat=='bathy_z2000':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_I_2000'          , '10m',facecolor='none',edgecolor='k')
        elif cfeat=='bathy_z3000':
            feature = cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_H_3000'          , '10m',facecolor='none',edgecolor='k')
        else:
            print 'feature unknown : '+cfeat
            sys.exit(42)
        ax.add_feature(feature,linewidth=0.5)

# ============================ CMAP ====================================
def get_lvl(bnds):
    if len(bnds)==3 :
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvlint = bnds[2]
        lvlmin = round(lvlmin/lvlint)*lvlint
        lvlmax = round(lvlmax/lvlint)*lvlint
        lvl= np.arange(lvlmin,lvlmax+0.000001,lvlint)
    elif len(bnds)==2:
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvl=np.linspace(lvlmin, lvlmax, num=20)
    else:
        lvlmin = bnds[0]
        lvlmax = bnds[-1]
        lvl=bnds[:]
    return lvl

def get_cmap(cpal, bnds, cext='neither', cbad='w'):
    if bnds:
        lvl=get_lvl(bnds)
        lvlmin=lvl[0] ; lvlmax=lvl[-1]
    else:
        print ' Need definition of levels (min,max) at least.'
        sys.exit(42)

    nintlvl=len(lvl)-1
    if cext=='neither':
        ntotlvl=len(lvl)-1 ; imin=0 ; imax=ntotlvl
    elif cext=='both':
        ntotlvl=len(lvl)+1 ; imin=1 ; imax=ntotlvl-1
    elif cext=='max':
        ntotlvl=len(lvl)   ; imin=0 ; imax=ntotlvl-1
    elif cext=='min':
        ntotlvl=len(lvl)   ; imin=1 ; imax=ntotlvl
    else:
        print 'colorbar extension should be neither, both, max or min'
        sys.exit(42)
    cmap = plt.get_cmap(cpal,ntotlvl)
    cmap = cmap(np.arange(ntotlvl)) ; cunder=cmap[0]; cover=cmap[-1]
    cmap = cmap[imin:imax]
    cmap=colors.LinearSegmentedColormap.from_list("cmap", cmap, nintlvl)
    cmap.set_bad(cbad, 1.0)
    cmap.set_under(cunder)
    cmap.set_over(cover)
    return cmap,lvl,lvlmin,lvlmax

# ============================ LEGEND ==================================
def get_corner(ax):
    x0=ax.get_position().x0
    x1=ax.get_position().x1
    y0=ax.get_position().y0
    y1=ax.get_position().y1
    return x0,x1,y0,y1

def add_legend(lh,ll,ncol=4,lframe=False,loc='bottom'):
   if loc=='bottom':
       lax=plt.axes([0.0, 0.0, 1.0, 0.05])
   else:
       print 'legend location not yet supported, supported location are: bottom,'
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
    return x0, y0, x1, y1

def add_colorbar(plt,cb,x0,y0,x1,y1,lvl=None,cunit='',cfmt='%5.2f',cext='neither',fontsize=14,cboffset=0.02,cbw=0.02):
    cax  = plt.axes([x1+cboffset, y0, cbw, y1-y0])
    cbar = plt.colorbar(cb, cax=cax, format=cfmt, extend=cext)
    cbar.ax.tick_params(labelsize=fontsize)
    cbar.ax.set_title(cunit)
    if lvl is not None:
        cbar.set_ticks(lvl)

def add_title(plt,ctitle,x0,y0,x1,y1):
    cax  = plt.axes([x0, y1, x1-x0, 1-y1])
    cax.text(0.5,0.5,ctitle,horizontalalignment='center',verticalalignment='center',fontsize=16)
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
def plot_section_line(plt,cfile_lst):
    for fsection in cfile_lst:
        print fsection
        lat,lon=get_latlon(fsection)
        plt.plot(lon.squeeze(),lat.squeeze(),'k-',linewidth=2.0,transform=ccrs.PlateCarree())

# ================================= NETCDF ===============================
def get_name(regex,varlst):
    revar = re.compile(r'\b%s\b'%regex,re.I)
    cvar  = filter(revar.match, varlst)
    if (len(cvar) > 1):
        print regex+' name list is longer than 1 or 0; error'
        print cvar
        print cvar[0]+' is selected'
    if (len(cvar) == 0):
        print 'no match between '+regex+' and :'
        print varlst
        sys.exit(42)
    return cvar[0]

def get_latlon_var(cfile):
    ncid   = nc.Dataset(cfile)
    clon=get_name("(nav_lon.*|lon|longitude|glamt)",ncid.variables.keys())
    clat=get_name("(nav_lat.*|lat|latitude|gphit)",ncid.variables.keys())
    ncid.close()
    return clat,clon

def get_latlon(cfile,offsety=None):
    clat,clon=get_latlon_var(cfile)
    lat2d =get_2d_data(cfile,clat,offsety=offsety)
    lon2d=get_2d_data(cfile,clon,offsety=offsety)
    delta_lon=np.abs(np.diff(lon2d))
    j_lst,i_lst=np.nonzero(delta_lon>180)
    for jj in j_lst: 
        if i_lst != [] :
            ii=i_lst[jj]
            lon2d[jj, ii+1:] += 360
    return lat2d,lon2d

def get_variable_shape(ncid,ncvar):
    redimt=re.compile(r"\b(t|tim|time_counter|time)\b", re.I)
    redimz=re.compile(r"\b(z|dep|depth|deptht)\b", re.I)
    redimy=re.compile(r"\b(y|y_grid_.+|latitude|lat)\b", re.I)
    redimx=re.compile(r"\b(x|x_grid_.+|lon|longitude|long)\b", re.I)
    dimlst = ncvar.dimensions
    if (len(ncvar.shape)==2) and redimx.match(dimlst[1]) and redimy.match(dimlst[0]):
        cshape='XY'
    elif (len(ncvar.shape)==3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape='XYT'
    elif (len(ncvar.shape)==3) and redimx.match(dimlst[2]) and redimy.match(dimlst[1]) and redimz.match(dimlst[0]):
        cshape='XYZ'
    elif (len(ncvar.shape)==4) and redimx.match(dimlst[3]) and redimy.match(dimlst[2]) and redimz.match(dimlst[1]) and redimt.match(dimlst[0]):
        cshape='XYZT'
    else:
        print 'cshape undefined, error'
        print dimlst
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
        print 'dimension direction unknown, need to be x, y, z or k'
        sys.exit(42)

    cdim=filter(redim.match, ncid.dimensions.keys());
    if (len(cdim) > 1):
        print regex+' name list is longer than 1; error'
        print cdim
        sys.exit(42)
    elif (len(cdim) == 0):
        print cdir+' dim in '+cfile+' is 0.'
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
    print ' reading '+cvar+' in '+cfile+' ...'
    if not offsety:
        nx,ny,nz,nt=get_dims(cfile)
        offsety=ny

    ncid   = nc.Dataset(cfile)
    var    = ncid.variables[cvar]
    shape = get_variable_shape(ncid,var)

    if shape=='XY' :
        print ' 2d variable XY'
        if (klvl > 0) :
            print 'error klvl larger than 0 (klvl = '+str(klvl)+')'
            sys.exit(42)
        if (ktime > 0) :
            print 'error ktime larger than 0 (ktime = '+str(ktime)+')'
            sys.exit(42)
        dat2d=var[0:offsety,:]
    elif shape=='XYT' :
        print ' 3d variable XYT'
        if (klvl > 0) :
            print 'error klvl larger than 0 (klvl = '+str(klvl)+')'
            sys.exit(42)
        if (ktime > 0) :
            print 'error ktime larger than 0 (ktime = '+str(ktime)+')'
            sys.exit(42)
        dat2d=var[ktime,0:offsety,:]
    elif shape=='XYZ' :
        print ' 3d variable XYZ'
        if (ktime > 0) :
            print 'error ktime larger than 0 (ktime = '+str(ktime)+')'
            sys.exit(42)
        dat2d=var[klvl,0:offsety,:]
    elif len(var.shape)==4 :
        print ' 4d variable XYZT'
        dat2d=var[ktime,klvl,0:offsety,:]
    else:
        print cvar+' contains '+str(len(var.shape))+' dimensions'
        print 'dimension names are '+var.dimensions
        print ' shape unknown, exit '
        sys.exit(42)
    ncid.close()
    return dat2d

def get_k_level(zlvl, cfile):
    # open netcdf
    ncid = nc.Dataset(cfile)
    cvar = get_name("(z|dep|depth|deptht)",ncid.variables.keys())
    zdep = ncid.variables[cvar][:].squeeze()

    z0=zlvl[0]
    err0=99999.
    for jk in range(0,len(zdep)):
        if np.abs(zdep[jk]-z0) < err0:
            jk0=jk
            err0=np.abs(zdep[jk]-z0)
    print 'the closest level to the requiered depth is: '
    print jk0,zdep[jk0]
    return jk0

