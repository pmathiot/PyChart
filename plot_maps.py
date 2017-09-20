import numpy as np
from numpy import ma
from netCDF4 import Dataset
import argparse
import cartopy
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import matplotlib.colors as colors
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import sys
from cartopy.feature import LAND
matplotlib.use('GTKAgg') 

def output_argument_lst(cfile, arglst):
    fid = open(cfile,"w") 
    fid.write(' python2.7'+' '.join(arglst))
    fid.close()

def def_output_name(arg_output,cext,cvar,jk,cproj):
    if arg_output:
        coutput=arg_output[0]
    else:
        coutput=cext+cvar+'_lev'+str(jk)+'_'+cproj
    return coutput

def def_cmap_lvl(bnds):
# get color limit
    if bnds:
        if len(bnds)==3 :
            lvlmin = bnds[0]
            lvlmax = bnds[1]
            lvlint = bnds[2]
            lvl= np.arange(lvlmin,lvlmax+0.000001,lvlint)
        else:
            lvl=bnds[:]
    else:
        print ' Need definition of levels (min,max,int) at least.'
        sys.exit(42)

    return lvl


def def_cmap(cm,lvl):
# get color bar
    nlvl=len(lvl)
    if cm:
        cmap  = plt.get_cmap(cm[0],nlvl-1)
    else:
        cmap  = plt.get_cmap('RdBu_r',nlvl-1)

    if lbad:
        cmap.set_bad('0.75', 1.0)
    else:
        cmap.set_under('0.75', 1.0)

    norm = colors.BoundaryNorm(boundaries=lvl, ncolors=nlvl-1)

    return cmap,norm

def get_land_features():
# get isf groiunding line, ice shelf front and coastline
    isf_features   = cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', '50m',facecolor='none',edgecolor='k')
    coast_features = cartopy.feature.NaturalEarthFeature('physical', 'coastline'                  , '50m',facecolor='0.75',edgecolor='k')

    return isf_features, coast_features

def get_plt_bound(ax_lst,nplt):
# get plot corner position
    bc_lst = [None] * nplt
    x0=1.0; x1=0.0; y0=1.0; y1=0.0
    for iplt in range(0,nplt):
        ax_lst[iplt].apply_aspect()
        bc_lst[iplt] = ax[iplt].get_position()
        x0=np.min([x0,bc_lst[iplt].x0])
        x1=np.max([x1,bc_lst[iplt].x1])
        y0=np.min([y0,bc_lst[iplt].y0])
        y1=np.max([y1,bc_lst[iplt].y1])

    return x0, y0, x1, y1

def draw_colorbar(plt,cb,lvl,x0,y0,x1,y1):
# draw colorbar
    cax  = plt.axes([x1+0.02, y0, 0.02, y1-y0])
    cbar = plt.colorbar(cb, cax=cax, format='%4.2f', extend='both',)
    cbar.set_ticks(lvl)
    cbar.ax.tick_params(labelsize=14)

def write_figure_title(ctitle,x0,y0,x1,y1):
# put whole figure title
    tax = plt.axes([x0,y1,x1-x0,1-y1])
    tax.text(0.5,0.5,ctitle,
            horizontalalignment='center',
            verticalalignment  ='center',
            fontsize=16)
    tax.set_axis_off()

def get_k_level(klvl, zlvl, cfile, cvar):
    if klvl:
        jk0=klvl[0]
    elif zlvl:
    # open netcdf
        ncid = Dataset(cfile)
        zdep = ncid.variables[cvar][:].squeeze()
 
        z0=zlvl[0]
        err0=99999.
        for jk in range(0,len(zdep)):
            if np.abs(zdep[jk]-z0) < err0:
                jk0=jk
                err0=np.abs(zdep[jk]-z0)
        print 'the closest level to the requiered depth is: '
        print jk0,zdep[jk0]
    else:
        jk0=0

    return jk0

           
# add write of the text file for the option

def get_argument():
# define argument
    parser = argparse.ArgumentParser()
    parser.add_argument("-f"  , metavar='file_name'     , help="names of input files"           , type=str  , nargs="+", required=True )
    parser.add_argument("-v"  , metavar='var_name'      , help="variable list"                  , type=str  , nargs="+", required=True )
    parser.add_argument("-r"  , metavar='file_ref'      , help="names of ref   files"           , type=str  , nargs=1  , required=False)
    parser.add_argument("-vr" , metavar='reference var_name', help="reference variable name"    , type=str  , nargs=1  , required=False)
    parser.add_argument("-ft" , metavar='figure title'  , help="title of the whole figure"      , type=str  , nargs=1  , required=False)
    parser.add_argument("-fid", metavar='runid'         , help="runids (title + mesh name)"     , type=str  , nargs="+", required=False)
    parser.add_argument("-rid", metavar='refid'         , help="refids (title + mesh name)"     , type=str  , nargs=1  , required=False)
    parser.add_argument("-c"  , metavar='color range'   , help="color range"                    , type=float, nargs="+", required=True )
    parser.add_argument("-s"  , metavar='subplot disposition' , help="subplot disposition (ixj)", type=str  , nargs=1  , required=True )
    parser.add_argument("-cm" , metavar='color map name', help="color mask name"                , type=str  , nargs=1  , required=False)
    parser.add_argument("-o"  , metavar='output name'   , help="output name"                    , type=str  , nargs=1  , required=False)
    parser.add_argument("-p"  , metavar='projection'    , help="projection"                     , type=str  , nargs=1  , required=False)
    parser.add_argument("-k"  , metavar='vertical level', help="level in fortran convention"    , type=int  , nargs=1  , required=False)
    parser.add_argument("-z"  , metavar='depth of the map', help="depth of the map"             , type=float, nargs=1  , required=False)
    parser.add_argument("--cntf", metavar='contour file' , help="contour file list"                , type=str  , nargs="+", required=False)
    parser.add_argument("--cntv", metavar='contour var ' , help="contour variable"                 , type=str  , nargs=1  , required=False)
    parser.add_argument("--cntreff", metavar='contour ref file' , help="contour reference file"    , type=str  , nargs=1  , required=False)
    parser.add_argument("--cntrefv", metavar='contour ref var ' , help="contour reference variable", type=str  , nargs=1  , required=False)
    parser.add_argument("--cntlev", metavar='contour level', help="contour level (1 value at this stage)", type=float  , nargs=1  , required=False)
    return parser.parse_args()

def def_projection(arg_proj):
    if arg_proj:
        proj_name=arg_proj[0].lower()
    else:
        proj_name='south_stereo'

    if proj_name=='south_stereo' : 
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)
        latlon_lim=[-180, 180, -90, -60]
        lbad=True
        joffset=-2
    if proj_name=='ant' : 
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)
        latlon_lim=[-180, 180, -90, -65]
        lbad=True
        joffset=-2
    if proj_name=='global' :
        proj=ccrs.PlateCarree()
        global_lim='T'
        lbad=False
        joffset=-1
    if proj_name=='natl'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-4.250e6,3.694e6,-1.546e6,6.398e6]
        lbad=True
        joffset=-2
    if proj_name=='greenland'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-1.124e6,0.897e6,1.648e6,5.198e6]
        lbad=True
        joffset=-1
    if proj_name=='ovf'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-3.553e5,2.141e6,9.915e5,3.4113e6]
        lbad=True
        joffset=-1
    if proj_name=='japan'         :
        proj=ccrs.LambertConformal(150, 30,cutoff=10)
        XY_lim=[-3.166e6,2.707e6,-1.008e6,4.865e6]
        lbad=True
        joffset=-1
    if proj_name=='global_robinson':
        proj=ccrs.Robinson()
        global_lim='T'
        lbad=False
        joffset=-1
    if proj_name=='global_mercator':
        proj=ccrs.Mercator(central_longitude=-90.0)
        global_lim='T'
        lbad=False
        joffset=-1
    if proj_name=='feroe'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[1.56e6,2.145e6,1.973e6,2.555e6]
        lbad=True
        joffset=-1
    if proj_name=='gulf_stream'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-4.250e6,1.115e5,-1.546e6,2.8155e6]
        lbad=True
        joffset=-1
    if proj_name=='ross':
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=-180.0)
        XY_lim=[-6.67e5,8.33e5,1.05e6,2.47e6]
        lbad=True
        joffset=-2
    return proj, XY_lim, lbad, joffset
# =======================================================================================================================================================
# get argument list
args=get_argument()

# get projection and extend
proj, XY_lim, lbad, joffset = def_projection(args.p)

# get file and title list and sanity check
cfile_lst  = args.f[:]
nfile = len(cfile_lst)
if args.fid:
    ctitle_lst = args.fid[:]
    if len(cfile_lst) != len(ctitle_lst):
        print 'title list and file list not the same length, exit'
        sys.exit(2)
else:
    ctitle_lst = ['']

# get k level
jk=get_k_level(args.k,args.z,args.f[0],'deptht')

# get var
cvar_lst = args.v[:]
nvar = len(cvar_lst)
cvar = cvar_lst[0]
if nvar == 1:
    cvar_lst=[cvar]*nfile

# deals with ref file
if args.rid:
    ref_title=' - '+args.rid[0]
else:
    ref_title=''

if args.vr:
    cvar_ref=args.vr[0]
else:
    cvar_ref=cvar

if args.r:
    ref_file    = args.r[0]
    ncid   = Dataset(ref_file)
    var    = ncid.variables[cvar_ref]
    if len(var.shape)==2 : 
        print ' 2d variable XY'
        ref2d=var[0:joffset,:] 
    elif len(var.shape)==3 :
        print ' 3d variable XYT'
        ref2d=var[0,0:joffset,:]
    elif len(var.shape)==4 :
        print ' 4d variable XYZT'
        ref2d=var[0,jk,0:joffset,:]
    else:
        print var.shape
        print cvar+' contains '+str(len(var.shape))+' dimensions' 
        print ' shape unknown, exit '
    outext='diff_'
    ncid.close()
else:
    outext=''
    ref2dm = 0.0

if args.cntreff:
    ref_file    = args.cntreff[0]
    ncid   = Dataset(ref_file)
    var    = ncid.variables[args.cntrefv[0]]
    if len(var.shape)==2 :
        print ' 2d variable XY'
        cntref2d=var[0:joffset,:]
    elif len(var.shape)==3 :
        print ' 3d variable XYT'
        cntref2d=var[0,0:joffset,:]
    elif len(var.shape)==4 :
        print ' 4d variable XYZT'
        cntref2d=var[0,jk,0:joffset,:]
    else:
        print var.shape
        print cvar+' contains '+str(len(var.shape))+' dimensions'
        print ' shape unknown, exit '
    outext='diff_'
    ncid.close()
else:
    outext=''
    cntref2dm = 0.0


# initialisation
nplt = nfile
ax = [None] * nplt 

# get whole figure title
if args.ft:
    fig_title=args.ft[0]
else:
    fig_title=cvar

# get color limit
vlevel = def_cmap_lvl(args.c)

# get color bar
cmap, norm = def_cmap(args.cm,vlevel)

# load land feature
isf_features, coast_features = get_land_features()

# get subplot disposition
csub = args.s[0].split('x')
nisplt= int(csub[ 0])
njsplt= int(csub[-1])
if nisplt*njsplt < nfile:
    print nisplt*njsplt,' panels'
    print nfile*nvar, ' plot asked'
    print ' number subplot lower than the number of plot asked (nfile*nvar) '
    sys.exit(1)

# define figure dimension
plt.figure(figsize=np.array([210,210*njsplt/nisplt]) / 25.4)

for ifile in range(0,nfile):

    if args.fid:
    # deals with mesh mask
        if args.fid[ifile] == 'OBS':
           ncid_msh = Dataset(cfile_lst[ifile])
           lat2d = ncid_msh.variables['lat'][:]
           lon2d = ncid_msh.variables['lon'][:]
           msk = 1
        else:
           print 'mesh_hgr_'+args.fid[ifile]+'.nc'
           ncid_msh = Dataset('mesh_hgr_'+args.fid[ifile]+'.nc')
           lat2d  = ncid_msh.variables['gphit'][0:joffset,:]
           zlon2d = ncid_msh.variables['glamt'][0:joffset,:]
           lon2d=zlon2d.copy()
           for i,start in enumerate(np.argmax(np.abs(np.diff(zlon2d)) > 180, axis=1)):
               lon2d[i, start+1:] += 360
           
           print 'mesh_mask_'+args.fid[ifile]+'.nc'
           ncid_msk = Dataset('mesh_mask_'+args.fid[ifile]+'.nc')
           msk = ncid_msk.variables['tmask'][0,jk,0:joffset,:]
        ncid_msk.close()
    else:
        ncid_msh  = Dataset(cfile_lst[ifile])
        lat2d = ncid_msh.variables['nav_lat'][0:joffset,:]
        zlon2d = ncid_msh.variables['nav_lon'][0:joffset,:]
        lon2d=zlon2d.copy()
        for i,start in enumerate(np.argmax(np.abs(np.diff(zlon2d)) > 180, axis=1)):
            lon2d[i, start+1:] += 360
        msk = 1
    ncid_msh.close()


    # load input file
    print ' processing '+cfile_lst[ifile]+' '+cvar_lst[ifile]
    ncid   = Dataset(cfile_lst[ifile])
    var    = ncid.variables[cvar_lst[ifile]]
    if len(var.shape)==2 : 
        print ' 2d variable XY'
        var2d=var[0:joffset,:] 
    elif len(var.shape)==3 :
        print ' 3d variable XYT'
        var2d=var[0,0:joffset,:]
    elif len(var.shape)==4 :
        print ' 4d variable XYZT'
        var2d=var[0,jk,0:joffset,:]
    else:
        print var.shape
        print cvar_lst[ifile]+' contains '+str(len(var.shape))+' dimensions' 
        print ' shape unknown, exit '
        sys.exit(1)
    var2dm = ma.masked_where(msk*var2d==0.0,var2d)
    if args.r:
        ref2dm = ma.masked_where(msk*ref2d==0.0,ref2d)
    if (not lbad) :
        var2dm = var2dm.filled(-99)
    ncid.close()

    # define subplot
    ax[ifile] = plt.subplot(njsplt, nisplt, ifile+1, projection=proj, axisbg='0.75')

    # put proj, extend, grid ...
    if 'latlon_lim' in globals():
        ax[ifile].set_extent(latlon_lim, ccrs.PlateCarree())
    elif 'XY_lim' in globals():
        ax[ifile].set_xlim(XY_lim[0:2]); ax[ifile].set_ylim(XY_lim[2:4]) 
    elif global_lim:
        ax[ifile].set_global()
    else:
        print ' plot limit unknown, exit'
        sys.exit(1)

    ax[ifile].coastlines(resolution='50m',linewidth=0.5)
    #ax[ifile].add_feature(coast_features,linewidth=0.5) # issue with natl proj (why ????)
    ax[ifile].add_feature(isf_features,linewidth=0.5)
    ax[ifile].add_feature(cartopy.feature.LAKES,edgecolor='k',linewidth=0.5,facecolor='none')
    ax[ifile].gridlines()
    ax[ifile].set_title(ctitle_lst[ifile]+ref_title)

    # make plot
   # could be an option to not plot the map
    pcol = ax[ifile].pcolormesh(lon2d,lat2d,var2dm-ref2dm,cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),rasterized=True)

    # could be an option to add a contour over the map
    if args.cntf:
        cntfile_lst=args.cntf[:]
        cntvar=args.cntv[0]
        cntlev=args.cntlev[0]
        print ' processing '+cntfile_lst[ifile]+' '+cntvar
        ncid   = Dataset(cntfile_lst[ifile])
        var    = ncid.variables[cntvar   ]
        if len(var.shape)==2 :
            print ' 2d variable XY'
            var2d=var[0:joffset,:]
        elif len(var.shape)==3 :
            print ' 3d variable XYT'
            var2d=var[0,0:joffset,:]
        elif len(var.shape)==4 :
            print ' 4d variable XYZT'
            var2d=var[0,jk,0:joffset,:]
        else:
            print var.shape
            print cvar+' contains '+str(len(var.shape))+' dimensions'
            print ' shape unknown, exit '
            sys.exit(1)
        var2dm = ma.masked_where(var2d==0.0,var2d)
        lon2dm = ma.masked_where(var2d==0.0,lon2d)
        lat2dm = ma.masked_where(var2d==0.0,lat2d)
        if args.cntreff:
            cntref2dm = ma.masked_where(msk*cntref2d==0.0,cntref2d)
        if (not lbad) :
            var2dm = var2dm.filled(-99)
        ncid.close()
        ax[ifile].contour(lon2dm,lat2dm,var2dm,levels=[cntlev, cntlev],transform=ccrs.PlateCarree(),colors='0.25',linewidths=0.5)
        if args.cntreff:
            ax[ifile].contour(lon2dm,lat2dm,cntref2dm,levels=[cntlev, cntlev],transform=ccrs.PlateCarree(),colors='gray',linewidths=0.5)

# remove extra white space
hpx=0.06+0.035*njsplt
plt.subplots_adjust(left=0.01,right=0.89, bottom=0.01, top=0.89, wspace=0.1, hspace=hpx)

# put common colorbar
xl, yb, xr, yt = get_plt_bound(ax, nfile) # left, bottom, right, top

# add_colorbar
draw_colorbar(plt,pcol,vlevel,xl,yb,xr,yt)

# put whole figure title
write_figure_title(fig_title,xl,yb,xr,yt)

# get figure name
coutput_name=def_output_name(args.o[0],outext,cvar,jk,proj_name)

# argument lst output
output_argument_lst(coutput_name+'.txt',sys.argv)

# save figure
plt.savefig(coutput_name+'.png', format='png', dpi=150)

# show figure
plt.show()
