import numpy as np
from numpy import ma
from netCDF4 import Dataset
import argparse
import cartopy
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import matplotlib.colors as colors
import cartopy.crs as ccrs
import sys
from cartopy.feature import LAND
matplotlib.use('GTKAgg') 

# define argument
parser = argparse.ArgumentParser()
parser.add_argument("-f"  , metavar='file_name'     , help="names of input files"           , type=str  , nargs="+", required=True )
parser.add_argument("-v"  , metavar='var_name'      , help="variable list"                  , type=str  , nargs=1  , required=True )
parser.add_argument("-r"  , metavar='file_ref'      , help="names of ref   files"           , type=str  , nargs=1  , required=False)
parser.add_argument("-vr" , metavar='reference var_name', help="reference variable name"    , type=str  , nargs=1  , required=False)
parser.add_argument("-ft" , metavar='figure title'  , help="title of the whole figure"      , type=str  , nargs=1  , required=False)
parser.add_argument("-fid", metavar='runid'         , help="runids (title + mesh name)"     , type=str  , nargs="+", required=False)
parser.add_argument("-rid", metavar='refid'         , help="refids (title + mesh name)"     , type=str  , nargs=1  , required=False)
parser.add_argument("-c"  , metavar='color range'   , help="color range"                    , type=float, nargs=3  , required=True )
parser.add_argument("-s"  , metavar='subplot disposition' , help="subplot disposition (ixj)", type=str  , nargs=1  , required=True )
parser.add_argument("-cm" , metavar='color map name', help="color mask name"                , type=str  , nargs=1  , required=False)
parser.add_argument("-o"  , metavar='output name'   , help="output name"                    , type=str  , nargs=1  , required=False)
parser.add_argument("-p"  , metavar='projection'    , help="projection"                     , type=str  , nargs=1  , required=False)
parser.add_argument("-k"  , metavar='vertical level', help="level in fortran convention"    , type=int  , nargs=1  , required=False)
#parser.add_argument("-b"  , metavar='vertical level', help="level in fortran convention"    , type=int  , nargs=1  , required=False)
parser.add_argument("--cntf", metavar='contour file' , help="contour file list"                , type=str  , nargs="+", required=False)
parser.add_argument("--cntv", metavar='contour var ' , help="contour variable"                 , type=str  , nargs=1  , required=False)
parser.add_argument("--cntreff", metavar='contour ref file' , help="contour reference file"    , type=str  , nargs=1  , required=False)
parser.add_argument("--cntrefv", metavar='contour ref var ' , help="contour reference variable", type=str  , nargs=1  , required=False)
parser.add_argument("--cntlev", metavar='contour level', help="contour level (1 value at this stage)", type=float  , nargs=1  , required=False)
args = parser.parse_args()

# get projection and extend
if args.p:
    proj_name=args.p[0].lower()
else:
    proj_name='south_stereo'

if proj_name=='south_stereo' : 
    proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)
    latlon_lim=[-180, 180, -90, -60]
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
    joffset=-1
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
if args.k:
    jk=args.k[0]-1
else:
    jk=0

# get var
cvar_lst = args.v[:]
nvar = len(cvar_lst)

#loop over file
# for ivar in range(0,nvar):
ivar=0
cvar = cvar_lst[ivar]

# deals with ref file
if args.rid:
    ref_title=' - '+args.rid[0]
else:
    ref_title=''

if args.vr:
    cvar_ref=args.vr[ivar]
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
nplt = nfile * nvar
ax = [None] * nplt 
bc = [None] * nplt

plt.figure(figsize=np.array([210,210]) / 25.4)

# get whole figure title
if args.ft:
    fig_title=args.ft[0]
else:
    fig_title=cvar

# get color limit
if args.c:
    rmin = args.c[0] 
    rmax = args.c[1] 
    rint = args.c[2]
else:
    rmax = np.max(var2dm - ref2dm)
    rmin = np.min(var2dm - ref2dm)
    rint = 20

# get color bar
vlevel= np.arange(rmin,rmax+0.000001,rint)
if args.cm:
    cmap  = plt.get_cmap(args.cm[0],len(vlevel)-1)
else:
    cmap  = plt.get_cmap('RdBu_r',len(vlevel)-1)
if lbad:
    cmap.set_bad('0.75', 1.0)
else:
    cmap.set_under('0.75', 1.0)

# get subplot disposition
csub = args.s[0].split('x')
nisplt= int(csub[ 0])
njsplt= int(csub[-1])
if nisplt*njsplt < nfile*nvar:
    print ' number subplot lower than the number of plot asked (nfile*nvar) '
    sys.exit(1)

for ifile in range(0,nfile):

    if args.fid:
    # deals with mesh mask
        if args.fid[ifile] == 'OBS':
           ncid  = Dataset(cfile_lst[ifile])
           lat2d = ncid.variables['lat'][:]
           lon2d = ncid.variables['lon'][:]
           msk = 1
        else:
           ncid  = Dataset('mesh_mask_'+args.fid[ifile]+'.nc')
           lat2d = ncid.variables['gphit'][0,0:joffset,:]
           zlon2d = ncid.variables['glamt'][0,0:joffset,:]
           lon2d=zlon2d.copy()
           for i,start in enumerate(np.argmax(np.abs(np.diff(zlon2d)) > 180, axis=1)):
               lon2d[i, start+1:] += 360
           msk = ncid.variables['tmask'][0,jk,0:joffset,:]
        ncid.close()
    else:
        ncid  = Dataset(cfile_lst[ifile])
        lat2d = ncid.variables['nav_lat'][0:joffset,:]
        zlon2d = ncid.variables['nav_lon'][0:joffset,:]
        lon2d=zlon2d.copy()
        for i,start in enumerate(np.argmax(np.abs(np.diff(zlon2d)) > 180, axis=1)):
            lon2d[i, start+1:] += 360
        msk = 1
        ncid.close()


    # load input file
    print ' processing '+cfile_lst[ifile]+' '+cvar
    ncid   = Dataset(cfile_lst[ifile])
    var    = ncid.variables[cvar   ]
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
    ax[ifile].add_feature(cartopy.feature.LAKES,edgecolor='k',linewidth=0.5,facecolor='none')
    ax[ifile].gridlines()
    ax[ifile].set_title(ctitle_lst[ifile]+ref_title)

    # make plot
   # could be an option to not plot the map
    pcol = ax[ifile].pcolormesh(lon2d,lat2d,var2dm-ref2dm,vmin=rmin,vmax=rmax,cmap=cmap,transform=ccrs.PlateCarree(),rasterized=True)
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
        ax[ifile].contour(lon2dm,lat2dm,var2dm,levels=[cntlev, cntlev],transform=ccrs.PlateCarree(),colors='0.25',linewidth=0.5)
        if args.cntreff:
            ax[ifile].contour(lon2dm,lat2dm,cntref2dm,levels=[cntlev, cntlev],transform=ccrs.PlateCarree(),colors='gray',linewidth=0.5)

# remove extra white space
print njsplt
hpx=0.06+0.035*njsplt
plt.subplots_adjust(left=0.01,right=0.89, bottom=0.01, top=0.89, wspace=0.1, hspace=hpx)

# put common colorbar
x0=1.0; x1=0.0; y0=1.0; y1=0.0
for ifile in range(0,nfile):
    ax[ifile].apply_aspect()
    bc[ifile] = ax[ifile].get_position()
    x0=np.min([x0,bc[ifile].x0])
    x1=np.max([x1,bc[ifile].x1])
    y0=np.min([y0,bc[ifile].y0])
    y1=np.max([y1,bc[ifile].y1])

cax = plt.axes([x1+0.02, y0, 0.02, y1-y0])
cbar= plt.colorbar(pcol, cax=cax, extend='both',)
cbar.set_ticks(vlevel)
cbar.set_ticklabels(vlevel)
cbar.ax.tick_params(labelsize=14)

# put whole figure title
tax = plt.axes([x0,y1,x1-x0,1-y1])
tax.text(0.5,0.5,fig_title, 
        horizontalalignment='center',
        verticalalignment  ='center',
        fontsize=16)
tax.set_axis_off()

# get figure nname
if args.o:
    coutput_name=args.o[0]
else:
    coutput_name=outext+cvar+'_lev'+str(jk)+'_'+proj_name+'.png'

# save figure
plt.savefig(coutput_name, format='png', dpi=300)

# show figure
plt.show()

