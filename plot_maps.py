#!/usr/local/sci/bin/python2.7

import matplotlib
#matplotlib.use('GTKAgg') 
import numpy as np
from numpy import ma
import os
import argparse
import cartopy
import matplotlib.pyplot as plt
import matplotlib as matplotlib
import matplotlib.colors as colors
from matplotlib.colors import Normalize
import cartopy.crs as ccrs
import sys
import re
from cartopy.feature import LAND
sys.path.insert(0,'/home/h05/pmathiot/PYTHON/PYCHART_GIT/')
from lib_misc import *
#matplotlib.use('GTKAgg') 

def def_output_name(arg_output,cext,cvar,jk,cproj):
    if arg_output:
        coutput=arg_output
    else:
        coutput=cext+cvar+'_lev'+str(jk)+'_'+cproj
    return coutput

def sanity_check(args):
    # sanity check
    if args.spfid:
        if len(args.spfid) != len(args.mapf):
            print('title list and file list not the same length, exit')
            sys.exit(2)

    if args.cntf:
        if len(args.cntf) != len(args.spfid):
            print('Length of contour file list and title list is not the same length, exit')
            sys.exit(42)

def get_jk(jk0=None,z0=None,cfile=None):
    jk=0
    if cfile:
        if jk0 and z0:
            print('ERROR, jk and z define, check script argument list')
            sys.exit(42)
        if jk0:
            jk=jk0[0]
        elif z0:
            jk=get_k_level(z0,cfile)
    return jk

def get_var_lst(cvar,cfile):
    cvar_lst=cvar
    nvar=len(cvar_lst)
    nfile=len(cfile)
    if len(cvar) == 1:
        print('Length of variable list is 1, we assume it is the same variable for all input file')
        cvar_lst=[cvar[0]]*nfile
    elif nvar != nfile:
        print('length of variable list should be 1 or '+str(nfile))
        sys.exit(42)
    return cvar_lst

# get file and title list and sanity check
def get_file_and_varname(cdir,cfile_lst,cvar_lst):
# def file list
    cpath_lst  = [cdir + f for f in cfile_lst]
# get map var list 
    cvar_lst_ext=get_var_lst(cvar_lst,cpath_lst)
    return cpath_lst,cvar_lst_ext

# get subplot disposition
def get_subplot(csubplt,nplt):
    csub = csubplt.split('x')
    ni=int(csub[ 0])
    nj=int(csub[-1])
    if ni*nj < nplt:
        print(ni*nj,' panels')
        print(nplt, ' plot asked')
        print(' number subplot lower than the number of plot asked ')
        print(' add -s XxY option in the command line')
        sys.exit(42)
    return ni,nj


def get_argument():
# define argument
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir"     , metavar='data_dir'                 , help="data dir"                       , type=str  , nargs=1  , default=['./']    , required=False)
    parser.add_argument("--mapf"    , metavar='pcolor_file_names'        , help="names of input files"           , type=str  , nargs="+", required=True )
    parser.add_argument("--mapv"    , metavar='pcolor_var_names'         , help="variable list"                  , type=str  , nargs="+", required=True )
    parser.add_argument("--mapreff" , metavar='pcolor_ref_file_name'     , help="names of ref   files"           , type=str  , nargs=1  , required=False)
    parser.add_argument("--maprefv" , metavar='pcolor_ref_var_name'      , help="reference variable name"        , type=str  , nargs=1  , required=False)
    parser.add_argument("--mapsf"   , metavar='pcolor_scale_factor'      , help="map data scale factor"          , type=float, nargs=1  , default=[1.0]     , required=False)
    parser.add_argument("--mapjk"   , metavar='pcolor_jk_depth'          , help="level in fortran convention"    , type=int  , nargs=1  , required=False)
    parser.add_argument("--mapz"    , metavar='pcolor_z_depth'           , help="depth of the map"               , type=float, nargs=1  , required=False)
    parser.add_argument("--cbn"     , metavar='colormap_name'            , help="color map name"                 , type=str  , nargs=1  , default=['jet']   , required=False)
    parser.add_argument("--cblvl"   , metavar='colorbar_range'           , help="color range"                    , type=float, nargs="+", required=True )
    parser.add_argument("--cbu"     , metavar='colorbar_unit'            , help="colorbar unit"                  , type=str  , nargs=1  , default=['']      , required=False)
    parser.add_argument("--cbfmt"   , metavar='colorbar_fmt'             , help="colorbar format"                , type=str  , nargs=1  , default=['%5.2f'] , required=False)
    parser.add_argument("--cbext"   , metavar='colorbar_extend'          , help="colorbar extend"                , type=str  , nargs=1  , default=['both']  , required=False)
    parser.add_argument("--ft"      , metavar='figure_title'             , help="title of the whole figure"      , type=str  , nargs=1  , default=['']      , required=False)
    parser.add_argument("--spfid"   , metavar='runid'                    , help="runids (title + mesh name)"     , type=str  , nargs="+", default=['']      , required=False)
    parser.add_argument("--sprid"   , metavar='refid'                    , help="refids (title + mesh name)"     , type=str  , nargs=1  , default=['']      , required=False)
    parser.add_argument("--mask"    , metavar='mask file name'           , help="mask file name"                 , type=str  , nargs="+", required=False)
    parser.add_argument("--mesh"    , metavar='mesh file name'           , help="mesh file name"                 , type=str  , nargs="+", required=False)
    parser.add_argument("--sp"      , metavar='subplot disposition'      , help="subplot disposition (ixj)"      , type=str  , nargs=1  , default=['1x1']   , required=False)
    parser.add_argument("-o"        , metavar='output name'              , help="output name"                    , type=str  , nargs=1  , default=['output'], required=False)
    parser.add_argument("-p"        , metavar='projection'               , help="projection"                     , type=str  , nargs=1  , default=['global'], required=False)
    parser.add_argument("--cntf"    , metavar='contour file'             , help="contour file list"              , type=str  , nargs="+", required=False)
    parser.add_argument("--cntv"    , metavar='contour var '             , help="contour variable"               , type=str  , nargs=1  , required=False)
    parser.add_argument("--cntreff" , metavar='contour ref file'         , help="contour reference file"         , type=str  , nargs=1  , required=False)
    parser.add_argument("--cntrefv" , metavar='contour ref var '         , help="contour reference variable"     , type=str  , nargs=1  , required=False)
    parser.add_argument("--cntsf"   , metavar='contour data scale factor', help="contour data scale factor"      , type=float, nargs=1  , default=[1.0]    , required=False)
    parser.add_argument("--cntjk"   , metavar='contour jk level'         , help="contour jk level "              , type=int  , nargs=1  , required=False)
    parser.add_argument("--cntz"    , metavar='contour depth in m'       , help="contour depth in m"             , type=float, nargs=1  , required=False)
    parser.add_argument("--cntlvl"  , metavar='contour line level'       , help="contour line level"             , type=float, nargs="+", required=False)
    parser.add_argument("--bathyf"  , metavar='bathy file'               , help="bathy file"                     , type=str  , nargs="+", required=False)
    parser.add_argument("--bathyv"  , metavar='bathy var '               , help="contour variable"               , type=str  , nargs=1  , required=False)
    parser.add_argument("--bathylvl", metavar='contour line level'       , help="contour line level"             , type=float, nargs="+", required=False)
    parser.add_argument("--secf"    , metavar='section line file list '  , help="section file list describing section to plot", type=str, nargs="+", required=False)
    return parser.parse_args()

def def_projection(proj_name):
    if proj_name=='ortho_natl' :
        proj=ccrs.Orthographic(central_longitude=-60.0, central_latitude=45.0)
        joffset=-1
        XY_lim=[-180, 180, -90, -60]
    elif proj_name=='south_stereo' : 
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)
        XY_lim=[-180, 180, -90, -60]
        joffset=-2
    elif proj_name=='south_ocean' : 
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)
        XY_lim=[(-180, 180, -90, -45),ccrs.PlateCarree()]
        joffset=-2
    elif proj_name=='ant' : 
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0)
        XY_lim=[-180, 180, -90, -65]
        joffset=-2
    elif proj_name=='arctic' : 
        proj=ccrs.Stereographic(central_latitude=90.0, central_longitude=0.0)
        XY_lim=[-180, 180, 60, 90]
        joffset=-2
    elif proj_name=='global' :
        proj=ccrs.PlateCarree()
        XY_lim=['global']
        joffset=-1
    elif proj_name=='natl'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-4.039e6,2.192e6,-1.429e6,4.805e6]
        joffset=-2
    elif proj_name=='greenland'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-1.124e6,0.897e6,1.648e6,5.198e6]
        joffset=-1
    elif proj_name=='ovf'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-3.553e5,2.141e6,9.915e5,3.4113e6]
        joffset=-2
    elif proj_name=='ovf_larger'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-2.5e5,2.45e6,0.9e6,3.6e6]
        joffset=-2
    elif proj_name=='irminger'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[-2.164e5,1.395e6,1.635e6,3.265e6]
        joffset=-2
    elif proj_name=='japan'         :
        proj=ccrs.LambertConformal(150, 30,cutoff=10)
        XY_lim=[-3.166e6,2.707e6,-1.008e6,4.865e6]
        joffset=-1
    elif proj_name=='global_robinson':
        proj=ccrs.Robinson()
        XY_lim=['global']
    elif proj_name=='global_mercator':
        proj=ccrs.Mercator(central_longitude=-90.0)
        XY_lim=['global']
        joffset=-1
    elif proj_name=='feroe'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[1.56e6,2.145e6,1.973e6,2.555e6]
        joffset=-1
    elif proj_name=='gulf_stream'         :
        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
        XY_lim=[(-4.250e6,1.115e5,-1.546e6,2.8155e6), proj]
        joffset=-1
    elif proj_name=='ross':
        proj=ccrs.Stereographic(central_latitude=-90.0, central_longitude=-180.0)
        XY_lim=[(-6.67e5,8.33e5,1.05e6,2.47e6), proj]
        joffset=-2
    else:
        print('projection '+proj_name+' unknown')
        print('should be ross, gulf_stream, feroe, global_mercator, global_robinson, japan, ovf, greenland, natl, global, south_stereo, ant')
        sys.exit(42)
    return proj, XY_lim, joffset
# =======================================================================================================================================================

def main():

    # get argument list
    args=get_argument()
    
    # get file and title list and sanity check
    if args.mapf:
        cmaprunfile,cmaprunvar=get_file_and_varname(args.dir[0],args.mapf[:],args.mapv[:])
        mapjk=get_jk(args.mapjk,args.mapz,cmaprunfile[0])
    
    sanity_check(args)
   
    # get title list for each subplot 
    if args.spfid:
        ctitle_lst = args.spfid[:]
        
    # get projection and extend
    proj, XY_lim, joffset = def_projection(args.p[0])
    
    # deals with ref file
    creft=''
    mapref2d=0.0
    if args.mapreff:
        # def file list
        if args.maprefv:
            cmapreffile,cmaprefvar=get_file_and_varname(args.dir[0], args.mapreff[:],args.maprefv[:])
        else:
            cmapreffile = args.mapreff[:]
            cmaprefvar  = cmaprunvar[:]
        nmapreffile = len(cmapreffile)
    
        if args.sprid:
            creft=' - '+args.sprid[0]
        print(nmapreffile,' ref',cmapreffile[0],args.mapreff )
        if nmapreffile==1 :
            mapref2d=get_2d_data(cmapreffile[0],cmaprefvar[0],klvl=mapjk,offsety=joffset)
        else:
            print('more than 1 ref file, not yet implememnnted')
            sys.exit(42)
    
    if args.cntf:
        # get file and title list and sanity check
        ccntrunfile,ccntrunvar=get_file_and_varname(args.dir[0], args.cntf[:],args.cntv[:])
        cntjk=get_jk(args.cntjk,args.cntz,ccntrunfile[0])
        cntlvl=get_lvl(args.cntlvl)
        cntclr=[None]*len(cntlvl)
        for ii, val in enumerate(cntlvl):
            if val<0:
                cntclr[ii]='0.75'
            if val >=0:
                cntclr[ii]='k'
    
        if args.cntreff:
            cntref2d=get_2d_data(args.cntreff[0],args.cntrefv[0],klvl=cntjk,offsety=joffset)
        else:
            cntref2dm = 0.0
    
    
    
    # get whole figure title
    fig_title=cmaprunvar[0]
    if args.ft:
        fig_title=args.ft[0]
    
    # get map colorbar
    cextend=args.cbext[0]
    cmap, norm, maplvl, rmin, rmax = get_cmap(args.cbn[0],args.cblvl,cext=cextend)
    
    # get map/cnt scale factor
    map_sf=args.mapsf[0]
    cnt_sf=args.cntsf[0]
    
    # get subplot disposition
    if args.mapf:
        nplt=len(cmaprunfile)
    elif args.cntf:
        nplt=len(ccntrunfile)
    nisplt,njsplt = get_subplot(args.sp[0],nplt)
    
    # supplot title prefix
    csubplt_title=['']
    if nplt > 1:
        csubplt_title=['a) ','b) ','c) ','d) ','e) ','f) ']
    
    # initialisation
    ax = [None] * nplt
    
    # define figure dimension
    plt.figure(figsize=np.array([297,297*njsplt/nisplt]) / 25.4)
    #lt.figure(figsize=np.array([210,210*njsplt/nisplt]) / 25.4)
    
    for ifile in range(0,nplt):
    
        # load input file
        mapvar2d=get_2d_data(cmaprunfile[ifile],cmaprunvar[ifile],klvl=mapjk,offsety=joffset)

        # initialisation msk
        msk = 1.0
        if args.mask:
        # deals with mesh mask 
            cmsk=args.mask[ifile]
            if os.path.isfile(cmsk):
               print('open '+cmsk)
               msk = get_2d_data(cmsk,'tmask',klvl=mapjk,offsety=joffset)
               msk = ma.masked_where(msk==0.0,msk)
    
        if args.mesh:
            cmeshf=get_var_lst(args.mesh,args.mapf)
            lat2d,lon2d=get_latlon(cmeshf[ifile],joffset) 
        else:
            lat2d,lon2d=get_latlon(cmaprunfile[ifile],joffset) 
    
        mapvar2d = ma.masked_where(mapvar2d*msk==0.0,mapvar2d)
    
        # mask data
        if args.mapreff:
            mapref2d = ma.masked_array(mapref2d, mapvar2d.mask)
    
        # define subplot
        ax[ifile] = plt.subplot(njsplt, nisplt, ifile+1, projection=proj)
    
        # put proj, extend, grid ...
        if XY_lim[0] == 'global':
             ax[ifile].set_global()
        else:
            ax[ifile].set_extent(XY_lim[0], XY_lim[1])
    
        add_land_features(ax[ifile],['isf','lakes','land'])
        ax[ifile].gridlines(linewidth=1, color='k', linestyle='--')
        ax[ifile].set_title(csubplt_title[ifile]+ctitle_lst[ifile]+creft)
    
        # make plot
        ncrs=1
        if args.mapf:
        # could be an option to not plot the map
            print('plot pcolormesh ...')
            maptoplot2d=(mapvar2d-mapref2d)*map_sf
            pcol = ax[ifile].pcolormesh(lon2d[::ncrs,::ncrs],lat2d[::ncrs,::ncrs],maptoplot2d[::ncrs,::ncrs],cmap=cmap,norm=norm,vmin=rmin,vmax=rmax,transform=ccrs.PlateCarree(),rasterized=True)
    
        # add contour if ask
        if args.cntf:
            var2d=get_2d_data(ccntrunfile[ifile],ccntrunvar[ifile],klvl=cntjk,offsety=joffset)
            var2dm = ma.masked_where(var2d==0.0,var2d)
            if args.cntreff:
                cntref2dm = ma.masked_where(msk*cntref2d==0.0,cntref2d)
    
            print('plot contour ...')
            cnttoplot=(var2dm-cntref2dm)*cnt_sf
            ax[ifile].contour(lon2d[::ncrs,::ncrs],lat2d[::ncrs,::ncrs],cnttoplot[::ncrs,::ncrs],levels=cntlvl,transform=ccrs.PlateCarree(),colors=cntclr,linewidths=1)
            
        # add bathy line if ask
        if args.bathyf:
            cbathyf=get_var_lst(args.bathyf,args.mapf)
            cbathyv=get_var_lst(args.bathyv,args.mapf)
            bathy2d=get_2d_data(cbathyf[ifile],cbathyv[ifile],offsety=joffset)
            bathy2dm = ma.masked_where(bathy2d==0.0,bathy2d)
            print('plot bathymetry ...')
            ax[ifile].contour(lon2d,lat2d,bathy2dm,levels=args.bathylvl[:],transform=ccrs.PlateCarree(),colors='0.5',linewidths=0.5)
    
        # add section line if ask
        if args.secf:
            print('plot section line')
            plot_section_line(ax[ifile],args.secf[:])
    
    # remove extra white space
    hpx=0.06+0.035*njsplt
    plt.subplots_adjust(left=0.01,right=0.88, bottom=0.01, top=0.89, wspace=0.1, hspace=hpx)
    
    # get_figure_corner position
    xl, yb, xr, yt = get_plt_bound(ax, nplt) # left, bottom, right, top
    
    # add common colorbar
    if args.mapf:
        add_colorbar(plt,pcol,xl,yb,xr,yt,lvl=maplvl[:],cunit=args.cbu[0],cfmt=args.cbfmt[0],cext=cextend)
    
    # put whole figure title
    add_title(plt,fig_title,xl,yb,xr,yt)
    
    # get figure name
    if args.o:
        coutput_name=args.o[0]
    else:
        coutput_name=def_output_name(args,outxt,cvar,mapjk,args.p[0])
    
    # argument lst output
    output_argument_lst(coutput_name+'.txt',sys.argv)
    
    # save figure
    plt.savefig(coutput_name+'.png', format='png', dpi=150)
    
    # show figure
    print('show figure ...')
    plt.show()

if __name__ == '__main__':
    main()

