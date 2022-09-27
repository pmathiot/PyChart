#!/usr/bin/python
import sys
import argparse

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

import lib_misc as libpc

import cmocean as cmo

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

    if args.mapjk and args.mapz:
        print('ERROR, --mapjk and --mapz defined together, check script argument list')
        sys.exit(42)

def get_jk(jk0=None,z0=None,cfile=None):
    jk=0
    if cfile:
        if jk0:
            jk=jk0[0]
        elif z0:
            jk,_=libpc.get_k_level(z0,cfile)
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
    parser.add_argument("--dir"     , metavar='data_dir'                 , help="data dir"                       , \
                                      type=str  , nargs=1  , default=['./']    , required=False)

    parser.add_argument("--mapf"    , metavar='pcolor_file_names'        , help="names of input files"           , \
                                      type=str  , nargs="+", required=True )
    parser.add_argument("--mapv"    , metavar='pcolor_var_names'         , help="variable list"                  , \
                                      type=str  , nargs="+", required=True )
    parser.add_argument("--mapreff" , metavar='pcolor_ref_file_name'     , help="names of ref   files"           , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--maprefv" , metavar='pcolor_ref_var_name'      , help="reference variable name"        , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--maprefop", metavar='pcolor_ref_operation'     , help="operation made for copmarison"  , \
                                      type=str  , nargs=1  , default=['-']     , choices=['-','/'], required=False)
    parser.add_argument("--mapsf"   , metavar='pcolor_scale_factor'      , help="map data scale factor"          , \
                                      type=float, nargs="+"  , default=[1.0]     , required=False)
    parser.add_argument("--mapjk"   , metavar='pcolor_jk_depth'          , help="level in fortran convention"    , \
                                      type=int  , nargs=1  , required=False)
    parser.add_argument("--mapz"    , metavar='pcolor_z_depth'           , help="depth of the map"               , \
                                      type=float, nargs=1  , required=False)

    parser.add_argument("--cbn"     , metavar='colormap_name'            , help="color map name"                 , \
                                      type=str  , nargs=1  , default=['viridis']   , required=False)
    parser.add_argument('--cbcmo'   , help='use cmocean colorbar'     , action="store_true", default=False, required=False)
    parser.add_argument("--cblvl"   , metavar='colorbar_range'           , help="color range"                    , \
                                      type=float, nargs="+", required=True )
    parser.add_argument("--cbnorm"  , metavar='colorbar_norm_method'     , help="color map method (LogNorm, Normalize, BoundaryNorm)", \
                                      type=str, nargs=1, default=['BoundaryNorm'], choices=['BoundaryNorm','LogNorm','Normalize'], required=False )
    parser.add_argument("--cbu"     , metavar='colorbar_unit'            , help="colorbar unit"                  , \
                                      type=str  , nargs=1  , default=['']      , required=False)
    parser.add_argument("--cbfmt"   , metavar='colorbar_fmt'             , help="colorbar format"                , \
                                      type=str  , nargs=1  , default=['%5.2f'] , required=False)
    parser.add_argument("--cbext"   , metavar='colorbar_extend'          , help="colorbar extend"                , \
                                      type=str  , nargs=1  , default=['both'], choices=['both', 'neither', 'max', 'min']  , required=False)

    parser.add_argument("--ft"      , metavar='figure_title'             , help="title of the whole figure"      , \
                                      type=str  , nargs=1  , default=['']      , required=False)
    parser.add_argument("--spfid"   , metavar='runid'                    , help="runids (title + mesh name)"     , \
                                      type=str  , nargs="+", default=['']      , required=False)
    parser.add_argument("--sprid"   , metavar='refid'                    , help="refids (title + mesh name)"     , \
                                      type=str  , nargs="+",                     required=False)
    parser.add_argument("--mask"    , metavar='mask file name'           , help="mask file name"                 , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--mesh"    , metavar='mesh file name'           , help="mesh file name"                 , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--llonce"  , metavar='read lat/lon each time'   , help="read lat/lon for each plot [0=no]", \
                                      type=int, nargs=1  , default=[0]   , choices=[0, 1]    , required=False)
    parser.add_argument("--sp"      , metavar='subplot disposition'      , help="subplot disposition (ixj)"      , \
                                      type=str  , nargs=1  , default=['1x1']   , required=False)
    parser.add_argument("--ploc"    , metavar='gridspec indices'         , help="0,0 : top left plot, 0,: : top line", \
                                      type=str  , nargs="+",                 required=False)
    parser.add_argument("-o"        , metavar='output name'              , help="output name"                    , \
                                      type=str  , nargs=1  , default=['figure'], required=False)
    parser.add_argument("-p"        , metavar='projection'               , help="projection"                     , \
                                      type=str  , nargs=1  , default=['global'], required=False)
    parser.add_argument("--crs"     , metavar='sampling value'           , help="sampling value (every ncrs pts)", \
                                      type=int  , nargs=1  , default=[1],        required=False)

    parser.add_argument("--cntf"    , metavar='contour file'             , help="contour file list"              , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--cntv"    , metavar='contour var '             , help="contour variable"               , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--cntreff" , metavar='contour ref file'         , help="contour reference file"         , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--cntrefv" , metavar='contour ref var '         , help="contour reference variable"     , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--cntsf"   , metavar='contour data scale factor', help="contour data scale factor"      , \
                                      type=float, nargs=1  , default=[1.0]    , required=False)
    parser.add_argument("--cntjk"   , metavar='contour jk level'         , help="contour jk level "              , \
                                      type=int  , nargs=1  , required=False)
    parser.add_argument("--cntz"    , metavar='contour depth in m'       , help="contour depth in m"             , \
                                      type=float, nargs=1  , required=False)
    parser.add_argument("--cntlvl"  , metavar='contour line level'       , help="contour line level"             , \
                                      type=float, nargs="+", required=False)
    parser.add_argument("--cntrefop", metavar='contour_ref_operation'   , help="operation made for copmarison"  , \
                                      type=str  , nargs=1  , default=['-']     , choices=['-','/'], required=False)
    parser.add_argument("--bathyf"  , metavar='bathy file'               , help="bathy file"                     , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--bathyv"  , metavar='bathy var '               , help="contour variable"               , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--bathylvl", metavar='contour line level'       , help="contour line level"             , \
                                      type=float, nargs="+", required=False)
    parser.add_argument("--secf"    , metavar='section line file list '  , help="section file list describing section to plot", \
                                      type=str, nargs="+", required=False)
    parser.add_argument("--joffset" , metavar='offset on j'              , help="do not read the top j lines, it could be needed for some grid (ORCA like for example) and some projection", \
                                      type=int  , nargs=1  , default=[0],required=False)
    return parser.parse_args()

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
                            [(-180,  180,  60,  90),ccrs.PlateCarree()] ],
           'ross'        :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=-180.0), \
                            [(-6.67e5,8.33e5,1.05e6,2.47e6), 'cproj' ] ],
           'amundsen'    :[ ccrs.Stereographic(central_latitude=-90.0, central_longitude=0.0), \
                            [( -99, -130, -70, -78),ccrs.PlateCarree()] ],
           'tip_ant_pen' :[ ccrs.Stereographic(central_latitude= -90.0, central_longitude=0.0)   , \
                            [ (-3.407e6,-1.69e6,1.012e6,2.71e6),'cproj' ] ],
                            #[ (-3.307e6,-1.99e6,1.112e6,2.41e6),'cproj' ] ],
           'global'         :[ ccrs.PlateCarree()                    , ['global'] ],
           'global_robinson':[ ccrs.Robinson(central_longitude=0)    , ['global'] ],
           'global_mercator':[ ccrs.Mercator(central_longitude=-90.0), ['global'] ],
           'natl'      : [ ccrs.LambertConformal(-40, 45,cutoff=20), [(-4.039e6,2.192e6,-1.429e6,4.805e6), 'cproj' ] ],
           'greenland' : [ ccrs.LambertConformal(-40, 45,cutoff=20), [(-1.124e6,0.897e6,1.648e6,5.198e6), 'cproj' ] ]
          }
#    elif proj_name=='natl' :
#        proj=ccrs.LambertConformal(-40, 45,cutoff=20)
#        XY_lim=[-4.039e6,2.192e6,-1.429e6,4.805e6]
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
# =======================================================================================================================================================

def main():

    # get argument list
    args=get_argument()

    sanity_check(args)

    # get projection and extend
    proj, XY_lim = def_projection(args.p[0])

    # deals with ref file
    mapref2d=0.0
    joffset=args.joffset[0]

    # get file and title list and sanity check
    # to be simplify as cnt and map use same code (use of class a class ?)
    if args.mapf:
        cmaprunfile,cmaprunvar=get_file_and_varname(args.dir[0],args.mapf[:],args.mapv[:])
        mapjk=get_jk(args.mapjk,args.mapz,cmaprunfile[0])

        if args.mapreff:
            cmapreffile,cmaprefvar=get_file_and_varname(args.dir[0],args.mapreff[:],args.maprefv[:])

    if args.cntf:
        # get file
        ccntrunfile,ccntrunvar=get_file_and_varname(args.dir[0], args.cntf[:],args.cntv[:])
        cntjk=get_jk(args.cntjk,args.cntz,ccntrunfile[0])

        if args.cntreff:
            ccntreffile,ccntrefvar=get_file_and_varname(args.dir[0],args.cntreff[:],args.cntrefv[:])

    # define contour lvl
    if args.cntf:
        cntlvl=libpc.get_lvl(args.cntlvl)
        cntclr=[None]*len(cntlvl)
        for ii, val in enumerate(cntlvl):
            if val<0:
                cntclr[ii]='0.75'
            if val >=0:
                cntclr[ii]='k'

    # get map colorbar
    cextend=args.cbext[0]
    cmap, norm = libpc.get_cmap(args.cbn[0],args.cblvl,args.cbnorm[0], cext=cextend, cmo=args.cbcmo)

    # get map/cnt scale factor
    map_sf=args.mapsf
    cnt_sf=args.cntsf[0]

    # get subplot disposition
    if args.mapf:
        nplt=len(cmaprunfile)
    elif args.cntf:
        nplt=len(ccntrunfile)
    nisplt,njsplt = get_subplot(args.sp[0],nplt)

    # title list
    csptitle = libpc.get_subplt_title(args,nplt)

    # initialisation
    ax = [None] * nplt
    lpltloc=[None]*nplt

    # define figure dimension
    fig=plt.figure(figsize=np.array([297,297*njsplt/nisplt]) / 25.4)

    # define grid specification
    gs = fig.add_gridspec(njsplt, nisplt)

    for ifile in range(0,nplt):
        print('')
        print('start plotting subplot {}/{}'.format(ifile+1,nplt))
        print('')
        # define subplot
        if args.ploc :
            lpltloc[ifile]=eval('gs['+args.ploc[ifile]+']')
        else:
            lpltloc[ifile] = gs[ifile//nisplt,ifile%nisplt]

        # deal with mesh
        if args.mesh:
            cmeshf=get_var_lst(args.mesh,args.mapf)
            cllfile=cmeshf[ifile]
        else:
            cllfile=cmaprunfile[ifile]

        print(ifile,args.llonce[0])
        if (ifile == 0):
            lat2d,lon2d=libpc.get_latlon(cllfile,joffset)
        if (args.llonce[0]==0 and ifile > 0):
            print('toto')
            lat2d,lon2d=libpc.get_latlon(cllfile,joffset)

        # define subplot
        ax[ifile] = fig.add_subplot(lpltloc[ifile], projection=proj)

        # put proj, extend, grid ...
        if XY_lim[0] == 'global':
            ax[ifile].set_global()
        else:
            ax[ifile].set_extent(XY_lim[0], XY_lim[1])

        libpc.add_land_features(ax[ifile],['isf','lakes','land'])
        ax[ifile].gridlines(linewidth=1, color='k', linestyle='--')
        ax[ifile].set_title(csptitle[ifile],fontsize=18)

        # make plot
        ncrs=args.crs[0]

        # add map if ask
        if args.mapf:

            # deals with mask
            msk = 1.0
            if args.mask:
                cmsk=args.mask[ifile]
                print('open '+cmsk)
                mapjk=get_jk(args.mapjk,args.mapz,cmaprunfile[ifile])
                msk = libpc.get_2d_data(cmsk,'tmask',klvl=mapjk,offsety=joffset)
                msk = np.ma.masked_where(msk==0.0,msk)

            print('plot pcolormesh ...')
            mapjk=get_jk(args.mapjk,args.mapz,cmaprunfile[ifile])
            mapvar2d  = libpc.get_2d_data(cmaprunfile[ifile],cmaprunvar[ifile],klvl=mapjk,offsety=joffset)
            mapvar2dm = np.ma.masked_where(mapvar2d*msk==0.0,mapvar2d)
            if args.mapreff:
                # def file list
                mapref2d=libpc.get_2d_data(cmapreffile[ifile],cmaprefvar[ifile],klvl=mapjk,offsety=joffset)
                mapref2dm = np.ma.masked_where(msk*mapref2d==0.0,mapref2d)
            else:
                mapref2dm = 0.0

            if args.maprefop[0] == '-':
                maptoplot2d=(mapvar2dm-mapref2dm)*map_sf[ifile]
            elif args.maprefop[0] == '/':
                maptoplot2d=(mapvar2dm/mapref2dm)
            pcol = ax[ifile].pcolormesh(lon2d[::ncrs,::ncrs],lat2d[::ncrs,::ncrs],maptoplot2d[::ncrs,::ncrs], \
                                        cmap=cmap,norm=norm,transform=ccrs.PlateCarree(),rasterized=True)

        # add contour if ask
        if args.cntf:
            print('plot contour ...')
            # need to be simplify as same code as map (input runf,runv,reff,refv,jk,offset,msk,cnt_sf)
            cntvar2d  = libpc.get_2d_data(ccntrunfile[ifile],ccntrunvar[ifile],klvl=cntjk,offsety=joffset)
            cntvar2dm = np.ma.masked_where(cntvar2d==0.0,cntvar2d)
            if args.cntreff:
                cntref2d=libpc.get_2d_data(ccntreffile[ifile],ccntrefvar[ifile],klvl=cntjk,offsety=joffset)
                cntref2dm = np.ma.masked_where(msk*cntref2d==0.0,cntref2d)
            else:
                cntref2dm = 0.0

            if args.cntrefop[0] == '-':
                cnttoplot2d=(cntvar2dm-cntref2dm)*cnt_sf
            elif args.cntrefop[0] == '/':
                cnttoplot2d=(cntvar2dm/cntref2dm)

            ax[ifile].contour(lon2d[::ncrs,::ncrs],lat2d[::ncrs,::ncrs],cnttoplot2d[::ncrs,::ncrs], \
                              levels=cntlvl,transform=ccrs.PlateCarree(),colors=cntclr,linewidths=1)

        # add bathy line if ask
        if args.bathyf:
            cbathyf=get_var_lst(args.bathyf,args.mapf)
            cbathyv=get_var_lst(args.bathyv,args.mapf)
            bathy2d=libpc.get_2d_data(cbathyf[ifile],cbathyv[ifile],offsety=joffset)
            bathy2dm = np.ma.masked_where(bathy2d==0.0,bathy2d)

            print('plot bathymetry ...')
            ax[ifile].contour(lon2d[::ncrs,::ncrs],lat2d[::ncrs,::ncrs],bathy2dm[::ncrs,::ncrs], \
                              levels=args.bathylvl[:],transform=ccrs.PlateCarree(),colors='0.5',linewidths=0.5)

        # add section line if ask
        if args.secf:
            print('plot section line')
            libpc.plot_section_line(ax[ifile],args.secf[:])

    print('')

    # remove extra white space
    hpx=0.06+0.035*njsplt
    fig.subplots_adjust(left=0.01,right=0.88, bottom=0.02, top=0.89, wspace=0.1, hspace=hpx)

    # get_figure_corner position
    corner_coord = libpc.get_plt_bound(ax, nplt) # left, bottom, right, top

    # add common colorbar
    if args.mapf:
        print('add colorbar')
        libpc.add_colorbar(pcol,corner_coord,cunit=args.cbu[0],cfmt=args.cbfmt[0],cext=cextend)

    # put whole figure title
    libpc.add_title(args.ft[0],corner_coord)

    # argument lst output
    print('save figure')
    libpc.save_output(args.o[0],fig)

    # show figure
    plt.show()

if __name__ == '__main__':
    main()
