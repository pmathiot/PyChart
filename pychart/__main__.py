import argparse
import matplotlib.pyplot as plt
from pychart.io_utils import load_data
from pychart.plot import plot_cartesian, PlotData
from pychart.figure import FigureBuilder  # <- new class file replacing buid_figure_layout, etc.
from pychart import cb
import time


def fix_list(arg_list, nfile, name):
    """
    Ensure argument list matches nfile length.
    If one element is given, replicate it nfile times.
    If None or empty:
      - if required=True, raise error
      - else fill with default (replicated nfile times)
    """
    # Normalize None → []
    if arg_list is None:
        arg_list = [None]

    # If single value → replicate
    if len(arg_list) == 1:
        arg_list = [arg_list[0]] * nfile
    elif len(arg_list) != nfile:
        raise ValueError(
            f"Number of values for {name} ({len(arg_list)}) must be 1 or equal to number of files ({nfile})."
        )
    return arg_list

# Add a print plot summary and a print of a YML to read ?
# need to think about it

def main():
    parser = argparse.ArgumentParser(description="PyChart Command-Line Tool")
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
    parser.add_argument("--maprefjt"   , metavar='pcolor_jt_ref_file'    , help="time frame in fortran convention", \
                                      type=int  , nargs='+'  , default=[1], required=False)
    parser.add_argument("--maprefop", metavar='pcolor_ref_operation'     , help="operation made for copmarison"  , \
                                      type=str  , nargs=1  , default=[None]     , choices=[None,'-','/'], required=False)
    parser.add_argument("--maprefsf", metavar='pcolor_scale_factor'      , help="map data scale factor"          , \
                                      type=float, nargs='+', required=False)
    parser.add_argument("--mapsf"   , metavar='pcolor_scale_factor'      , help="map data scale factor"          , \
                                      type=float, nargs="+"  , default=[1.0]     , required=False)
    parser.add_argument("--mapjt"   , metavar='pcolor_jt'                , help="time frame in fortran convention", \
                                      type=int  , nargs='+'  , default=[1], required=False)

    group_map = parser.add_mutually_exclusive_group(required=False)
    group_map.add_argument("--mapz"    , metavar='pcolor_z_depth'           , help="depth of the map"               , \
                                      type=float, nargs=1  , default=[1.0], required=False)
    group_map.add_argument("--mapjk"   , metavar='pcolor_jk_depth'          , help="level in fortran convention"    , \
                                      type=int  , nargs=1  , default=[1], required=False)

    parser.add_argument("--cbn"     , metavar='colormap_name'            , help="color map name"                 , \
                                      type=str  , nargs=1  , default=['viridis']   , required=False)
    parser.add_argument('--cbcmo'   , help='use cmocean colorbar'     , action="store_true", default=False, required=False)
    parser.add_argument("--cblvl"   , metavar='colorbar_range'           , help="color range"                    , \
                                      type=float, nargs="+", required=True )
    parser.add_argument("--cbnorm"  , metavar='colorbar_norm_method'     , help="color map method (LogNorm, Normalize, BoundaryNorm, TwoSlopeNorm)", \
                                      type=str, nargs=1, default=['BoundaryNorm'], choices=['BoundaryNorm','LogNorm','Normalize','TwoSlopeNorm'], required=False )
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
    parser.add_argument("--ploc"    , metavar='gridspec indices'         , help="[0,2,0,3] : position x0,x1,y0,y1 in the plot grid", \
                                      type=str  , nargs="+",                 required=False)
    parser.add_argument("-o"        , metavar='output name'              , help="output name"                    , \
                                      type=str  , nargs=1  , default=['figure'], required=False)
    parser.add_argument("-p"        , metavar='projection'               , help="projection"                     , \
                                      type=str  , nargs="+"  , default=['ant'], required=False)
    parser.add_argument("--crs"     , metavar='sampling value'           , help="sampling value (every ncrs pts)", \
                                      type=int  , nargs=1  , default=[1],        required=False)
    parser.add_argument("--debug"   , metavar='box index [imin, imax, jmin, jmax]', help=" box index [imin, imax, jmin, jmax]", \
                                      type=int  , nargs=4  , required=False)

    parser.add_argument("--cntf"    , metavar='contour file'             , help="contour file list"              , \
                                      type=str  , nargs="+", required=False)
    parser.add_argument("--cntv"    , metavar='contour var '             , help="contour variable"               , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--cntreff" , metavar='contour ref file'         , help="contour reference file"         , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--cntrefv" , metavar='contour ref var '         , help="contour reference variable"     , \
                                      type=str  , nargs=1  , required=False)
    parser.add_argument("--cntrefop", metavar='contour_ref_operation'     , help="operation made for comparison"  , \
                                      type=str  , nargs=1  , default=[None]     , choices=[None,'-','/'], required=False)
    parser.add_argument("--cntsf"   , metavar='contour data scale factor', help="contour data scale factor"      , \
                                      type=float, nargs=1  , default=[1.0]    , required=False)
    parser.add_argument("--cntrefsf"   , metavar='contour reference data scale factor', help="contour reference data scale factor"      , \
                                      type=float, nargs=1  , default=[1.0]    , required=False)
    parser.add_argument("--cntjt"   , metavar='contour_data_jt'                , help="time frame in fortran convention", \
                                      type=int  , nargs='+'  , default=[1], required=False)

    group_cnt = parser.add_mutually_exclusive_group(required=False)
    group_cnt.add_argument("--cntjk"   , metavar='contour jk level'         , help="contour jk level "              , \
                                      type=int  , nargs=1  , required=False)
    group_cnt.add_argument("--cntz"    , metavar='contour depth in m'       , help="contour depth in m"             , \
                                      type=float, nargs=1  , required=False)

    parser.add_argument("--cntlvl"  , metavar='contour line level'       , help="contour line level"             , \
                                      type=float, nargs="+", required=False)
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
    args = parser.parse_args()

    # --- fix list consistency for map arguments ---
    nfile_map = len(args.mapf) if args.mapf else 0
    if nfile_map > 0:
        args.mapreff = fix_list(args.mapreff, nfile_map, "map reference files")
        args.mapv = fix_list(args.mapv, nfile_map, "map variables")
        args.maprefv = fix_list(args.maprefv, nfile_map, "map reference variables")
        args.maprefop = fix_list(args.maprefop, nfile_map, "map reference operation")
        args.mapsf = fix_list(args.mapsf, nfile_map, "map scale factors")
        args.maprefsf = fix_list(args.maprefsf, nfile_map, "map reference scale factors")
        args.mapjt = fix_list(args.mapjt, nfile_map, "map time frames")
        args.mapjk = fix_list(args.mapjk, nfile_map, "map levels")
        args.mapz = fix_list(args.mapz, nfile_map, "map depths")

    # --- fix list consistency for contour arguments ---
    nfile_cnt = len(args.cntf) if args.cntf else 0
    if nfile_cnt > 0:
        args.cntf = fix_list(args.cntf, nfile_cnt, "contour files")
        args.cntreff = fix_list(args.cntreff, nfile_cnt, "contour reference files")
        args.cntv = fix_list(args.cntv, nfile_cnt, "contour variables")
        args.cntrefv = fix_list(args.cntrefv, nfile_cnt, "contour reference variables")
        args.cntrefop = fix_list(args.cntrefop, nfile_map, "cnt reference operation")
        args.cntsf = fix_list(args.cntsf, nfile_cnt, "contour scale factors")
        args.cntrefsf = fix_list(args.cntrefsf, nfile_cnt, "contour reference scale factors")
        args.cntjt = fix_list(args.cntjt, nfile_cnt, "contour time frames")
        args.cntjk = fix_list(args.cntjk, nfile_cnt, "contour levels")
        args.cntz = fix_list(args.cntz, nfile_cnt, "contour depths")

    if args.ploc:
        nloc = len(args.ploc)
        if nloc != nfile_map and nloc != nfile_cnt:
            raise ValueError(
                f"Number of values for ploc ({nloc}) must be equal to number of map files ({nfile_map}) or contour files ({nfile_cnt})."
            )
        
    args.spfid = fix_list(args.spfid, nfile_map, "run name for subplot titles")
    args.sprid = fix_list(args.sprid, nfile_map, "ref name for subplot titles")

    # --- fix list consistency for fig arguments ---
    args.p = fix_list(args.p, nfile_map, "projection")

    # --- Build configuration dictionary ---
    pychart_config = {
        "figure": {
            "title": args.ft[0],
            "projection": args.p,
            "sp": args.sp[0],
            "output": args.o[0],
            "crs": args.crs[0],
            "ploc": args.ploc,
            "spfid": args.spfid,
            "sprid": args.sprid,
        },
        "map": {
            "files": args.mapf,
            "vars": args.mapv,
            "refs": args.mapreff,
            "ref_vars": args.maprefv,
            "scale": args.mapsf,
            "ref_scale": args.maprefsf,
            "jt": args.mapjt,
            "jk": args.mapjk,
            "z": args.mapz,
            "op": args.maprefop,
        },
        "cnt": {
            "files": args.cntf,
            "vars": args.cntv,
            "refs": args.cntreff,
            "ref_vars": args.cntrefv,
            "scale": args.cntsf,
            "ref_scale": args.cntrefsf,
            "op": args.cntrefop,
            "jt": args.cntjt,
            "z": args.cntz,
            "jk": args.cntjk,
            "levels": args.cntlvl,
        },
        "cb": {
            "colormap": args.cbn[0],
            "levels": args.cblvl,
            "norm": args.cbnorm[0],
            "extend": args.cbext[0],
            "units": args.cbu[0],
            "fmt": args.cbfmt[0],
            "cmocean": args.cbcmo,
        },
    }

    # --- Build the figure using the new class ---
    figure_cfg = pychart_config["figure"]

    fb = FigureBuilder(config=figure_cfg, projections_yaml="projections.yml")

    fig, axes = fb.build_layout()                 # builds figure + subplots + titles

    map_config = pychart_config["map"]
    cnt_config = pychart_config["cnt"]
    cb_config = pychart_config["cb"]

    # get map colorbar
    map_cb = cb.cb(
        cb_config["colormap"], 
        cb_config["norm"],
        cb_config["units"],
        cb_config["fmt"],
        cb_config["extend"],
        cb_config["levels"],
        cmo=cb_config["cmocean"]
        )

    print(figure_cfg["spfid"])

    # --- Loop through subplots ---
    for iax, ax in enumerate(axes):
        # MAP
        if iax < len(map_config["files"]):
            map_data = PlotData(
                file=map_config["files"][iax],
                var=map_config["vars"][iax],
                jk=int(map_config["jk"][iax]) if map_config["jk"][iax] is not None else 1,
                kt=int(map_config["jt"][iax]) if map_config["jt"][iax] is not None else 1,
                fileref=map_config["refs"][iax],
                varref=map_config["ref_vars"][iax],
                jkref=int(map_config["jk"][iax]) if map_config["jk"][iax] is not None else 1,
                ktref=int(map_config["jt"][iax]) if map_config["jt"][iax] is not None else 1,
                sf=float(map_config["scale"][iax]) if map_config["scale"][iax] is not None else 1.0,
                sfref=float(map_config["ref_scale"][iax]) if map_config["ref_scale"][iax] is not None else 1.0,
                trun=figure_cfg["spfid"][iax],
                tref=figure_cfg["sprid"][iax],
                refop=map_config["op"][iax]
            )

            if map_data.type == 'tri_unstructured':
                map_data = map_data.to_triunstructured()
            elif map_data.type == 'ico_unstructured':
                map_data = map_data.to_icounstructured()
            else:
                map_data = map_data.to_structured()

            print(map_data.trun, map_data.tref, map_data.refop)

            map_data.get_coords()
            map_data.get_data()
            map_data.compute_data()
            pcol = map_data.plot_map(ax, map_cb)  # returns QuadMesh or similar for colorbar
            map_data.add_title(ax)

        # CONTOUR
        if iax < len(cnt_config["files"]):
            cnt_data = PlotData(
                file=cnt_config["files"][iax],
                var=cnt_config["vars"][iax],
                jk=int(cnt_config["jk"][iax]) if cnt_config["jk"][iax] is not None else 1,
                kt=int(cnt_config["jt"][iax]) if cnt_config["jt"][iax] is not None else 1,
                fileref=cnt_config["refs"][iax],
                varref=cnt_config["ref_vars"][iax],
                jkref=int(cnt_config["jk"][iax]) if cnt_config["jk"][iax] is not None else 1,
                ktref=int(cnt_config["jt"][iax]) if cnt_config["jt"][iax] is not None else 1,
                sf=float(cnt_config["scale"][iax]) if cnt_config["scale"][iax] is not None else 1.0,
                sfref=float(cnt_config["ref_scale"][iax]) if cnt_config["ref_scale"][iax] is not None else 1.0,
                refop=cnt_config["op"][iax] 
            )

            if cnt_data.type == 'tri_unstructured':
                cnt_data = cnt_data.to_triunstructured()
            elif cnt_data.type == 'ico_unstructured':
                cnt_data = cnt_data.to_icounstructured()
            else:
                cnt_data = cnt_data.to_structured()
            
            cnt_data.get_coords()
            cnt_data.get_data()
            cnt_data.compute_data()
            cntlvl=cb.get_lvl(args.cntlvl)
            cnt_data.plot_cnt(ax, levels=cntlvl, colors='k', linewidths=1)

        else:
            ax.set_visible(False)

    # --- Add colorbar using the class method ---
    fb.add_colorbar(pcol,map_cb)

    print(f"Saving figure to {figure_cfg['output']}")
    # make layout tight
    #plt.tight_layout()
    plt.show()
    fig.savefig(figure_cfg["output"], dpi=150)


if __name__ == "__main__":
    main()