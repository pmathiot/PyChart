import argparse
import matplotlib.pyplot as plt
from .plot import add_map_plot, add_cnt_plot
from .figure import FigureBuilder
from .log import set_debug, info, debug
from .cb import cb

"""
cli.py

This module provides the command-line interface (CLI) for the PyChart application. It includes functions for parsing arguments, performing sanity checks, and building the configuration for generating plots.

Functions:
- parse_args(): Parses command-line arguments and returns them as a namespace.
- fix_list(arg_list, nfile, name): Ensures argument lists match the required length.
- sanity_checks_args(args): Performs sanity checks and fixes inconsistencies in the parsed arguments.
- get_config(args): Builds the configuration dictionary from the parsed arguments.
- main(): The main entry point for the CLI, orchestrating the figure generation process.
"""

def parse_args():
    """
    Parse command-line arguments for the PyChart application.

    Returns:
    - argparse.Namespace: Parsed arguments as a namespace.
    """
    parser = argparse.ArgumentParser(
        description=(
            "PyChart Command-Line Tool\n\n"
            "Minimal required arguments for a basic call:\n"
            "  --mapf FILES --mapv VARS\n"
            "Other arguments are optional and grouped by function below."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    # ------------------------- Required Arguments -------------------------
    required_group = parser.add_argument_group("Required arguments for minimal call")
    required_group.add_argument(
        "--mapf", metavar='FILES', type=str, nargs='+', required=True,
        help="Input data file(s) for the main variable"
    )
    required_group.add_argument(
        "--mapv", metavar='VARS', type=str, nargs='+', required=True,
        help="Variable name(s) in the input files"
    )

    # ------------------------- Map Options -------------------------
    map_group = parser.add_argument_group("Map options")
    map_group.add_argument("--dir", metavar='DATA_DIR', type=str, nargs=1, default=['./'],
                           help="Directory containing data files")
    map_group.add_argument("--mapreff", metavar='REF_FILES', type=str, nargs='+',
                           help="Reference file(s) for comparison")
    map_group.add_argument("--maprefv", metavar='REF_VARS', type=str, nargs='+',
                           help="Reference variable(s)")
    map_group.add_argument("--maprefop", metavar='REF_OP', type=str, nargs=1, default=[None],
                           choices=[None,'-','/'], help="Operation for comparison with reference")
    map_group.add_argument("--mapsf", metavar='SCALE', type=float, nargs='+', default=[1.0],
                           help="Scale factor for main data")
    map_group.add_argument("--maprefsf", metavar='REF_SCALE', type=float, nargs='+', default=[1.0],
                           help="Scale factor for reference data")
    map_group.add_argument("--mapjt", metavar='TIME', type=int, nargs='+', default=[1],
                           help="Time frame(s) in Fortran convention")
    map_group.add_argument("--maprefjt", metavar='REF_TIME', type=int, nargs='+', default=[1],
                           help="Reference time frame(s) in Fortran convention")

    # Mutually exclusive: depth or level
    map_depth_group = parser.add_mutually_exclusive_group()
    map_depth_group.add_argument("--mapz", metavar='DEPTH', type=float, nargs=1, default=[1.0],
                                 help="Depth of the map (in meters)")
    map_depth_group.add_argument("--mapjk", metavar='LEVEL', type=int, nargs=1, default=[1],
                                 help="Level in Fortran convention")

    # ------------------------- Colormap Options -------------------------
    cb_group = parser.add_argument_group("Colormap options")
    cb_group.add_argument("--cbn", metavar='COLORMAP', type=str, nargs=1, default=['viridis'],
                          help="Colormap name (matplotlib)")
    cb_group.add_argument("--cbcmo", action="store_true", default=False,
                          help="Use cmocean colormap instead of matplotlib")
    cb_group.add_argument("--cblvl", metavar='LEVELS', type=float, nargs='+', default=[None],
                          help="Colorbar levels or range")
    cb_group.add_argument("--cbnorm", metavar='NORM', type=str, nargs=1,
                          choices=['BoundaryNorm','LogNorm','Normalize','TwoSlopeNorm'], default=['BoundaryNorm'],
                          help="Normalization method for color mapping")
    cb_group.add_argument("--cbu", metavar='UNIT', type=str, nargs=1, default=[''],
                          help="Colorbar unit")
    cb_group.add_argument("--cbfmt", metavar='FORMAT', type=str, nargs=1, default=['%5.2f'],
                          help="Colorbar format")
    cb_group.add_argument("--cbext", metavar='EXTEND', type=str, nargs=1,
                          choices=['both', 'neither', 'max', 'min'], default=['both'],
                          help="Colorbar extend")

    # ------------------------- Contour Options -------------------------
    cnt_group = parser.add_argument_group("Contour options")
    cnt_group.add_argument("--cntf", metavar='FILES', type=str, nargs='+',
                           help="Contour data files")
    cnt_group.add_argument("--cntv", metavar='VARS', type=str, nargs=1,
                           help="Contour variable")
    cnt_group.add_argument("--cntreff", metavar='REF_FILES', type=str, nargs=1,
                           help="Reference contour file")
    cnt_group.add_argument("--cntrefv", metavar='REF_VARS', type=str, nargs=1,
                           help="Reference contour variable")
    cnt_group.add_argument("--cntrefop", metavar='REF_OP', type=str, nargs=1, default=[None],
                           choices=[None,'-','/'], help="Operation for contour comparison")
    cnt_group.add_argument("--cntsf", metavar='SCALE', type=float, nargs=1, default=[1.0],
                           help="Contour data scale factor")
    cnt_group.add_argument("--cntrefsf", metavar='REF_SCALE', type=float, nargs=1, default=[1.0],
                           help="Contour reference scale factor")
    cnt_group.add_argument("--cntjt", metavar='TIME', type=int, nargs='+', default=[1],
                           help="Contour time frames")

    # Mutually exclusive: contour depth or level
    cnt_depth_group = parser.add_mutually_exclusive_group()
    cnt_depth_group.add_argument("--cntz", metavar='DEPTH', type=float, nargs=1,
                                 help="Contour depth (meters)")
    cnt_depth_group.add_argument("--cntjk", metavar='LEVEL', type=int, nargs=1,
                                 help="Contour depth level")

    cnt_group.add_argument("--cntlvl", metavar='LEVELS', type=float, nargs='+',
                           help="Contour line levels")

    # ------------------------- Figure Layout Options -------------------------
    fig_group = parser.add_argument_group("Figure layout options")
    fig_group.add_argument("--ft", metavar='TITLE', type=str, nargs=1, default=[''],
                           help="Figure title")
    fig_group.add_argument("--spfid", metavar='RUNID', type=str, nargs='+', default=[''],
                            help="Run ID(s) for main data (used in figure title)")
    fig_group.add_argument("--sprid", metavar='REFID', type=str, nargs='+', default=[''],
                            help="Reference run ID(s) for comparison (used in figure title)")
    fig_group.add_argument("--sp", metavar='SUBPLOTS', type=str, nargs=1, default=['1x1'],
                           help="Subplot layout (ixj)")
    fig_group.add_argument("--ploc", metavar='GRID', type=str, nargs='+',
                           help="Custom subplot location in gridspec: [x0,x1,y0,y1]")
    fig_group.add_argument("-o", metavar='OUTPUT', type=str, nargs=1, default=['figure'],
                           help="Output file name")
    fig_group.add_argument("-p", metavar='PROJ', type=str, nargs='+', default="noproj",
                            help="Projection name(s)")
    fig_group.add_argument("--projbox", metavar='BOX', type=int, nargs=4, default=None,
                            help="Box index for subregion: imin imax jmin jmax")
    fig_group.add_argument("--crs", metavar='N', type=int, nargs=1, default=[1],
                           help="Sampling value (every n-th point)")
    fig_group.add_argument("--joffset", metavar='OFFSET', type=int, nargs=1, default=[None],
                           help="Offset on j (do not read top j lines, useful for ORCA grids)")

    # ------------------------- Bathymetry and Sections -------------------------
    bathy_group = parser.add_argument_group("Bathy & Sections")
    bathy_group.add_argument("--bathyf", metavar='FILES', type=str, nargs='+', help="Bathy file(s)")
    bathy_group.add_argument("--bathyv", metavar='VARS', type=str, nargs=1, help="Bathy variable")
    bathy_group.add_argument("--bathylvl", metavar='LEVELS', type=float, nargs='+', help="Bathy contour levels")
    bathy_group.add_argument("--secf", metavar='FILES', type=str, nargs='+', help="Section line files to plot")

    # ------------------------- Misc -------------------------
    misc_group = parser.add_argument_group("Miscellaneous")
    misc_group.add_argument("--mask", metavar='FILES', type=str, nargs='+', help="Mask file(s)")
    misc_group.add_argument("--mesh", metavar='FILES', type=str, nargs='+', help="Mesh file(s)")
    misc_group.add_argument("--llonce", metavar='0/1', type=int, nargs=1, default=[0],
                            choices=[0,1], help="Read lat/lon for each plot (1=yes)")
    misc_group.add_argument("--debug", action="store_true", help="Enable debug mode with detailed output")

    return parser.parse_args()

def fix_list(arg_list, nfile, name):
    """
    Ensure argument list matches the required length.

    If one element is given, replicate it nfile times. If None or empty:
    - If required=True, raise an error.
    - Otherwise, fill with default (replicated nfile times).

    Parameters:
    - arg_list (list): The argument list to validate.
    - nfile (int): The required length of the list.
    - name (str): The name of the argument (for error messages).

    Returns:
    - list: The validated and adjusted argument list.

    Raises:
    - ValueError: If the argument list length is invalid.
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

def sanity_checks_args(args):
    """
    Perform sanity checks and fix inconsistencies in the parsed arguments.

    Parameters:
    - args (argparse.Namespace): Parsed arguments as a namespace.

    Returns:
    - argparse.Namespace: The updated arguments with fixed inconsistencies.
    """
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

    return args

def get_config(args):
    """
    Build the configuration dictionary from the parsed arguments.

    Parameters:
    - args (argparse.Namespace): Parsed arguments as a namespace.

    Returns:
    - dict: The configuration dictionary for the PyChart application.
    """
    # --- Build configuration dictionary ---
    pychart_config = {
        "figure": {
            "title": args.ft[0],
            "projection": args.p,
            "projbox": args.projbox,
            "sp": args.sp[0],
            "output": args.o[0],
            "crs": args.crs[0],
            "ploc": args.ploc,
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
            "spfid": args.spfid,
            "sprid": args.sprid,
            "offsety": args.joffset[0],
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
            "offsety": args.joffset[0],
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
    return pychart_config

def main():
    """
    The main entry point for the PyChart CLI.

    Parses arguments, performs sanity checks, builds the configuration, and generates the figure.
    """
    print("")

    # prepare config dictionary from input arguments
    args = parse_args()

    # set debug mode if needed
    set_debug(args.debug)

    # sanity checks on input arguments
    args = sanity_checks_args(args)

    # prepare config dictionary from input arguments
    pychart_config = get_config(args)

    # --- Build the figure using the new class ---
    figure_config = pychart_config["figure"]
    map_config = pychart_config["map"]
    cnt_config = pychart_config["cnt"]
    cb_config = pychart_config["cb"]

    fb = FigureBuilder(config=figure_config)
    print(fb)                         # print figure summary
    
    fb.build_layout()                 # builds figure + subplots + titles

    # --- Create colorbar instance ---
    map_cb = cb(cb_config)
    print(map_cb)
    
    # --- Loop through subplots ---
    for iax, ax in enumerate(fb.axes):
        info(f"Processing subplot {iax+1}/{len(fb.axes)}")
        print()
        # MAP
        if iax < len(map_config["files"]):
            debug(map_config)
            pcol = add_map_plot(map_config, map_cb, iax, ax, fb.proj[iax])

        # CONTOUR
        if (cnt_config["files"] and (iax < len(cnt_config["files"]))):
            debug(cnt_config["files"])
            add_cnt_plot(cnt_config, iax, ax)

    # --- Add colorbar using the class method ---
    fb.add_colorbar(pcol, map_cb)

    fb.save_figure()

    info(f"Display figure ...")
    plt.show()