# PyChart
PyChart is a tool to map from netcdf output

## Example output
Below you can find an exemple of the kind of output from PyChart:

<img src="./FIGURES/bsf_OPM006.png">

The command line to produce this plot is:

```
python ~/GIT/PyChart/plot_maps.py \
   --ft 'BSF (OPM006)' \
   --spfid '1979-1988'                            '1989-1998'                           '1999-2008'                            '2009-2018'                         \
   --mapf eORCA025.L121-OPM006_10y_y1979_psi.nc eORCA025.L121-OPM006_10y_y1989_psi.nc eORCA025.L121-OPM006_10y_y1999_psi.nc eORCA025.L121-OPM006_10y_y2009_psi.nc  \
   --cntf eORCA025.L121-OPM006_10y_y1979_psi.nc eORCA025.L121-OPM006_10y_y1989_psi.nc eORCA025.L121-OPM006_10y_y1999_psi.nc eORCA025.L121-OPM006_10y_y2009_psi.nc  \
   --mapv sobarstf \
   --cntv sobarstf \
   --cblvl  0 100 10 \
   --cntlvl 0 100 10 \
   --cbu Sv \
   --cbn magma_r \
   --mapsf 0.000001 \
   --cntsf 0.000001 \
   --cbext both \
   --mesh   mesh.nc     \
   --bathyf mesh.nc     \
   --bathyv bathy_metry \
   --bathylvl 1000 2000 3000 4000 \
   -p south_ocean \
   --sp 2x2 \
   -o bsf_OPM006
```

## Requirement

PyChart is tested with this conda environement (environment.yml):
```
name: PyChart
channels:
  - defaults
dependencies:
  - python=3.7.4
  - cartopy
  - gsw
  - scipy
  - netcdf4
  - dask
  - xarray
  - jupyter
  - seawater
```

This environement can be created via this command:
```
conda env create -f environment.yml
```
with environment.yml being the environement file described above.

## Usage

```
usage: plot_maps.py [-h] [--dir data_dir] --mapf pcolor_file_names
                    [pcolor_file_names ...] --mapv pcolor_var_names
                    [pcolor_var_names ...] [--mapreff pcolor_ref_file_name]
                    [--maprefv pcolor_ref_var_name]
                    [--mapsf pcolor_scale_factor] [--mapjk pcolor_jk_depth]
                    [--mapz pcolor_z_depth] [--cbn colormap_name] --cblvl
                    colorbar_range [colorbar_range ...] [--cbu colorbar_unit]
                    [--cbfmt colorbar_fmt] [--cbext colorbar_extend]
                    [--ft figure_title] [--spfid runid [runid ...]]
                    [--sprid refid]
                    [--mask mask file name [mask file name ...]]
                    [--mesh mesh file name [mesh file name ...]]
                    [--sp subplot disposition] [-o output name]
                    [-p projection] [--cntf contour file [contour file ...]]
                    [--cntv contour var] [--cntreff contour ref file]
                    [--cntrefv contour ref var]
                    [--cntsf contour data scale factor]
                    [--cntjk contour jk level] [--cntz contour depth in m]
                    [--cntlvl contour line level [contour line level ...]]
                    [--bathyf bathy file [bathy file ...]]
                    [--bathyv bathy var]
                    [--bathylvl contour line level [contour line level ...]]
                    [--secf section line file list  [section line file list  ...]]

optional arguments:
  -h, --help            show this help message and exit
  --dir data_dir        data dir
  --mapf pcolor_file_names [pcolor_file_names ...]
                        names of input files
  --mapv pcolor_var_names [pcolor_var_names ...]
                        variable list
  --mapreff pcolor_ref_file_name
                        names of ref files
  --maprefv pcolor_ref_var_name
                        reference variable name
  --mapsf pcolor_scale_factor
                        map data scale factor
  --mapjk pcolor_jk_depth
                        level in fortran convention
  --mapz pcolor_z_depth
                        depth of the map
  --cbn colormap_name   color map name
  --cblvl colorbar_range [colorbar_range ...]
                        color range
  --cbu colorbar_unit   colorbar unit
  --cbfmt colorbar_fmt  colorbar format
  --cbext colorbar_extend
                        colorbar extend
  --ft figure_title     title of the whole figure
  --spfid runid [runid ...]
                        runids (title + mesh name)
  --sprid refid         refids (title + mesh name)
  --mask mask file name [mask file name ...]
                        mask file name
  --mesh mesh file name [mesh file name ...]
                        mesh file name
  --sp subplot disposition
                        subplot disposition (ixj)
  -o output name        output name
  -p projection         projection
  --cntf contour file [contour file ...]
                        contour file list
  --cntv contour var    contour variable
  --cntreff contour ref file
                        contour reference file
  --cntrefv contour ref var 
                        contour reference variable
  --cntsf contour data scale factor
                        contour data scale factor
  --cntjk contour jk level
                        contour jk level
  --cntz contour depth in m
                        contour depth in m
  --cntlvl contour line level [contour line level ...]
                        contour line level
  --bathyf bathy file [bathy file ...]
                        bathy file
  --bathyv bathy var    contour variable
  --bathylvl contour line level [contour line level ...]
                        contour line level
  --secf section line file list  [section line file list  ...]
                        section file list describing section to plot
```
