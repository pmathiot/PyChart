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

usage: plot_maps.py [-h] [--dir data dir] --mapf file_name [file_name ...]
                    --mapv var_name [var_name ...] [--mapreff file_ref]
                    [--maprefv reference var_name]
                    [--mapsf map data scale factor] [--mapjk vertical level]
                    [--mapz depth of the map] [--cbn color map name] --cblvl
                    color range [color range ...] [--cbu color map unit]
                    [--cbfmt color bar fmt] [--cbext color bar extend]
                    [--ft figure title] [--spfid runid [runid ...]]
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
