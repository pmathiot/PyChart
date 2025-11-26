import numpy as np
import yaml
import cartopy.crs as ccrs
import xarray as xr
from matplotlib.collections import PolyCollection
from scipy.interpolate import griddata
from .io_utils import get_2d_data, get_latlon_var, get_2d_data
from .log import info, debug
import matplotlib.tri as tri
from cartopy.crs import Stereographic, NorthPolarStereo, SouthPolarStereo

def plot_cartesian(ax, data, title="Plot"):
    """
    Plot data on a Cartesian grid.

    Parameters:
    - ax (matplotlib.axes.Axes): The axis to plot on.
    - data (numpy.ndarray): The data to plot.
    - title (str, optional): Title of the plot. Default is "Plot".

    Returns:
    - matplotlib.image.AxesImage: The image object created by imshow.
    """
    im = ax.imshow(data, cmap="viridis", origin="lower")
    ax.set_title(title)
    return im

def add_map_plot(map_config, map_cb, iax, ax, proj=None):
    """
    Add a map plot to the given axis.

    Parameters:
    - map_config (dict): Configuration for the map plot.
    - map_cb (object): Colorbar configuration.
    - iax (int): Index of the axis.
    - ax (matplotlib.axes.Axes): The axis to plot on.
    - proj (cartopy.crs.Projection, optional): Projection to use. Default is None.

    Returns:
    - matplotlib.collections.QuadMesh: The plotted map object.
    """
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
        trun=map_config["spfid"][iax],
        tref=map_config["sprid"][iax],
        refop=map_config["op"][iax]
    )
    print("File to map:")
    print(map_data)
    print()

    if (not proj) and (map_data.type != 'structured'):
        raise ValueError("--noproj option is only valid for structured grids.")

    if map_data.type == 'tri_unstructured':
        map_data = map_data.to_triunstructured()
    elif map_data.type == 'ico_unstructured':
        map_data = map_data.to_icounstructured()
    else:
        map_data = map_data.to_structured()

    map_data.get_coords()
    map_data.get_data()
    map_data.compute_data()

    # update lvls for colorbar if not provided
    if map_cb.lvls[0] is [None]:
        info(" Setting colorbar levels from data min/max ...")
        map_cb.get_lvl([map_data.data.min(), map_data.data.max()])
        map_cb.compute_norm()
        info(f" Colorbar min/max set to {map_cb.lvls[0]:g} / {map_cb.lvls[-1]:g}")

    pcol = map_data.plot_map(ax, map_cb, proj)  # returns QuadMesh or similar for colorbar
    map_data.add_title(ax)

    return  pcol

def add_cnt_plot(cnt_config, iax, ax):
    """
    Add a contour plot to the given axis.

    Parameters:
    - cnt_config (dict): Configuration for the contour plot.
    - iax (int): Index of the axis.
    - ax (matplotlib.axes.Axes): The axis to plot on.
    """
    cnt_data = PlotData(
        file=cnt_config["files"][iax],
        var=cnt_config["vars"][iax],
        lvls=cnt_config["levels"],
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
    cntlvl=cb.get_lvl(cnt_data.lvls)
    cnt_data.plot_cnt(ax, levels=cntlvl, colors='k', linewidths=1)

class PlotData:
    """
    Class to handle plotting data and metadata.

    Attributes:
    - file (str): Main data file.
    - var (str): Variable name in the file.
    - lvls (list, optional): Contour levels.
    - jk (int): Depth level.
    - kt (int): Time frame.
    - fileref (str, optional): Reference file.
    - varref (str, optional): Reference variable.
    - jkref (int, optional): Reference depth level.
    - ktref (int, optional): Reference time frame.
    - refop (str, optional): Reference operation.
    - data (numpy.ndarray, optional): Loaded data.
    - dataref (numpy.ndarray, optional): Loaded reference data.
    - lon (numpy.ndarray, optional): Longitude coordinates.
    - lat (numpy.ndarray, optional): Latitude coordinates.
    - sf (float): Scale factor for the main data.
    - sfref (float): Scale factor for the reference data.
    - trun (str, optional): Title for the run ID.
    - tref (str, optional): Title for the variable.
    - type (str): Data type ('structured' or 'unstructured').
    """
    def __str__(self):
        lines = []
        lines.append("  Simulation field:")
        lines.append(f"    files    : {self.file}")
        lines.append(f"    vars     : {self.var}")
        if self.lvls is not None:
            lines.append(f"    levels   : {self.lvls}")
        lines.append(f"    jk       : {self.jk}")
        lines.append(f"    kt       : {self.kt}")
        lines.append(f"    scale    : {self.sf}")
        lines.append(f"    type     : {self.type}")
        lines.append(f"    name     : {self.trun}")

        # Reference field section
        if self.fileref is not None:
            lines.append(" ")
            lines.append("  Reference field:")
            lines.append(f"    fileref  : {self.fileref}")
            lines.append(f"    varref   : {self.varref}")
            lines.append(f"    jkref    : {self.jkref}")
            lines.append(f"    ktref    : {self.ktref}")
            lines.append(f"    refop    : {self.refop}")
            lines.append(f"    scale    : {self.sfref}")
            lines.append(f"    name     : {self.tref}")

        return "\n".join(lines)
    
    def __init__(self, file, var, lvls=None, jk=1, kt=1, fileref=None, varref=None, jkref=None, ktref=None, sf=1.0, sfref=1.0, trun=None, tref=None, refop=None):
        """
        Initialize the PlotData object.

        Parameters:
        - file (str): Main data file.
        - var (str): Variable name in the file.
        - lvls (list, optional): Contour levels. Default is None.
        - jk (int, optional): Depth level. Default is 1.
        - kt (int, optional): Time frame. Default is 1.
        - fileref (str, optional): Reference file. Default is None.
        - varref (str, optional): Reference variable. Default is None.
        - jkref (int, optional): Reference depth level. Default is None.
        - ktref (int, optional): Reference time frame. Default is None.
        - sf (float, optional): Scale factor for the main data. Default is 1.0.
        - sfref (float, optional): Scale factor for the reference data. Default is 1.0.
        - trun (str, optional): Title for the run ID. Default is None.
        - tref (str, optional): Title for the variable. Default is None.
        - refop (str, optional): Reference operation. Default is None.
        """
        self.file = file          # Main data file
        self.var = var            # Variable name in the file
        self.lvls = lvls
        self.jk = jk              # Depth level (if applicable)
        self.kt = kt              # Time frame (if applicable)
        self.fileref = fileref    # Reference file (optional)
        self.varref = varref      # Reference variable (optional)
        self.jkref = jkref        # Reference depth level (if applicable)
        self.ktref = ktref        # Reference time frame (if applicable)
        self.refop = refop        # Reference operation (if applicable)
        self.data = None          # Loaded data
        self.dataref = None       # Loaded reference data
        self.lon = None           # Longitude coordinates
        self.lat = None           # Latitude coordinates
        self.sf = sf            # Scale factor for the main data
        self.sfref = sfref         # Scale factor for the reference data
        self.trun = trun    # Title for the run ID
        self.tref = tref      # Title for the variable

        self.type = self._init_type()  # Data type: 'structured' or 'unstructured'

    def _init_type(self):
        """
        Detect the grid type based on dimension names and sizes.

        Returns:
        - str: Grid type ('structured', 'tri_unstructured', or 'ico_unstructured').
        """
        ds = xr.open_dataset(self.file)
        dimnames = [d.lower() for d in ds.sizes.keys()]
        grid_type = "structured"

        # Look for any dimension that includes "vertex"
        vertex_dims = [d for d in dimnames if "vertex" in d]
        debug(f"Detected vertex dimensions: {vertex_dims}, {ds.sizes[vertex_dims[0]]}")
        if vertex_dims:
            # Take the first one (usually only one)
            vertex_dim = vertex_dims[0]
            nvertex = ds.sizes[vertex_dim]

            if nvertex == 3:
                grid_type = "tri_unstructured"
            elif nvertex == 6:
                grid_type = "ico_unstructured"
            elif nvertex == 4:
                grid_type = "structured"
            else:
                raise ValueError("This unstructured grid type is unsupported (not 3 or 6 vertex)")

        return grid_type

    def add_title(self, ax):
        """
        Add a title to the axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to add the title to.
        """
        if self.tref and self.refop:
            ctitle=f"{self.trun} {self.refop} {self.tref}"
        elif self.trun:
            ctitle=f"{self.trun}"
        else:
            ctitle=f"{self.var}"
        ax.set_title(ctitle,fontsize=18)

    def compute_data(self):
        """
        Compute the data to be plotted, applying the operation with reference data if providinfo(f" Computing data ...")
ed.
        info(f" Computing data ...")
"""
        info(f" Computing data ...")
        if self.dataref is not None:
            if self.refop == "-":
                self.data = (self.data * self.sf) - (self.dataref * self.sfref)
            elif self.refop == "/":
                self.data = (self.data * self.sf) / (self.dataref * self.sfref)
        else:
            self.data = self.data * self.sf

    def to_structured(self):
        """
        Convert the data to a structured grid format.

        Returns:
        - StructuredPlotData: The structured plot data object.
        """
        return StructuredPlotData(
            file=self.file,
            var=self.var,
            lvls = self.lvls,
            jk=self.jk,
            kt=self.kt,
            fileref=self.fileref,
            varref=self.varref,
            jkref=self.jkref,
            ktref=self.ktref,
            sf=self.sf,
            sfref=self.sfref,
            trun=self.trun,
            tref=self.tref,
            refop=self.refop
        )
    
    def to_triunstructured(self):
        """
        Convert the data to a triangular unstructured grid format.

        Returns:
        - TriUnStructuredPlotData: The triangular unstructured plot data object.
        """
        return TriUnStructuredPlotData(
            file=self.file,
            var=self.var,
            lvls=self.lvls,
            jk=self.jk,
            kt=self.kt,
            fileref=self.fileref,
            varref=self.varref,
            jkref=self.jkref,
            ktref=self.ktref,
            sf=self.sf,
            sfref=self.sfref,
            trun=self.trun,
            tref=self.tref,
            refop=self.refop
        )
    
    def to_icounstructured(self):
        """
        Convert the data to an icosahedral unstructured grid format.

        Returns:
        - IcoUnStructuredPlotData: The icosahedral unstructured plot data object.
        """
        return IcoUnStructuredPlotData(
            file=self.file,
            var=self.var,
            lvls=self.lvls,
            jk=self.jk,
            kt=self.kt,
            fileref=self.fileref,
            varref=self.varref,
            jkref=self.jkref,
            ktref=self.ktref,
            sf=self.sf,
            sfref=self.sfref,
            trun=self.trun,
            tref=self.tref,
            refop=self.refop
        )

class StructuredPlotData(PlotData):
    """
    Class for handling structured grid plot data.
    """
    def __init__(self, file, var, lvls=None, jk=None, kt=None, fileref=None, varref=None, jkref=None, ktref=None, sf=1.0, sfref=1.0, trun=None, tref=None, refop=None):
        """
        Initialize the StructuredPlotData object.

        Parameters:
        - file (str): Main data file.
        - var (str): Variable name in the file.
        - lvls (list, optional): Contour levels. Default is None.
        - jk (int, optional): Depth level. Default is None.
        - kt (int, optional): Time frame. Default is None.
        - fileref (str, optional): Reference file. Default is None.
        - varref (str, optional): Reference variable. Default is None.
        - jkref (int, optional): Reference depth level. Default is None.
        - ktref (int, optional): Reference time frame. Default is None.
        - sf (float, optional): Scale factor for the main data. Default is 1.0.
        - sfref (float, optional): Scale factor for the reference data. Default is 1.0.
        - trun (str, optional): Title for the run ID. Default is None.
        - tref (str, optional): Title for the variable. Default is None.
        - refop (str, optional): Reference operation. Default is None.
        """
        super().__init__(file, var, lvls, jk, kt, fileref, varref, jkref, ktref, sf, sfref, trun, tref, refop)

    def get_data(self, joffset=-2):
        """
        Load the main data and 
        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
reference data (
        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
if applicable).

        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
        """
        self.data = get_2d_data(self.file, self.var, klvl=self.jk, ktime=self.kt, offsety=joffset)
        if self.fileref and self.varref:
            self.dataref = get_2d_data(self.fileref, self.varref, klvl=self.jkref, ktime=self.ktref, offsety=joffset)

    def get_coords(self, mesh_file=None, joffset=0):
        """
        Load latitude and longitude coordinates.

        Parameters:
        - mesh_file (str, optional): Path to the mesh file. Default is None.
        - joffset (int, optional): Offset for the y-dimension. Default is 0.
        """
        if mesh_file:
            self.lon, self.lat = self.get_latlon()
        else:
            self.lon, self.lat = self.get_latlon()

    def get_latlon(self):
        """
        Retrieve latitude and longitude data from the file.

        Returns:
        - tuple: Longitude and latitude arrays.
        """
        clat,clon=get_latlon_var(self.file)
        lat2d=get_2d_data(self.file,clat,offsety=-2)
        lon2d=get_2d_data(self.file,clon,offsety=-2)
        lon2d[lon2d>=180] = lon2d[lon2d>=180.] - 360.
        delta_lon=np.abs(np.diff(lon2d))
        for i, start in enumerate(np.argmax(delta_lon > 180, axis=1)):
            lon2d[i, start+1:] += 360
        return lon2d, lat2d

    def plot_map(self, ax, map_cb, proj=None, **kwargs):
        """
        Plot data on a map using the given axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to plot on.
        - map_cb (object): Colorbar configuration.
        - proj (cartopy.crs.Projection, optional): Projection to use. Default is None.

        Returns:
        - matplotlib.collections.QuadMesh: The plotted map object.
        """
        if self.data is None or self.lon is None or self.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        
        info(" Plotting pcolor map ...")
        debug(f" Plotting pcolor map with shape {self.data.shape}, lon shape {self.lon.shape}, lat shape {self.lat.shape}")
        print("")
        if proj:
            # Use transform only if a projection is defined
            pcm = ax.pcolormesh(self.lon, self.lat, self.data, cmap=map_cb.cmap, norm=map_cb.norm, transform=ccrs.PlateCarree(), **kwargs)
        else:
            # Plot without transform for Cartesian coordinates
            pcm = ax.pcolormesh(self.data, cmap=map_cb.cmap, norm=map_cb.norm, **kwargs)
        return pcm

    def plot_cnt(self, ax, levels=10, **kwargs):
        """
        Plot contour lines on the given axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to plot on.
        - levels (int or list, optional): Contour levels. Default is 10.
        """
        if self.data is None or self.lon is None or self.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        cs = ax.contour(self.lon, self.lat, self.data, levels=levels, transform=ccrs.PlateCarree(), **kwargs)
        return cs
    
class TriUnStructuredPlotData(PlotData):
    """
    Class for handling triangular unstructured grid plot data.
    """
    def __init__(self, file, var, lvls=None, jk=None, kt=None, fileref=None, varref=None, jkref=None, ktref=None, sf=1.0, sfref=1.0, trun=None, tref=None, refop=None):
        """
        Initialize the TriUnStructuredPlotData object.

        Parameters:
        - file (str): Main data file.
        - var (str): Variable name in the file.
        - lvls (list, optional): Contour levels. Default is None.
        - jk (int, optional): Depth level. Default is None.
        - kt (int, optional): Time frame. Default is None.
        - fileref (str, optional): Reference file. Default is None.
        - varref (str, optional): Reference variable. Default is None.
        - jkref (int, optional): Reference depth level. Default is None.
        - ktref (int, optional): Reference time frame. Default is None.
        - sf (float, optional): Scale factor for the main data. Default is 1.0.
        - sfref (float, optional): Scale factor for the reference data. Default is 1.0.
        - trun (str, optional): Title for the run ID. Default is None.
        - tref (str, optional): Title for the variable. Default is None.
        - refop (str, optional): Reference operation. Default is None.
        """
        super().__init__(file, var, lvls, jk, kt, fileref, varref, jkref, ktref, sf, sfref, trun, tref, refop)

    def get_data(self, joffset=-2):
        """
        Load the main data and reference data (if applicable).
        """
        ds = xr.open_dataset(self.file)
        self.data = ds[self.var].isel(time=self.kt - 1).values
        if self.fileref and self.varref:
            ref_ds = xr.open_dataset(self.fileref)
            self.dataref = ref_ds[self.varref].isel(time=self.ktref - 1).values

    def get_coords(self, mesh_file=None, joffset=0):
        """
        Load latitude and longitude coordinates.
        """
        ds = xr.open_dataset(self.file)
        self.lon = ds['lon'].values
        self.lat = ds['lat'].values
        self.tri = ds['triangles'].values
        ds.close()

    def plot_map(self, ax, map_cb, proj=None, **kwargs):
        """
        Plot data on a map using the given axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to plot on.
        - map_cb (object): Colorbar configuration.
        - proj (cartopy.crs.Projection, optional): Projection to use. Default is None.

        Returns:
        - matplotlib.collections.QuadMesh: The plotted map object.
        """
        if self.data is None or self.lon is None or self.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        debug(f" Plotting data on map with shape {self.data.shape}, lon shape {self.lon.shape}, lat shape {self.lat.shape}")
        # Precompute projection of nodes
        xy = ax.projection.transform_points(ccrs.Geodetic(), self.lon, self.lat)
        x, y = xy[:, 0], xy[:, 1]

        # Check if projection is a global longitude-based projection
        is_global = not isinstance(ax.projection, (Stereographic, NorthPolarStereo, SouthPolarStereo))

        if is_global:
            # Only for global projections, compute seam mask
            central_lon = ax.projection.proj4_params.get('lon_0', 0.0)
            lon_shifted = ((self.lon - central_lon + 180) % 360) - 180
            lontri = lon_shifted[self.tri]
            mask = np.ptp(lontri, axis=1) > 180
        else:
            # For stereo projections, no seam masking needed
            mask = np.zeros(len(self.tri), dtype=bool)

        triang = tri.Triangulation(x, y, self.tri)
        triang.set_mask(mask)

        #pcm = ax.tripcolor(triang, self.data, cmap=map_cb.cmap, norm=map_cb.norm, rasterized=True,  **kwargs)
        pcm = ax.tripcolor(x, y, self.tri[mask==False], self.data[mask==False], cmap=map_cb.cmap, norm=map_cb.norm, rasterized=True,  **kwargs)
        #pcm = ax.tripcolor(self.lon, self.lat, self.tri, self.data, cmap=map_cb.cmap, norm=map_cb.norm, transform=ccrs.PlateCarree(), rasterized=True,  **kwargs)
        return pcm

    def plot_cnt(self, ax, levels=10, **kwargs):
        """
        Plot contour lines on the given axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to plot on.
        - levels (int or list, optional): Contour levels. Default is 10.
        """
        if self.data is None or self.lon is None or self.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        debug(f" Plotting data on cnt with shape {self.data.shape}, lon shape {self.lon.shape}, lat shape {self.lat.shape}")

        if len(self.data) != len(self.lon) :
            data_nodal = self.elemental_to_nodal()
        else :
            data_nodal = self.data

        # Transform lon/lat to the axis projection coordinates
        xy = ax.projection.transform_points(ccrs.Geodetic(), self.lon, self.lat)
        x, y = xy[:, 0], xy[:, 1]

        # Check if projection is a global longitude-based projection
        is_global = not isinstance(ax.projection, (Stereographic, NorthPolarStereo, SouthPolarStereo))

        if is_global:
            # Only for global projections, compute seam mask
            central_lon = ax.projection.proj4_params.get('lon_0', 0.0)
            lon_shifted = ((self.lon - central_lon + 180) % 360) - 180
            lontri = lon_shifted[self.tri]
            mask = np.ptp(lontri, axis=1) > 180
        else:
            # For stereo projections, no seam masking needed
            mask = np.zeros(len(self.tri), dtype=bool)

        triang = tri.Triangulation(x, y, self.tri)
        triang.set_mask(mask)

        # Plot tricontour
        #cs = ax.tricontour(triang, data_nodal, levels=levels, **kwargs)
        cs = ax.tricontour(x, y, self.tri[mask==False], data_nodal, levels=levels, **kwargs)
        #cs = ax.tricontour(self.lon, self.lat, self.tri, data_nodal, transform=ccrs.PlateCarree(), levels=levels, **kwargs)
        return cs
    
    def elemental_to_nodal(self):
        """
        Convert elemental (triangle-centered) data to nodal (vertex) data.

        Returns:
        - numpy.ndarray: Data values at each node (averaged from connected triangles).
        """
        print("Converting elemental data to nodal data...")
        n_nodes = self.tri.max() + 1
        nodal_sum = np.zeros(n_nodes)
        nodal_count = np.zeros(n_nodes)

        for i, tri in enumerate(self.tri):
            for node in tri:
                nodal_sum[node] += self.data[i]
                nodal_count[node] += 1

        nodal_data = nodal_sum / nodal_count
        return nodal_data

class IcoUnStructuredPlotData(PlotData):
    """
    Class for handling icosahedral unstructured grid plot data.
    """
    def __init__(self, file, var, lvls, jk=None, kt=None, fileref=None, varref=None, jkref=None, ktref=None, sf=1.0, sfref=1.0, trun=None, tref=None, refop=None):
        """
        Initialize the IcoUnStructuredPlotData object.

        Parameters:
        - file (str): Main data file.
        - var (str): Variable name in the file.
        - lvls (list, optional): Contour levels. Default is None.
        - jk (int, optional): Depth level. Default is None.
        - kt (int, optional): Time frame. Default is None.
        - fileref (str, optional): Reference file. Default is None.
        - varref (str, optional): Reference variable. Default is None.
        - jkref (int, optional): Reference depth level. Default is None.
        - ktref (int, optional): Reference time frame. Default is None.
        - sf (float, optional): Scale factor for the main data. Default is 1.0.
        - sfref (float, optional): Scale factor for the reference data. Default is 1.0.
        - trun (str, optional): Title for the run ID. Default is None.
        - tref (str, optional): Title for the variable. Default is None.
        - refop (str, optional): Reference operation. Default is None.
        """
        super().__init__(file, var, lvls, jk, kt, fileref, varref, jkref, ktref, sf, sfref, trun, tref, refop)

    def get_data(self, joffset=-2):
        """
        Load the mrefer 
        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
nceand refer 
        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
nce data (if 
        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
applicable).

        Parameters:
        - joffset (int, optional): Offset for the y-dimension. Default is -2.
        """
        ds = xr.open_dataset(self.file)
        self.data = ds[self.var].isel(time=self.kt - 1).values
        if self.fileref and self.varref:
            ref_ds = xr.open_dataset(self.fileref)
            self.dataref = ref_ds[self.varref].isel(time=self.ktref - 1).values

    def get_coords(self):
        """
        Load latitude and longitude coordinates.
        """
        ds = xr.open_dataset(self.file)
        self.lon = ds['lon'].values
        self.lat = ds['lat'].values
        self.bnds_lon = ds['bounds_lon'].values
        self.bnds_lat = ds['bounds_lat'].values

    def plot_map(self, ax, map_cb, proj=None, **kwargs):
        """
        Plot data on a map using the given axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to plot on.
        - map_cb (object): Colorbar configuration.
        - proj (cartopy.crs.Projection, optional): Projection to use. Default is None.

        Returns:
        - matplotlib.collections.PolyCollection: The plotted map object.
        """
        if self.data is None or self.lon is None or self.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        debug(f" Plotting data on map with shape {self.data.shape}, lon shape {self.lon.shape}, lat shape {self.lat.shape}")
 
        # --- 1️⃣ Project polygon vertices once ---
        xy = ax.projection.transform_points(ccrs.Geodetic(), self.bnds_lon, self.bnds_lat)
        x = xy[:, :, 0]
        y = xy[:, :, 1]

        # --- 2️⃣ Get visible geographic extent ---
        lon_min, lon_max, lat_min, lat_max = ax.get_extent(crs=ccrs.PlateCarree())
        info(f"Visible extent: lon=[{lon_min:.1f}, {lon_max:.1f}], lat=[{lat_min:.1f}, {lat_max:.1f}]")

        # --- 3️⃣ Fast visibility pre-filter ---
        poly_lon_min = self.bnds_lon.min(axis=1)
        poly_lon_max = self.bnds_lon.max(axis=1)
        poly_lat_min = self.bnds_lat.min(axis=1)
        poly_lat_max = self.bnds_lat.max(axis=1)

        visible_mask = (
            (poly_lon_max >= lon_min) & (poly_lon_min <= lon_max) &
            (poly_lat_max >= lat_min) & (poly_lat_min <= lat_max)
        )
        debug(f"→ Keeping {np.sum(visible_mask)} of {len(visible_mask)} polygons visible")

        # --- 4️⃣ Compute seam mask only for global projections ---
        is_global = not isinstance(ax.projection, (Stereographic, NorthPolarStereo, SouthPolarStereo))
        if is_global:
            central_lon = ax.projection.proj4_params.get('lon_0', 0.0)
            bnds_lon_shifted = ((self.bnds_lon - central_lon + 180) % 360) - 180
            seam_mask = np.ptp(bnds_lon_shifted, axis=1) > 180
        else:
            seam_mask = np.zeros_like(visible_mask, dtype=bool)

        # --- 5️⃣ Combine both masks ---
        final_mask = visible_mask & (~seam_mask)
        debug(f"→ Final polygons kept: {np.sum(final_mask)}")

        # --- 6️⃣ Build visible polygons and data ---
        visible_indices = np.where(final_mask)[0]
        polys_visible = [list(zip(x[i, :], y[i, :])) for i in visible_indices]
        data_visible = self.data[visible_indices]

        # --- 7️⃣ Create PolyCollection ---
        collection = PolyCollection(
            polys_visible,
            array=data_visible,
            cmap=map_cb.cmap,
            norm=map_cb.norm,
            edgecolor='k',
            linewidth=0.1,
            rasterized=True
        )
        ax.add_collection(collection)
        return collection
    
    def plot_cnt(self, ax, levels=10, **kwargs):
        """
        Plot contour lines on the given axis.

        Parameters:
        - ax (matplotlib.axes.Axes): The axis to plot on.
        - levels (int or list, optional): Contour levels. Default is 10.
        """
        if self.data is None or self.lon is None or self.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        debug(f" Plotting data on cnt with shape {self.data.shape}, lon shape {self.lon.shape}, lat shape {self.lat.shape}")

        # 1️⃣ Get geographic extent in lat/lon
        lon_min, lon_max, lat_min, lat_max = ax.get_extent(crs=ccrs.PlateCarree())

        # 2️⃣ Filter points inside visible geographic box
        visible_mask = (
            (self.lon >= lon_min) & (self.lon <= lon_max) &
            (self.lat >= lat_min) & (self.lat <= lat_max)
        )
        lon_vis = self.lon[visible_mask]
        lat_vis = self.lat[visible_mask]
        data_vis = self.data[visible_mask]

        # 3️⃣ Project filtered points
        xy = ax.projection.transform_points(ccrs.PlateCarree(), lon_vis, lat_vis)
        x, y = xy[:,0], xy[:,1]

        # 4️⃣ Build a regular grid in projected coordinates
        x_grid = np.linspace(x.min(), x.max(), 500)   # smaller number may be enough
        y_grid = np.linspace(y.min(), y.max(), 500)
        x2d, y2d = np.meshgrid(x_grid, y_grid)

        # 5️⃣ Interpolate and contour
        grid_data = griddata((x, y), data_vis, (x2d, y2d), method='cubic')
        cs = ax.contour(x2d, y2d, grid_data, levels=levels, **kwargs)

        return cs
