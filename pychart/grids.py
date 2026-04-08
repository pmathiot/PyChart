# grid.py
from __future__ import annotations

import abc
import numpy as np
import xarray as xr
import matplotlib.tri as mtri
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs
from cartopy.crs import Stereographic, NorthPolarStereo, SouthPolarStereo
from scipy.interpolate import griddata
from typing import Optional, Tuple, Sequence

from .io_utils import get_latlon, get_latlon_var, get_2d_data, get_name
from .log import info, debug

# -----------------------------
# Grid strategy base class
# -----------------------------

class GridStrategy(abc.ABC):
    """Abstract base class describing a grid strategy.

    Subclasses must implement methods for loading coords, loading data and plotting.
    """

    @abc.abstractmethod
    def load_coords(self, cfg: "DataConfig") -> np.ndarray:
        """Return lon, lat arrays appropriate for plotting (1D or 2D as needed)."""

    @abc.abstractmethod
    def load_data(self, cfg: "DataConfig") -> np.ndarray:
        """Return data arrays."""

    @abc.abstractmethod
    def plot_map(self, ax, pd: "PlotData", map_cb, proj=None, **kwargs):
        """Plot the data on axis and return the artist used for colorbar linking."""

    @abc.abstractmethod
    def plot_contour(self, ax, pd: "PlotData", levels=10, **kwargs):
        """Plot contours for the field and return contour set."""

# -----------------------------
# Structured grid
# -----------------------------

class StructuredGrid(GridStrategy):
    def load_coords(self, cfg: "DataConfig") -> Tuple[np.ndarray, np.ndarray]:
        lon2d, lat2d = get_latlon(cfg.file, offsety=cfg.joffset)
        return lon2d, lat2d

    def load_data(self, cfg: "DataConfig") -> np.ndarray:
        data = get_2d_data(cfg.file, cfg.var, klvl=cfg.jk, ktime=cfg.kt, offsety=cfg.joffset)
        return data

    def plot_map(self, ax, pd: "PlotData", map_cb, proj=None, **kwargs):
        data = pd.data_to_plot
        if data is None or pd.lon is None or pd.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")

        info("Plotting structured pcolor map ...")
        debug(f"Plotting pcolor map with shape {data.shape}, lon shape {pd.lon.shape}, lat shape {pd.lat.shape}")

        if proj is not None:
            pcm = ax.pcolormesh(pd.lon, pd.lat, data,
                                cmap=map_cb.cmap, norm=map_cb.norm, transform=ccrs.PlateCarree(),
                                **kwargs)
        else:
            pcm = ax.pcolormesh(data, cmap=map_cb.cmap, norm=map_cb.norm, **kwargs)
        return pcm

    def plot_contour(self, ax, pd: "PlotData", levels=10, **kwargs):
        data = pd.data_to_plot
        if data is None or pd.lon is None or pd.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")
        cs = ax.contour(pd.lon, pd.lat, data, levels=levels, transform=ccrs.PlateCarree(), **kwargs)
        return cs

# -----------------------------
# Triangular unstructured grid
# -----------------------------

class TriGrid(GridStrategy):
    """Triangular unstructured grid strategy."""

    def load_coords(self, cfg: "DataConfig") -> Tuple[np.ndarray, np.ndarray]:
        ds = xr.open_dataset(cfg.file)
        try:
            clat, clon = get_latlon_var(cfg.file)
            lon = ds[clon].values
            lat = ds[clat].values
        finally:
            ds.close()
        return lon, lat

    def load_data(self, cfg: "DataConfig") -> np.ndarray:
        ds = xr.open_dataset(cfg.file)
        try:
            data = ds[cfg.var].isel(time=cfg.kt - 1).values
        finally:
            ds.close()
        return data

    def _load_extra(self, cfg: "PlotConfig") -> np.ndarray:
        ds = xr.open_dataset(cfg.file)
        try:
            ctriangles = get_name("(.*_face_nodes)", ds.variables.keys())
            triangles = ds[ctriangles].values
        finally:
            ds.close()
        return triangles

    def plot_map(self, ax, pd: "PlotData", map_cb, proj=None, **kwargs):
        if pd.data_to_plot is None or pd.lon is None or pd.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")

        if getattr(pd, 'tri', None) is None:
            pd.tri = self._load_extra(pd.cfg.run)

        xy = ax.projection.transform_points(ccrs.Geodetic(), pd.lon, pd.lat)
        x, y = xy[:, 0], xy[:, 1]

        is_global = not isinstance(ax.projection, (Stereographic, NorthPolarStereo, SouthPolarStereo))
        if is_global:
            central_lon = ax.projection.proj4_params.get('lon_0', 0.0)
            lon_shifted = ((pd.lon - central_lon + 180) % 360) - 180
            seam_mask = np.ptp(lon_shifted[pd.tri], axis=1) > 180
        else:
            seam_mask = np.zeros(len(pd.tri), dtype=bool)

        triang = mtri.Triangulation(x, y, pd.tri)
        triang.set_mask(seam_mask)

        plot_data = pd.data_to_plot

        pcm = ax.tripcolor(x, y, pd.tri[seam_mask==False], plot_data[seam_mask==False], cmap=map_cb.cmap, norm=map_cb.norm, rasterized=True, **kwargs)
        return pcm

    def plot_contour(self, ax, pd: "PlotData", levels=10, **kwargs):
        if pd.data_to_plot is None or pd.lon is None or pd.lat is None:
            raise ValueError("Data and coordinates must be loaded before plotting.")

        if getattr(pd, 'tri', None) is None:
            pd.tri = self._load_extra(pd.cfg.run)

        xy = ax.projection.transform_points(ccrs.Geodetic(), pd.lon, pd.lat)
        x, y = xy[:, 0], xy[:, 1]

        is_global = not isinstance(ax.projection, (Stereographic, NorthPolarStereo, SouthPolarStereo))
        if is_global:
            central_lon = ax.projection.proj4_params.get('lon_0', 0.0)
            lon_shifted = ((pd.lon - central_lon + 180) % 360) - 180
            seam_mask = np.ptp(lon_shifted[pd.tri], axis=1) > 180
        else:
            seam_mask = np.zeros(len(pd.tri), dtype=bool)

        triang = mtri.Triangulation(x, y, pd.tri)
        triang.set_mask(seam_mask)

        data_nodal = self.elemental_to_nodal(pd) if len(pd.data_to_plot) == len(pd.tri) else pd.data_to_plot

        cs = ax.tricontour(x, y, pd.tri[seam_mask==False], data_nodal, levels=levels, **kwargs)
        return cs

    def elemental_to_nodal(self,pd):
        """ 
        Convert elemental (triangle-centered) data to nodal (vertex) data.

        Returns:
        - numpy.ndarray: Data values at each node (averaged from connected triangles).
        """ 
        n_nodes = pd.tri.max() + 1  # because node is 0 based
        nodal_sum = np.zeros(n_nodes)
        nodal_count = np.zeros(n_nodes)

        for i, tri in enumerate(pd.tri):
            for node in tri:
                nodal_sum[node] += pd.data_to_plot[i]
                nodal_count[node] += 1

        nodal_data = nodal_sum / nodal_count
        return nodal_data


# -----------------------------
# Icosaheral grid
# -----------------------------

class IcoGrid(GridStrategy):
    """Icosaheral grid strategy using per-cell polygon bounds."""

    def load_coords(self, cfg: "PlotConfig") -> Tuple[np.ndarray, np.ndarray]:
        ds = xr.open_dataset(cfg.file)
        try:
            lon = ds['lon'].values
            lat = ds['lat'].values
            bnds_lon = ds['bounds_lon'].values
            bnds_lat = ds['bounds_lat'].values
        finally:
            ds.close()
        return lon, lat

    def load_data(self, cfg: "PlotConfig") -> np.ndarray:
        ds = xr.open_dataset(cfg.file)
        try:
            data = ds[cfg.var].isel(time=cfg.kt - 1).values
        finally:
            ds.close()
        return data

    def _load_bounds(self, cfg: "PlotConfig") -> Tuple[np.ndarray, np.ndarray]:
        ds = xr.open_dataset(cfg.file)
        try:
            bnds_lon = ds['bounds_lon'].values
            bnds_lat = ds['bounds_lat'].values
        finally:
            ds.close()
        return bnds_lon, bnds_lat

    def plot_map(self, ax, pd: "PlotData", map_cb, proj=None, **kwargs):
        if pd.data is None:
            raise ValueError("Data must be loaded before plotting.")

        if getattr(pd, 'bnds_lon', None) is None or getattr(pd, 'bnds_lat', None) is None:
            pd.bnds_lon, pd.bnds_lat = self._load_bounds(pd.cfg)

        xy = ax.projection.transform_points(ccrs.Geodetic(), pd.bnds_lon, pd.bnds_lat)
        x = xy[:, :, 0]
        y = xy[:, :, 1]

        lon_min, lon_max, lat_min, lat_max = ax.get_extent(crs=ccrs.PlateCarree())
        poly_lon_min = pd.bnds_lon.min(axis=1)
        poly_lon_max = pd.bnds_lon.max(axis=1)
        poly_lat_min = pd.bnds_lat.min(axis=1)
        poly_lat_max = pd.bnds_lat.max(axis=1)
        visible_mask = (
            (poly_lon_max >= lon_min) & (poly_lon_min <= lon_max) &
            (poly_lat_max >= lat_min) & (poly_lat_min <= lat_max)
        )

        is_global = not isinstance(ax.projection, (Stereographic, NorthPolarStereo, SouthPolarStereo))
        if is_global:
            central_lon = ax.projection.proj4_params.get('lon_0', 0.0)
            bnds_lon_shifted = ((pd.bnds_lon - central_lon + 180) % 360) - 180
            seam_mask = np.ptp(bnds_lon_shifted, axis=1) > 180
        else:
            seam_mask = np.zeros_like(visible_mask, dtype=bool)

        final_mask = visible_mask & (~seam_mask)
        debug(f"→ Final polygons kept: {np.sum(final_mask)} of {len(final_mask)}")

        indices = np.where(final_mask)[0]
        polys_visible = [list(zip(x[i, :], y[i, :])) for i in indices]
        data_visible = pd.data[indices]

        collection = PolyCollection(polys_visible, array=data_visible, cmap=map_cb.cmap, norm=map_cb.norm,
                                    edgecolor='k', linewidth=0.1, rasterized=True)
        ax.add_collection(collection)
        return collection

    def plot_contour(self, ax, pd: "PlotData", levels=10, **kwargs):
        if pd.lon is None or pd.lat is None:
            raise ValueError("Coordinates must be loaded before contour plotting.")

        lon = pd.lon
        lat = pd.lat
        data = pd.data

        lon_min, lon_max, lat_min, lat_max = ax.get_extent(crs=ccrs.PlateCarree())
        visible_mask = (
            (lon >= lon_min) & (lon <= lon_max) &
            (lat >= lat_min) & (lat <= lat_max)
        )

        lon_vis = lon[visible_mask]
        lat_vis = lat[visible_mask]
        data_vis = data[visible_mask]

        if lon_vis.size < 3:
            raise ValueError("Not enough points to contour in the visible area")

        xy = ax.projection.transform_points(ccrs.PlateCarree(), lon_vis, lat_vis)
        x, y = xy[:, 0], xy[:, 1]

        nx, ny = 300, 300
        x_grid = np.linspace(x.min(), x.max(), nx)
        y_grid = np.linspace(y.min(), y.max(), ny)
        x2d, y2d = np.meshgrid(x_grid, y_grid)

        grid_data = griddata((x, y), data_vis, (x2d, y2d), method='cubic')
        cs = ax.contour(x2d, y2d, grid_data, levels=levels, **kwargs)
        return cs
