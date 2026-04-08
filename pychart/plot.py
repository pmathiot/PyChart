# Refactored PyChart plotting module using DataConfig
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Optional, Sequence, Any
import numpy as np
import xarray as xr

from .log import info, debug
from .grids import StructuredGrid, TriGrid, IcoGrid, GridStrategy
from .fig_utils import get_lvls

# -----------------------------
# DataConfig: describes a single dataset (run or ref)
# -----------------------------

@dataclass
class DataConfig:
    def __str__(self):
        lines = []
        lines.append("DataConfig:")
        lines.append(f"    files    : {self.file}")
        lines.append(f"    vars     : {self.var}")
        lines.append(f"    jk, kt   : {self.jk}, {self.kt}")
        lines.append(f"    joffset  : {self.joffset}")
        lines.append(f"    scale    : {self.scale}")
        lines.append(f"    title    : {self.title}")
        return "\n".join(lines)
    
    file: str
    var: str
    jk: int = 1
    kt: int = 1
    scale: float = 1.0
    joffset: Optional[int] = None
    title: Optional[str] = None

    # internal cached data (not part of init / repr)
    _data: Optional[np.ndarray] = dataclasses.field(default=None, init=False, repr=False)

    def load(self, grid: GridStrategy) -> np.ndarray:
        """Load and cache the raw data array using provided grid strategy.
        This method will call `grid.load_data(self)`. It is idempotent.
        """
        if self._data is None:
            self._data = grid.load_data(self)
        else:
            debug(f"Using cached data for {self.file}:{self.var}")
        return self._data

    @property
    def data(self) -> np.ndarray:
        """Return cached data. Raise if not loaded yet (explicit is better).
        Use `load(grid)` to populate the cache before accessing.
        """
        if self._data is None:
            raise RuntimeError("Data not loaded: call load(grid) first")
        return self._data

# -----------------------------
# PlotConfig: a plot may involve a run + optional reference
# -----------------------------

@dataclass
class PlotConfig:
    def __str__(self):
        lines = []
        lines.append("PlotConfig:")

        # --- Simulation field
        run_str = str(self.run).replace("\n", "\n    ")
        lines.append("  Simulation field:")
        lines.append(f"    {run_str}")
        lines.append(f"")

        # --- Reference field
        if self.ref is not None:
            ref_str = str(self.ref).replace("\n", "\n    ")
            lines.append("  Reference field:")
            lines.append(f"    {ref_str}")
            lines.append(f"")
            lines.append(f"  Operation: {self.refop}")
            lines.append(f"")

        lines.append(f"  Line/Color Levels: {self.lvls}")
        return "\n".join(lines)
    
    run: DataConfig
    ref: Optional[DataConfig] = None
    refop: Optional[str] = None
    lvls: Optional[Sequence[float]] = None

# -----------------------------
# PlotData orchestration
# -----------------------------

class PlotData:
    """Orchestrates coordinates and data used for plotting.

    Attributes:
    cfg: PlotConfig
    grid: GridStrategy
    lon, lat: coordinates
    data_to_plot: final array to be used by plotting backends (after scaling/refop)
    """
    def __str__(self):
        lines = []
        lines.append("PlotData:")
        lines.append(f"    {self.cfg}")
        lines.append("")

        # --- Grid
        grid_name = self.grid.__class__.__name__
        lines.append(f"  Grid type   : {grid_name}")
        lines.append("")
        
        # --- Loaded status
        lines.append(f"  lon/lat     : {'loaded' if self.lon is not None else 'not loaded'}")
        lines.append(f"  data        : {'loaded' if self.cfg.run.data is not None else 'not loaded'}")
        if self.cfg.ref is not None:
            lines.append(f"  dataref     : {'loaded' if self.cfg.ref.data is not None else 'not loaded'}")
        else:
            lines.append(f"  dataref     : not loaded")

        return "\n".join(lines)
    
    def __init__(self, cfg: PlotConfig, grid: GridStrategy):
        self.cfg = cfg
        self.grid = grid

        self.lon: Optional[np.ndarray] = None
        self.lat: Optional[np.ndarray] = None

        # final data for plotting (run scaled / run - ref / run / ref)
        self.data_to_plot: Optional[np.ndarray] = None

    def load(self):
        """Load coordinates and underlying raw datasets (lazy DataConfig.load used).
        Assumes run (and ref if provided) share same grid.
        """
        info("Loading coordinates and data ...")

        # coords from run (ref supposed to share same grid)
        self.lon, self.lat = self.grid.load_coords(self.cfg.run)

        # Load run data
        self.cfg.run.load(self.grid)

        # Load reference data if present
        if self.cfg.ref is not None:
            # defensive check: ensure ref has values
            if not getattr(self.cfg.ref, 'file', None) or not getattr(self.cfg.ref, 'var', None):
                info("Reference config present but missing file/var — skipping ref load")
            else:
                self.cfg.ref.load(self.grid)

    def compute(self):
        """Compute `data_to_plot` from run and optional reference according to refop.
        Applies scales (DataConfig.scale) to each input.
        """
        print('')
        info("Computing data ...")

        if self.cfg.ref is not None and self.cfg.refop is not None:
            run_arr = self.cfg.run.data
            ref_arr = self.cfg.ref.data
            op = self.cfg.refop

            if op == "-":
                self.data_to_plot = (self.cfg.run.scale * run_arr) - (self.cfg.ref.scale * ref_arr)
            elif op == "/":
                # avoid division by zero
                denom = (self.cfg.ref.scale * ref_arr)
                with np.errstate(divide='ignore', invalid='ignore'):
                    self.data_to_plot = (self.cfg.run.scale * run_arr) / denom
            else:
                raise ValueError(f"Unsupported refop: {op}")
        else:
            # simple scaled run data
            self.data_to_plot = self.cfg.run.data * self.cfg.run.scale

    def plot_map(self, ax, map_cb, proj=None, **kwargs):
        """Delegate plotting to grid implementation using data_to_plot."""
        if self.data_to_plot is None:
            raise RuntimeError("data_to_plot not computed — call load() and compute() first")
        print('')
        info("Plotting map ...")
        return self.grid.plot_map(ax, self, map_cb, proj=proj, **kwargs)

    def plot_contour(self, ax, levels=10, **kwargs):
        if self.data_to_plot is None:
            raise RuntimeError("data_to_plot not computed — call load() and compute() first")
        return self.grid.plot_contour(ax, self, levels=levels, **kwargs)

    def add_title(self, ax):
        if self.cfg.ref is not None and self.cfg.refop is not None:
            left = self.cfg.run.title or self.cfg.run.var
            right = self.cfg.ref.title or self.cfg.ref.var
            ax.set_title(f"{left} {self.cfg.refop} {right}", fontsize=16)
        else:
            ax.set_title(self.cfg.run.title or self.cfg.run.var, fontsize=16)

# -----------------------------
# Grid detection
# -----------------------------

def detect_grid_strategy(ncfile: str) -> GridStrategy:
    ds = xr.open_dataset(ncfile)
    try:
        dimnames = [d.lower() for d in ds.sizes.keys()]
        vertex_dims = [d for d in dimnames if "vertex" in d]

        if vertex_dims:
            nvertex = ds.sizes[vertex_dims[0]]
            if nvertex == 3:
                return TriGrid()
            if nvertex == 6:
                return IcoGrid()
            if nvertex == 4:
                return StructuredGrid()
            raise ValueError("Unsupported unstructured grid")

        return StructuredGrid()
    finally:
        ds.close()

# -----------------------------
# Factory
# -----------------------------

def create_plotdata_from_config(cfg: PlotConfig) -> PlotData:
    grid = detect_grid_strategy(cfg.run.file)
    return PlotData(cfg, grid)

# -----------------------------
# Public API
# -----------------------------

def add_map_plot(map_config: dict, map_cb: Any, iax: int, ax, proj=None):
    # build DataConfig for run
    run_cfg = DataConfig(
        file=map_config["files"][iax],
        var=map_config["vars"][iax],
        jk=int(map_config["jk"][iax]) if map_config["jk"][iax] else 1,
        kt=int(map_config["jt"][iax]) if map_config["jt"][iax] else 1,
        scale=float(map_config.get("scale", [1.0])[iax]),
        joffset=map_config.get("offsety", None),
        title=map_config.get("spfid", [None])[iax],
    )

    # reference config
    ref_cfg = None
    if map_config.get("refs")[iax]:
        ref_cfg = DataConfig(
            file=map_config["refs"][iax],
            var=map_config["ref_vars"][iax],
            jk=int(map_config["jk"][iax]) if map_config["jk"][iax] else 1,
            kt=int(map_config["jt"][iax]) if map_config["jt"][iax] else 1,
            scale=float(map_config.get("ref_scale", [1.0])[iax]),
            joffset=map_config.get("offsety", None),
            title=map_config.get("sprid", [None])[iax],
        )

    cfg = PlotConfig(
        run=run_cfg,
        ref=ref_cfg,
        refop=map_config.get("op", [None])[iax],
        lvls=getattr(map_cb, "lvls", None)
    )
    print('--------- Add Colormap ---------')
    print(cfg)
    print('')

    pd = create_plotdata_from_config(cfg)

    if proj is None and isinstance(pd.grid, (TriGrid, IcoGrid)):
        raise ValueError("--noproj is only valid for structured grids")

    pd.load()
    pd.compute()

    # handle levels on the callback if needed
    lvls = getattr(map_cb, "lvls", None)
    if lvls is None or any(v is None for v in lvls):
        # compute from final data_to_plot
        data_min = float(np.nanmin(pd.data_to_plot))
        data_max = float(np.nanmax(pd.data_to_plot))
        map_cb.get_lvl([data_min, data_max])
        map_cb.compute_norm()
        print('')
        if map_cb.cnorm == "Normalize":
            info(f" Auto-computed min/max: {data_min} to {data_max}")            
        else:
            info(f" Auto-computed levels: {map_cb.lvls}")
        print('')
    
    artist = pd.plot_map(ax, map_cb, proj=proj)
    pd.add_title(ax)
    return artist


def add_cnt_plot(cnt_config: dict, iax: int, ax):
    run_cfg = DataConfig(
        file=cnt_config["files"][iax],
        var=cnt_config["vars"][iax],
        jk=int(cnt_config["jk"][iax]) if cnt_config["jk"][iax] else 1,
        kt=int(cnt_config["jt"][iax]) if cnt_config["jt"][iax] else 1,
        scale=float(cnt_config.get("scale", [1.0])[iax]),
        title=None,
    )

    ref_cfg = None
    if cnt_config.get("refs")[iax]:
        ref_cfg = DataConfig(
            file=cnt_config["refs"][iax],
            var=cnt_config["ref_vars"][iax],
            jk=int(cnt_config["jk"][iax]) if cnt_config["jk"][iax] else 1,
            kt=int(cnt_config["jt"][iax]) if cnt_config["jt"][iax] else 1,
            scale=float(cnt_config.get("ref_scale", [1.0])[iax]),
        )

    cfg = PlotConfig(
        run=run_cfg,
        ref=ref_cfg,
        refop=cnt_config.get("op", [None])[iax],
        lvls=cnt_config.get("levels")
    )

    if cfg.lvls:
        cfg.lvls = get_lvls(cfg.lvls)
    else:
        cfg.lvls = np.linspace(float(pd.data.min()), float(pd.data.max()), 10)
    
    print('')
    print('--------- Add Contours ---------')
    print(cfg)
    print('')

    pd = create_plotdata_from_config(cfg)

    pd.load()

    pd.compute()

    return pd.plot_contour(ax, levels=cfg.lvls, colors='k', linewidths=1)
