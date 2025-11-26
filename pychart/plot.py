# Refactored PyChart plotting module using DataConfig
from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import Optional, Sequence, Any
import numpy as np
import xarray as xr

from .log import info, debug
from .grids import StructuredGrid, TriGrid, IcoGrid, GridStrategy

# -----------------------------
# DataConfig: describes a single dataset (run or ref)
# -----------------------------

@dataclass
class DataConfig:
    file: str
    var: str
    jk: int = 1
    kt: int = 1
    scale: float = 1.0
    title: Optional[str] = None

# -----------------------------
# PlotConfig: a plot may involve a run + optional reference
# -----------------------------

@dataclass
class PlotConfig:
    run: DataConfig
    ref: Optional[DataConfig] = None
    refop: Optional[str] = None
    lvls: Optional[Sequence[float]] = None

# -----------------------------
# PlotData orchestration
# -----------------------------

class PlotData:
    def __init__(self, cfg: PlotConfig, grid: GridStrategy):
        self.cfg = cfg
        self.grid = grid

        self.lon: Optional[np.ndarray] = None
        self.lat: Optional[np.ndarray] = None
        self.data: Optional[np.ndarray] = None
        self.dataref: Optional[np.ndarray] = None
        self.tri: Optional[np.ndarray] = None

    def load(self):
        info(" Loading coordinates and data ...")

        # load coords from run (ref supposed to share same grid)
        print(self.cfg.run)
        self.lon, self.lat = self.grid.load_coords(self.cfg.run)

        # load main field
        print('TOTO : ',self.cfg.run)
        self.data = self.grid.load_data(self.cfg.run)

        # load reference
        if self.cfg.ref.file and self.cfg.ref.var:
            self.dataref = self.grid.load_data(self.cfg.ref)

    def compute(self):
        info(" Computing data ...")

        if self.cfg.ref and self.cfg.refop:
            op = self.cfg.refop
            if op == "-":
                self.data = (
                    self.cfg.run.scale * self.data
                    - self.cfg.ref.scale * self.dataref
                )
            elif op == "/":
                self.data = (
                    self.cfg.run.scale * self.data
                    / (self.cfg.ref.scale * self.dataref)
                )
            else:
                raise ValueError(f"Unsupported refop: {op}")
        else:
            self.data = self.data * self.cfg.run.scale

    def plot_map(self, ax, map_cb, proj=None, **kwargs):
        return self.grid.plot_map(ax, self, map_cb, proj=proj, **kwargs)

    def plot_contour(self, ax, levels=10, **kwargs):
        return self.grid.plot_contour(ax, self, levels=levels, **kwargs)

    def add_title(self, ax):
        if self.cfg.ref and self.cfg.refop:
            ax.set_title(f"{self.cfg.run.title} {self.cfg.refop} {self.cfg.ref.title}")
        else:
            ax.set_title(self.cfg.run.title or self.cfg.run.var)

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
        title=map_config.get("spfid", [None])[iax],
    )

    # reference config
    ref_cfg = None
    if map_config.get("refs"):
        ref_cfg = DataConfig(
            file=map_config["refs"][iax],
            var=map_config["ref_vars"][iax],
            jk=int(map_config["jk"][iax]) if map_config["jk"][iax] else 1,
            kt=int(map_config["jt"][iax]) if map_config["jt"][iax] else 1,
            scale=float(map_config.get("ref_scale", [1.0])[iax]),
            title=map_config.get("sprid", [None])[iax],
        )

    cfg = PlotConfig(
        run=run_cfg,
        ref=ref_cfg,
        refop=map_config.get("op", [None])[iax],
        lvls=getattr(map_cb, "lvls", None)
    )

    info("File to map:")
    info(str(cfg))

    pd = create_plotdata_from_config(cfg)

    if proj is None and isinstance(pd.grid, (TriGrid, IcoGrid)):
        raise ValueError("--noproj is only valid for structured grids")

    pd.load()
    pd.compute()

    lvls = getattr(map_cb, "lvls", None)
    if lvls is None or any(v is None for v in lvls):
        map_cb.get_lvl([float(np.nanmin(pd.data)), float(np.nanmax(pd.data))])
        map_cb.compute_norm()

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
    if cnt_config.get("refs"):
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

    pd = create_plotdata_from_config(cfg)
    pd.load()
    pd.compute()

    if cfg.lvls:
        levels = cfg.lvls
    else:
        levels = np.linspace(float(pd.data.min()), float(pd.data.max()), 10)

    return pd.plot_contour(ax, levels=levels, colors='k', linewidths=1)
