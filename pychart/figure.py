import matplotlib.pyplot as plt
import numpy as np
import yaml
import cartopy.crs as ccrs
import cartopy.feature as cfeature

class FigureBuilder:
    def __init__(self, config, projections_yaml=None):
        self.config = config
        self.fig = None
        self.axes = []
        self.ploc = config.get("ploc", None)

        # Load projections if a YAML file is provided
        self.projections = {}
        if projections_yaml:
            with open(projections_yaml, "r") as f:
                self.projections = yaml.safe_load(f)["projections"]

        # Initialize projection and extent if the projection is available as key in config
        self.proj = []
        self.extent = []
        if "projection" in config and self.projections:
            self._init_projection(config["projection"])

    def _init_projection(self, proj_names):
        """Initialize projection and extent from the YAML configuration."""

        for iproj, proj_name in enumerate(proj_names):
            if proj_name not in self.projections:
                raise ValueError(f"Projection '{proj_name}' not found in the YAML file.")
            proj_data = self.projections[proj_name]
            proj_class = getattr(ccrs, proj_data["projection"])
            params = proj_data.get("params", {})
            extent = proj_data["extent"]
            if extent[-1] == "cproj":
                extent[-1] = self.proj
            self.extent.append(extent)
            self.proj.append(proj_class(**params))

    def build_layout(self):
        nisplt, njsplt = self._get_subplot(self.config["sp"])
        fig = plt.figure(figsize=np.array([297, 297*(njsplt+1)/nisplt]) / 25.4)
        gs = fig.add_gridspec(njsplt, nisplt)

        pltloc=[]
        if self.ploc:
            nplt = len(self.ploc)
            for iplt in range(nplt):
                # Convert the string to a list of integers
                pltpos = list(map(int, self.ploc[iplt].strip('[]').split(',')))
                pltloc.append(gs[pltpos[2]:pltpos[3]+1, pltpos[0]:pltpos[1]+1])
        else:
            pltloc.append(gs[iplt // nisplt, iplt % nisplt])
            nplt = len(pltloc)

        for iplt in range(nplt):
            print(f"Creating subplot {iplt+1}/{nplt} with projection {self.proj[iplt]}")
            print(pltloc[iplt])
            ax = fig.add_subplot(pltloc[iplt], projection=self.proj[iplt])
            extent= self.extent[iplt]
            if extent:
                if extent[0] == "global":
                    ax.set_global()
                else:
                    if extent[1] == "PlateCarree":
                        ax.set_extent(extent[0], crs=ccrs.PlateCarree())
                    else:
                        ax.set_extent(extent[0], extent[1])
            ax.gridlines()
            self.add_land_features(ax, ['isf', 'lakes', 'land'])
            self.axes.append(ax)

        self._add_title(self.config["title"])
            # remove extra white space
        hpx=0.06+0.035*njsplt
        fig.subplots_adjust(left=0.05,right=0.88, bottom=0.06, top=0.85, wspace=0.1, hspace=hpx)
        self.fig = fig
        
        return fig, self.axes

    def add_land_features(self, ax, features):
        feature_map = {
            'isf': cfeature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none'),
            'lakes': cfeature.NaturalEarthFeature('physical', 'lakes', '50m', facecolor='none'),
            'coast': cfeature.NaturalEarthFeature('physical', 'coastline', '50m', facecolor='none'),
            'land': cfeature.NaturalEarthFeature('physical', 'land', '50m', facecolor='none')
        }
        for f in features:
            ax.add_feature(feature_map[f], linewidth=0.5, edgecolor='k')

    def add_colorbar(self, pcol, map_cb, fontsize=16, cboffset=0.02, cbw=0.02):
        boxxy = self._get_figure_bounds()
        cax = plt.axes([boxxy[2]+cboffset, boxxy[1], cbw, boxxy[3]-boxxy[1]])
        cbar = plt.colorbar(pcol, cax=cax, format=map_cb.fmt, extend=map_cb.ext)
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.ax.set_title(map_cb.unit,fontsize=fontsize,y=1.0)
        return cbar

    # --- helper methods ---

    def _get_subplot(self, csubplt):
        ni, nj = map(int, csubplt.split('x'))
        return ni, nj

    def _get_figure_bounds(self):
        x0, x1, y0, y1 = 1.0, 0.0, 1.0, 0.0
        for ax in self.axes:
            ax.apply_aspect()
            box = ax.get_position()
            x0, x1 = min(x0, box.x0), max(x1, box.x1)
            y0, y1 = min(y0, box.y0), max(y1, box.y1)
        return [x0, y0, x1, y1]

    def _add_title(self, title, yoffset=0.05, height=0.1):
        boxxy = self._get_figure_bounds()
        cax = plt.axes([boxxy[0], boxxy[3] + yoffset, boxxy[2]-boxxy[0], height])
        cax.text(0.5, 0.5, title, ha='center', va='bottom', fontsize=20)
        cax.axis('off')
