import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
from .log import info

import cmocean

# ============================ CMAP ====================================

class cb:
    """
    A class to represent a colormap.

    Attributes:
    - cmap (cmap object): Colormap used in pcolormesh.
    - norm (colors.Normalize object): Normalization used in pcolormesh.
    - lvls (list): Discrete levels defining the colormap.
    - unit (str): Units for the colorbar.
    - fmt (str): Format for the colorbar labels.
    - ext (str): Extend parameter for the colorbar.
    - cnorm (str): Normalization type ('BoundaryNorm', 'LogNorm', etc.).

    Methods:
    - __str__: Returns a string representation of the colormap.
    - __init__: Initializes the colormap with the given configuration.
    - compute_norm: Computes the normalization based on the levels and colormap.
    - get_lvl: Computes an array of discrete levels for the colormap.
    """

    def __str__(self):
        # Determine the colormap name (matplotlib or cmocean)
        cmap_name = getattr(self.cmap, 'name', str(self.cmap))

        # Summarize level range
        if self.lvls[0] != None :
            lvl_min = self.lvls[0]
            lvl_max = self.lvls[-1]
            lvl_count = len(self.lvls)
            clvl_info = f" from {lvl_min:g} to {lvl_max:g} ({lvl_count} levels)"
        else:
            lvl_min = lvl_max = lvl_count = "N/A"
            clvl_info = "to be determined from data"
        
        return (
            f"Colorbar:\n"
            f"  Colormap     : {cmap_name}\n"
            f"  Norm         : {self.cnorm}\n"
            f"  Levels       : {clvl_info}\n"
            f"  Extend       : {self.ext}\n"
            f"  Units        : {self.unit}\n"
            f"  Format       : {self.fmt}\n"
        )

    def __init__(self, config):
        """
        Initialize the colormap with the given configuration.

        Parameters:
        - config (dict): Configuration dictionary containing:
            - units (str): Units for the colorbar.
            - fmt (str): Format for the colorbar labels.
            - extend (str): Extend parameter for the colorbar.
            - norm (str): Normalization type ('BoundaryNorm', 'LogNorm', etc.).
            - colormap (str): Name of the colormap.
            - cmocean (bool): Whether to use a cmocean colormap.
            - levels (list): Discrete levels for the colormap.
        """
        # get map colorbar
        self.unit=config["units"]
        self.fmt=config["fmt"]
        self.ext=config["extend"]
        self.cnorm=config["norm"]
        self.norm=None
        self.lvls=[None]

        if config["cmocean"]:
            self.cmap=eval(config["colormap"])
        else:
            self.cmap = plt.get_cmap(config["colormap"])

        if config["levels"] != [None] :
            self.get_lvl(config["levels"])
            self.compute_norm()
        else:
            info(" Colorbar levels will be determined from data min/max.")
            info("   and norm set to 'Normalize'")
            print('')
            self.cnorm='Normalize'


    def compute_norm(self):
        """
        Compute the normalization based on the levels and colormap.

        Sets the `norm` attribute to the appropriate normalization object.
        """
        dnorm={
               'BoundaryNorm':colors.BoundaryNorm(self.lvls, self.cmap.N, extend=self.ext),
               'LogNorm':colors.LogNorm(vmin=self.lvls[0],vmax=self.lvls[-1]),
               'Normalize':colors.Normalize(vmin=self.lvls[0],vmax=self.lvls[-1]),
               'TwoSlopeNorm':colors.TwoSlopeNorm(vmin=self.lvls[0],vmax=self.lvls[2],vcenter=self.lvls[1])
              }
        self.norm = dnorm[self.cnorm]

    def get_lvl(self, bnds):
        """
        Compute an array of discrete levels used to define the colormap.

        Parameters:
        - bnds (list): List of levels with length 2, 3, or more.

        Returns:
        - np.array: Array of levels.

        Raises:
        - SystemExit: If the levels are not properly defined.
        """
        if self.cnorm == 'TwoSlopeNorm':
                lvl=[bnds[0],bnds[2],bnds[1]]
        else:
            if len(bnds)==2:
                lvlmin = bnds[0]
                lvlmax = bnds[1]
                lvl=np.linspace(lvlmin, lvlmax, num=10)
            elif len(bnds)==3 :
                lvlmin = bnds[0]
                lvlmax = bnds[1]
                lvlint = bnds[2]
                lvlmax = lvlmin+round((lvlmax-lvlmin)/lvlint)*lvlint
                lvl= np.arange(lvlmin,lvlmax+0.000001,lvlint)
            elif len(bnds) > 3:
                lvl=np.array(bnds[:])
            else:
                print(' Need definition of levels (min,max) at least.')
                sys.exit(42)
        self.lvls=lvl