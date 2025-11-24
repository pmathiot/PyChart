import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import sys
from pychart.log import info

import cmocean

# ============================ CMAP ====================================

class cb:
    """
    A class to represent a colormap

    Attributes
    ----------
    cmap : cmap object
        colormap used in pcolormesh
    norm : colors.norm object
        normalisation used in pcolormesh

    Methods
    -------
    add_colorbar:
        add a colorbar to the plot
    """

    def __str__(self):
        # Determine the colormap name (matplotlib or cmocean)
        cmap_name = getattr(self.cmap, 'name', str(self.cmap))

        # Summarize level range
        if self.lvls:
            lvl_min = self.lvls[0]
            lvl_max = self.lvls[-1]
            lvl_count = len(self.lvls)
            clvl_info = f" from {lvl_min:g} to {lvl_max:g} ({lvl_count} levels)"
        else:
            lvl_min = lvl_max = lvl_count = "N/A"
            clvl_info = "to be determined from data"
            norm_type = "None"
        
        return (
            f"Colorbar:\n"
            f"  Colormap     : {cmap_name}\n"
            f"  Norm         : {self.cnorm}\n"
            f"  Levels       : {clvl_info}\n"
            f"  Extend       : {self.ext}\n"
            f"  Units        : {self.unit}\n"
            f"  Format       : {self.fmt}\n"
        )

    def __init__(self,config):
        # get map colorbar
        self.unit=config["units"]
        self.fmt=config["fmt"]
        self.ext=config["extend"]
        self.cnorm=config["norm"]
        self.norm=None
        self.lvls=None

        if config["levels"]:
            self.get_lvl(config["levels"])
            self.compute_norm()
        else:
            info(" Colorbar levels will be determined from data min/max.")
            info("   and norm set to 'Normalize'")
            print('')
            self.cnorm='Normalize'

        if config["cmocean"]:
            self.cmap=eval(config["colormap"])
        else:
            self.cmap = plt.get_cmap(config["colormap"])

    @property
    def cmap(self):
        return self._cmap

    @cmap.setter
    def cmap(self,cmap):
        self._cmap = cmap

    def compute_norm(self):
        dnorm={
               'BoundaryNorm':colors.BoundaryNorm(self.lvls, self.cmap.N, extend=self.ext),
               'LogNorm':colors.LogNorm(vmin=self.lvls[0],vmax=self.lvls[-1]),
               'Normalize':colors.Normalize(vmin=self.lvls[0],vmax=self.lvls[-1]),
               'TwoSlopeNorm':colors.TwoSlopeNorm(vmin=self.lvls[0],vmax=self.lvls[2],vcenter=self.lvls[1])
              }
        self.norm = dnorm[self.cnorm]

    def get_lvl(self,bnds):
        """
        compute an array of discrete levels used to define the colormap based on a list of levels of length 2, 3 or more.
        If the list length is 2, 10 equidistance levels are computed from list[0] to list[1]
                            3, X levels are computed from list[0] to list[1] by a step of list[2]
                            >3, the levels are the one specified in the input parameters

        Parameters
        ----------
        parameter 1: list
            list of levels (length 2, 3 or more)

        Returns
        -------
        output 1: np.array
            array of levels

        Raises
        ------
        No raise
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


