import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

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

    def __init__(self,cmap,cnorm,cunit,cfmt,cext,lvls):
        self._unit=cunit
        self._fmt=cfmt
        self._ext=cext

        if cnorm == 'TwoSlopeNorm':
            self._lvls=[lvls[0],lvls[2],lvls[1]]
        else:
            self._lvls=get_lvl(lvls)

        self.cmap = plt.get_cmap(cmap)
        self.norm = self._compute_norm(cnorm)

    @property
    def norm(self):
        return self._norm

    @norm.setter
    def norm(self,norm):
        self._norm = norm

    @property
    def cmap(self):
        return self._cmap

    @cmap.setter
    def cmap(self,cmap):
        self._cmap = cmap

    def _compute_norm(self,cnorm):
        dnorm={
               'BoundaryNorm':colors.BoundaryNorm(self._lvls, self.cmap.N, extend=self._ext),
               'LogNorm':colors.LogNorm(vmin=self._lvls[0],vmax=self._lvls[-1]),
               'Normalize':colors.Normalize(vmin=self._lvls[0],vmax=self._lvls[-1]),
               'TwoSlopeNorm':colors.TwoSlopeNorm(vmin=self._lvls[0],vmax=self._lvls[2],vcenter=self._lvls[1])
              }
        return dnorm[cnorm]

    def add_colorbar(self,pcol,boxxy,fontsize=16,cboffset=0.02,cbw=0.02):
        cax  = plt.axes([boxxy[2]+cboffset, boxxy[1], cbw, boxxy[3]-boxxy[1]])
        cbar = plt.colorbar(pcol, cax=cax, format=self._fmt,extend=self._ext)
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.ax.set_title(self._unit,fontsize=fontsize,y=1.0)
        return cax, cbar

    def __str__(self):
        return 'cb : {}, {}, {}'.format(self.unit,self.fmt,self.ext)

def get_lvl(bnds):
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

    return lvl


