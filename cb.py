import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# ============================ CMAP ====================================

class cb:

    def __init__(self,cmap,cnorm,cunit,cfmt,cext,lvls):
        self.unit=cunit
        self.fmt=cfmt
        self.ext=cext
        self.cmap=plt.get_cmap(cmap)

        self.__set_norm(lvls,cnorm)

    def __get_norm(self):
        return self.__norm

    def __set_norm(self,bnds,cnorm):

        if bnds:
            lvl=get_lvl(bnds)
        else:
            print(' Need definition of levels (min,max) at least.')
            sys.exit(42)

        if cnorm == 'BoundaryNorm':
            clrsnorm = colors.BoundaryNorm(lvl, self.cmap.N, extend=self.ext)
        elif cnorm == 'LogNorm':
            clrsnorm = colors.LogNorm(vmin=lvl[0],vmax=lvl[-1])
        elif cnorm == 'Normalize':
            clrsnorm = colors.Normalize(vmin=lvl[0],vmax=lvl[-1])
        self.__norm = clrsnorm
    
    norm = property(__get_norm, __set_norm)

    def add_colorbar(self,pcol,boxxy,cunit='',fontsize=16,cboffset=0.02,cbw=0.02):
        cax  = plt.axes([boxxy[2]+cboffset, boxxy[1], cbw, boxxy[3]-boxxy[1]])
        cbar = plt.colorbar(pcol, cax=cax, format=self.fmt,extend=self.ext)
        cbar.ax.tick_params(labelsize=fontsize)
        cbar.ax.set_title(self.unit,fontsize=fontsize,y=1.0)
        return cax, cbar

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
    if len(bnds)==3 :
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvlint = bnds[2]
        lvlmax = lvlmin+round((lvlmax-lvlmin)/lvlint)*lvlint
        lvl= np.arange(lvlmin,lvlmax+0.000001,lvlint)
    elif len(bnds)==2:
        lvlmin = bnds[0]
        lvlmax = bnds[1]
        lvl=np.linspace(lvlmin, lvlmax, num=10)
    else:
        lvl=np.array(bnds[:])
    return lvl


