import cartopy

class ax:
    def __init__(self,ctitle,cproj,lfeature,ploc,fig):
        self._feature=lfeature
        self._title=ctitle
        self._proj=None
        self._xylim=None
        self.box = [None, None, None, None]
        self.ax = fig.add_subplot(ploc, projection=self._proj)
        
    def _feature(self):
        dfeature={'isf':cartopy.feature.NaturalEarthFeature('physical', 'antarctic_ice_shelves_polys', '50m', facecolor='none'),
                  'lakes':cartopy.feature.NaturalEarthFeature('physical', 'lakes'                    , '50m', facecolor='none'),
                  'coast':cartopy.feature.NaturalEarthFeature('physical', 'coastline'                , '50m', facecolor='0.75'),
                  'land' :cartopy.feature.NaturalEarthFeature('physical', 'land'                     , '50m', facecolor='0.75'),
                  'bathy_z1000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_J_1000'  , '10m', facecolor='none'),
                  'bathy_z2000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_I_2000'  , '10m', facecolor='none'),
                  'bathy_z3000':cartopy.feature.NaturalEarthFeature('physical', 'bathymetry_H_3000'  , '10m', facecolor='none')
                 }

        for _,cfeat in enumerate(self._feature):
            self._ax.add_feature(dfeature[cfeat],linewidth=0.5,edgecolor='k')

    @property
    def ax(self):
        return self._ax

    @ax.setter
    def ax(self,pax):
        self._ax = pax

    def _extent(self):
        # put proj, extend, grid ...
        if self._xylim[0] == 'global':
            self._ax.set_global()
        else:
            self._ax.set_extent(self._xylim[0], self._xylim[1])

    def define_ax_properties(self)
        self._extent()
        self._feature()
        self._ax.gridlines(linewidth=1, color='k', linestyle='--')
        self._ax.set_title(self._title,fontsize=18)
        self._box=self._ax.get_position()

    @property
    def box(self):
        return self._box

    @box.setter
    def box(self,box):
        """ 
        find the coordinate of the axes

        Parameters
        ----------
        parameter 1: axe

        Returns
        -------
        output 1: list
            position box [x left, x right, y bottom, y top]

        Raises
        ------
        No raise
        """
        self._box=box
 

class fig:

    def __init__(ploc,ctitle,njsplt,nisplt,nplt) 
        # define figure dimension
        # define grid specification
        # define title
        self._fig = plt.figure(figsize=np.array([297,297*njsplt/nisplt]) / 25.4)
        self._title = ctitle
        self.pltloc = self.compute_pltloc(njsplt,nisplt,nplt)

    def compute_pltloc(self,njsplt,nisplt,nplt)

        pltloc=[None]*nplt

        gs = self._fig.add_gridspec(njsplt, nisplt)
        for iplt in range(0,nplt):
            if args.ploc :
                pltloc[iplt] = eval('gs['+args.ploc[ifile]+']')
            else:
                pltloc[iplt] = gs[iplt//nisplt,iplt%nisplt]

        return pltloc

# ============================ output argument list in txt file ================
def save_as(self,cfile):
    """
    save command line as string in a text file.
    save figure.

    Parameters
    ----------
    parameter 1: string
        file name used for the figure and command line record (without extension)

    Returns
    -------
    None

    Raises
    ------
    No raise
    """

    # save argument list
    fid = open(cfile+'.txt',"w")
    fid.write(' python '+subprocess.list2cmdline(sys.argv[0:])+'\n')
    fid.close()

    # save figure
    self._fig.savefig(cfile+'.png', format='png', dpi=150)

def add_title(self,boxxy):
    cax  = plt.axes([boxxy[0], boxxy[3], boxxy[2]-boxxy[0], 1-boxxy[3]])
    cax.text(0.5,0.5,ctitle,horizontalalignment='center',verticalalignment='bottom',fontsize=20)
    cax.axis('off')

def get_plt_bound(self,ax_lst):
    """
    find the box containing of a list of axes

    Parameters
    ----------
    parameter 1: list of axes

    Returns
    -------
    output 1: list
        position box [x left, x right, y bottom, y top]

    Raises
    ------
    No raise
    """
    x0=1.0; x1=0.0; y0=1.0; y1=0.0
    for iplt in range(0,len(ax_lst)):
        ax_lst[iplt].ax.apply_aspect()
        box = ax_lst[iplt].box
        x0=np.min([x0,box[0]])
        x1=np.max([x1,box[1]])
        y0=np.min([y0,box[2]])
        y1=np.max([y1,box[3]])
    return [x0, y0, x1, y1]


