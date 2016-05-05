from spm_data_layer import *
from spm_data_load import load_data as spm_load_data
from spm_data_inspect import inspect_channels as spm_inspect_channels
from spm_data_align import \
    align_correlate as spm_align_correlate, \
    align_forback as spm_align_forback, \
    align_offsets as spm_align_offsets

import tkinter, os, math, pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable    
#from scipy.ndimage import map_coordinates
from scipy import signal

"""module containing SPMdata class with auxiliary routines
"""

# TODO: do 2d rezu vzit v uvahu i rozliseni v z-smeru

def stack_spm(spm1, spm2, finer=True):
    """stack two SPMdata objects into a single SPMdata object
    'spmnew' and return it
    
    spm1, spm2 - input SPMdata structures
    finer - if data resolution for spm1 and spm2 differs,
        three situations may happen:
        if finer=None, then no resampling is performed;
        if finer=True, then data with lower resolution are
        upsampled so that they can be added to data with higher
        resolution;
        if finer=False, then data with higher resolution are
        downsampled so that they can be added to data with lower
        resolution;
        each layer contains attributes horstep and verstep, which
        represent data resolution in x- and y-direction; these
        attributes should not be modified, since they express
        resolution of the original data 
    """
     
    # channels to be in spmnew and number of layers in spmnew
    chanlist = list(set(
        list(spm1.channels.keys()) + list(spm2.channels.keys())))
    numlay = spm1.numlay + spm2.numlay
    
    # create spmnew
    spmnew = SPMdata()
    spmnew.layers = [None]*numlay
    
    # add all channels from chanlist into spmnew:
    for chan in chanlist:      
        # if chan in both...
        if chan in spm1.channels and chan in spm2.channels:        
            # old dimensions of arrays for chan
            n1xo, n1yo = spm1.arrshape[chan]
            n2xo, n2yo = spm2.arrshape[chan]

            # horizontal and vertical resolution
            x1step = spm1.xsteps[chan]
            x2step = spm2.xsteps[chan]
            
            y1step = spm1.ysteps[chan]
            y2step = spm2.ysteps[chan]
        
            # if no resampling...
            if finer is None:
                xstep = min(x1step, x2step)
                ystep = min(y1step, y2step)

                # new array dimensions
                n1x, n1y = n1xo, n1yo
                n2x, n2y = n2xo, n2yo
                                
                # arrays to be stacked
                arr1 = spm1.channels[chan]            
                arr2 = spm2.channels[chan]                
            else:
                # if upsampling...
                if finer:
                    xstep = min(x1step, x2step)
                    ystep = min(y1step, y2step)
                # if downsampling...
                else:
                    xstep = max(x1step, x2step)
                    ystep = max(y1step, y2step)
                
                # new array dimensions
                n1x = int(n1xo * x1step / xstep)
                n1y = int(n1yo * y1step / ystep)

                n2x = int(n2xo * x2step / xstep)
                n2y = int(n2yo * y2step / ystep)

                # resampled arrays to be stacked
                arr1 = gen_resample_data_3d(spm1.channels[chan],
                    n1xo, n1yo, spm1.numlay,
                    n1x,  n1y,  spm1.numlay)
                arr2 = gen_resample_data_3d(spm2.channels[chan],
                    n2xo, n2yo, spm2.numlay,
                    n2x,  n2y,  spm2.numlay)
                
            # we assume that spm1 and spm2 have the same units
            # for all channels they share
            spmnew.units[chan] = spm2.units[chan]
            spmnew.xsteps[chan] = xstep
            spmnew.ysteps[chan] = ystep   
                
            # allocate new array with appropriate dimensions
            xlim, ylim = max(n1x, n2x), max(n1y, n2y)
            arr = np.empty((numlay, xlim, ylim))
            arr.fill(np.nan)
                        
            # fill new array with valid data            
            arr[:len(spm1), :n1x, :n1y] = arr1.copy()        
            arr[len(spm1):, :n2x, :n2y] = arr2.copy()
            
        else:
            # if chan in spm1 only...
            if chan in spm1.channels:
                print("stack_spm: The second SPMdata structure does "
                    "not contain channel '{}'.".format(chan))    
                spmaux = spm1
                ranaux = slice(0, len(spm1))
            else:
                print("stack_spm: The first SPMdata structure does "
                    "not contain channel '{}'.".format(chan))
                spmaux = spm2
                ranaux = slice(len(spm1), numlay)
                
            # new dimensions are equal to old dimensions
            nx, ny = spmaux.arrshape[chan]
            
            # allocate new array with appropriate dimensions
            arr = np.empty((numlay, nx, ny))
            arr.fill(np.nan)
            
            # fill new array with valid data
            arr[ranaux, :nx, :ny] = spmaux.channels[chan]

            spmnew.units[chan] = spmaux.units[chan]
            spmnew.xsteps[chan] = spmaux.xsteps[chan]
            spmnew.ysteps[chan] = spmaux.ysteps[chan]
            
        # save new array into spmnew
        spmnew.setchan(chan, arr)
        
    # update data for layers from spm1
    for i, lay in enumerate(spm1.layers):
        spmnew.layers[i] = lay.copy(i, spmnew.channels)

    # update data for layers from spm2
    for i, lay in enumerate(spm2.layers):
        spmnew.layers[i + len(spm1)] = lay.copy(
            i + len(spm1), spmnew.channels)
    
    print(("stack_spm: New SPMdata object has {} layers and "
        "channels: {}").format(spmnew.numlay,
        list(spmnew.channels.keys())))
    
    return spmnew

def create_spm(paths=None, numlay=0):
    """create SPMdata structure with all necessary attributes
    
    paths - if None, then new SPMdata structure is generated with
        no data; otherwise 'paths' is passed to 'load_data' routine
    numlay - number of layers in the new structure; this parameter
        applies only when 'paths' is None
    """
    
    spmdata = SPMdata(paths)
    if paths is None:
        if numlay < 0:
            print(("create_spm: Invalid number"
                " of layers: {}.").format(numlay))
            del spmdata
            return
        spmdata.layers = [None]*numlay
        for i in range(numlay):
            spmdata.layers[i] = SPMdataLayer()
        # BEWARE: spmdata.layers = [SPMdataLayer()]*numlay cannot
        # be used since there is only ONE SPMdataLayer() allocated
        # and spmdata.layers then is a list of numlay  pointer to
        # the same SPMdataLayer() object
    
    return spmdata

def pickle_data(spmdata, filename):
    """export data to external file by mechanism provided by
    pickle
    
    spmdata - SPMdata structure to export
    filename - where to export the data
    """
    
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)


def unpickle_data(spmdatalist, filename):
    """import data from external file by mechanism provided by
    pickle
    
    spmdatalist - list whose first element will be assigned data
        from file
    filename - where to import data from
    """

    with open(filename, 'rb') as f:
        spmdatalist[0] = pickle.load(f)

class SPMdata:
    """class embracing all the relevant data for SPM images
    """
    
    def __init__(self, paths=None):
    
        # dictionary of channels
        # self.channels = {ELEV:self.elev (any scan), ...}
        self.channels = {}
        
        # list of layers
        # self.layers = [<SPMdataLayer object 1>, ...]
        self.layers = []
        
        # resolution of data in x- and y-direction for all channels
        self.ysteps, self.xsteps = {}, {}
        
        # dictionary of units for each channel
        # it is possible that each SPMdata structure has its
        # own ad hoc channels with their own names and units
        self.units = UNITS
        
        # BEWARE: self.layers[i].channels[chan] is a view of
        # self.channels[chan][i, 
        #                     xoffind:xoffind + xnum,
        #                     yoffind:yoffind + ynum
        #                     ];
        # in other words, self.layers[i].channels[chan] points
        # to the area of self.channels[chan][i] containing only
        # valid data; after aligning it might be the case that
        # valid data are padded by NANs in self.channels[chan][i]
        
        # optional text description of the data stored in the
        # structure
        self.description = "--No description given yet.--"
        
        if paths is not None:
            self.load_data(paths)
        
    def __iter__(self):
        """iterator over SPMdata structure is set to be iterator
        over its layers
        """
        
        return iter(self.layers)
    	
    def view_channel(self, chan, view, key=None, laylist=None):
        """put all layers in array corresponding to channel
        'chan' on top of each other, such that layers with
        smaller dimensions are seen in front of layers with
        larger dimensions
        
        chan - channel to be viewed
        view - Cut2D structure into which the view is saved
        key - key function for sorting
        laylist - list of layers to be viewed; if None, then
        all layers are used
        """
        
        if chan not in self.channels.keys():
            print(("view_channel: Invalid channel '{}'."
                ).format(chan))
            return
        
        if laylist is None: laylist = list(range(self.numlay))
        key = (lambda x: x.xnum) if key is None else key
        
        # allocate view array
        val = np.empty_like(self.channels[chan][0])
        val.fill(np.nan)

        # store data into view array
        for lay in sorted(self.layers, key=key):
            if chan not in lay.channels.keys(): continue
            val[lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum
            ] = lay.channels[chan]
        
        # save view array into Cut2D structure 'view'
        view.channels[chan] = val.copy()
        view.minmax[chan] = (np.nanmin(self.channels[chan]),
                             np.nanmax(self.channels[chan]))
    	
    def delete_layers(self, laylist):
        """delete layers specified by 'laylist'
        """
        
        # list of layers to be preserved and list of channels
        # to be deleted
        complist = [i for i in range(self.numlay)
            if i not in laylist]
        chandellist = []
        
        # leave only relevant layers
        self.layers = [lay for i, lay in enumerate(self.layers)
            if i in complist] 
               
        for chan, val in self.channels.items():
            # for each channel preserve only relevant layers
            arr = val[complist].copy()            
            
            # which layers contain channel 'chan'?
            nonmissing = [i for i, lay in enumerate(self.layers)
                if chan in lay.channels.keys()]            

            # if 'chan' is missing in all preserved layers,
            # delete 'chan'
            if not nonmissing:
                print(("delete_layers: Channel '{}' was present"
                    " only in deleted layers and will be deleted"
                    " as well.").format(chan))
                # since direct deletion would alter
                # self.channels.items() over which the iteration
                # is done, we store the channel name in chandellist
                chandellist.append(chan)
            else:
                self.setchan(chan, arr, updatelayers=True)
                
        # delete all useless channels
        for chan in chandellist:
            self.delchan(chan)
            
        print(("delete_layers: There remain {} layers and channels"
            " {}.").format(self.numlay, list(self.channels.keys())))

    def delete_channels(self, *cl):
        """delete channels specified by 'chanlist'; if there
        remain no channels in the SPMdata structure after deletion,
        layers are NOT deleted, i.e. they preserve their metadata
        """
        
        for chan in cl:
            self.delchan(chan)
        print("delete_channels: Channels '{}' deleted.".format(cl))
        
    def sort_spm(self, key=None):
        """sort layers
        
        key - sorting function; if None, then sorting according
            to heights is used
        """
        
        # if implicit sorting is used, check whether all layers
        # have defined attribute 'height'; similar discussion
        # should be made for any other 'key'
        if key is None:
            heilist = [lay.height is None for lay in self.layers]
            if any(heilist):
                print("sort_spm: Not all layers have defined"
                    " 'height' attribute. Halt.")
                return
        
        key = (lambda x: x.height) if key is None else key

        # mark each layer by its rank
        for i, lay in enumerate(self.layers):
            setattr(lay, "_rank", i)
        
        # sort layers
        self.layers.sort(key=key)
        
        # store sorting permutation
        indsort = [lay._rank for lay in self.layers]

        # delete auxiliary rank parameter
        for i, lay in enumerate(self.layers):
            delattr(lay, "_rank")        

        # reallocate arrays for all channels according
        # to sorting permutation
        for chan, val in self.channels.items():
            arr = val[indsort].copy()
            
            # update channel
            self.setchan(chan, arr, updatelayers=True)
        
    def rewrite_channel(self, channame, arr, checknan=True):
        """replace array corresponding to 'channame' by 'arr';
        all other metadata are preserved and layer pointers are
        automatically updated
        
        channame - channel to be updated
        arr - new array
        checknan - if True, then 'arr' is tested whether it
            contains valid data inside valid data frame;
            if False, then this test is not performed and channels
            are updated even by incomplete data, which may cause
            problems in aligning procedures
        """
        
        # check channel name
        chan = channame
        if chan not in self.channels.keys():
            print(("rewrite_channel: Invalid channel"
                " name '{}'.").format(chan))
            return
            
        # check whether in the valid data frames for each layer
        # there are NaNs; if so, halt with no change in SPMdata
        # structure
        for i, lay in enumerate(self.layers):
            if chan not in lay.channels.keys(): continue
            arrslice = arr[i,
                lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum]
                
            if checknan and np.any(np.isnan(arrslice.flat)):
                print("rewrite_channel: Nuexpected NaNs in a "
                    "valid data frame. Halt.")
                return

        # if no NaNs present in valid data frame, update pointers
        # in each layer and save new data            
        self.setchan(chan, arr, updatelayers=True)
        
    def check_valid_frames(self):
        """maintenance tool for analysing whether frames defined
        by lay.xoffind, lay.yoffind, lay.xnum, lay.ynum for each
        layer lay contain only valid data
        """
        
        # THIS PROCEDURE SHOULD BE AUGMENTED BY ADDITIONAL
        # FEATURES, NOW IT ONLY CHECKS WHETHER THERE ARE SOME
        # NANS IN VALID DATA FRAMES, BUT NOT CHECK WHETHER THERE
        # ARE VALID DATA OUTSIDE THE FRAME ETC.
        
        nonans = True
        for i, lay in enumerate(self.layers):
            for chan, val in lay.channels.items():
                if np.any(np.isnan(val.flat)):
                    print(("check_valid_frames: In layer no."
                        " {:3}, channel '{}' there are redundant"
                        " NaNs.").format(i, chan))
                    nonans = False
        if nonans:
            print("check_valid_frames: There are no redundant NaNs.")
        
    def add_layers(self, laylist, heilist=None):
        """insert additional layers full of NaNs into SPMdata
        structure
        
        laylist - list of layer indices to be inserted, e.g. if
            there are 5 layers initially and one wants to add one
            layer after the first layer in SPMdata and two layers
            at the end of SPMdata, one sets laylist=[2,6,7], since
            the final SPMdata looks like:
            [1(old), 2(new), 3(old), 4(old), 5(old), 6(new), 7(new)];
            it is assumed that laylist does not contain two
            identical indices; this assumption also applies to
            negative indices, i.e.
            if self.numlay == 6 and laylist = [2, -4],
            then we obtain wrong behaviour of the routine since
            ind1=2 denotes the same position as ind2 = -4;
        heilist -  list of heights for new layers;
            if None, then lay.height is left as None
        """
        
        # process list of layers
        laylist = set(laylist)
        numl = len(laylist) + self.numlay   
        laylist = [ind % numl for ind in laylist]
        laylist.sort()
        numl = len(laylist) + self.numlay        
        
        # if indices out of bounds, halt
        if max(laylist) >= numl:
            print("add_layers: Invalid list of layer indices.")
            return

        # process list of heights
        heilist = [None]*numl if heilist is None else heilist[:numl]
        heilist += [None]*(numl - len(heilist))

        # allocate new layers
        for i, ind in enumerate(laylist):
            # insert new layer in the position ind
            self.layers[ind:ind] = [SPMdataLayer()]
            
            self.layers[ind].height = heilist[i]

        # reallocate channel arrays
        for chan, val in self.channels.items():
            arr = np.empty((numl, val.shape[1], val.shape[2]))
            arr.fill(np.nan)
            
            complist = [i for i in range(numl) if i not in laylist]
            arr[complist] = val
            self.setchan(chan, arr, updatelayers=True)
        
    def add_channel(self, channame, arr, direction='forward',
        units=None, deep=True, xofflist=None, yofflist=None,
        xnum=None, ynum=None, xstep=1.0, ystep=1.0):
        """add a new channel to the structure
        
        arr - array to represent the channel, if arr corresponds
            to the number of layers greater than self.numlay,
            then arr is reallocated and padded with NaNs to match
            the self.numlay, if arr corresponds to the number of
            layers less than self.numlay, then only first
            self.numlay layers of arr is used
        channame - name of the channel, if channame is already
            present in the structure, the original array is
            overwritten
        direction - which scan direction the arr represents, it
            can be either 'forward' or 'backward';
            if the latter is used and channame does not ends with
            'b', this 'b' is appended to the channame
        units - in which units the data in arr should be
            interpreted;
            if None, then default UNKNOWN value is used
        deep - if True, then deep copy of arr is used
        xofflist, yofflist - x- and y-offinds for each layer;
            if None, then zero offinds are used
        xnum, ynum - lists of xnums and ynums for each layer;
            if None, then xnum and ynum of each layer is set to 
            arr.shape[1] and arr.shape[2], respectively; these
            parameters are used only when the SPMdata structure
            is empty before adding channame; if there are already
            some channels in SPMdata structure, then these
            parameters ARE NOT USED
        xstep, ystep - data resolution in x- and y-direction, i.e.
            how long physical dimension correspond to a single 
            index in the arrays for individual channels;
            if they are left unspecified, then these are by default
            set as 1 to the SPMdata structure
        """
        
        # process the channel name with respect to the direction
        if direction == 'backward':
            if not is_backward_scan(channame):
                print(("add_channel: Invalid channel name '{}'"
                    " for backward scanning.").format(channame))
                return
        elif direction == 'forward':
            if is_backward_scan(channame):
                print(("add_channel: Invalid channel name '{}'"
                    " for forward scanning.").format(channame))
                return
        else:
            print("add_channel: Invalid channel name.")
            return
        
        # if SPMdata is empty, set xnums and ynums for each layer
        if not self.channels:
            # initialize xnums and ynum and adjust them
            # to match number of layers
            if xnum is None:
                print("add_channel: No xnums specified.")
                xnum = [arr.shape[1]]*self.numlay
            else:
                xnum = list(xnum[:self.numlay])
                xnum += [arr.shape[1]]*(self.numlay - len(xnum))

            if ynum is None:
                print("add_channel: No ynums specified.")
                ynum = [arr.shape[2]]*self.numlay
            else:
                ynum = list(ynum[:self.numlay])
                ynum += [arr.shape[2]]*(self.numlay - len(ynum))            
            
            # update xnums and ynums for each layer
            for i, lay in enumerate(self.layers):
                lay.xnum, lay.ynum = xnum[i], ynum[i]
        
        # new channel array must be at least as big as valid data
        # frame
        xnum, ynum = max(self.xnums), max(self.ynums)
        if arr.shape[1] < xnum or arr.shape[2] < ynum:
            print(("add_channel: Input array (x, y)-shape must"
                " be at least ({}, {}). Provided shape is"
                " ({}, {}).").format(xnum, ynum,
                arr.shape[1], arr.shape[2]))          
            return      
        
        # process data array
        if arr.shape[0] >= self.numlay:
            array = arr[:self.numlay]
        else:
            array = np.empty(
                (self.numlay, arr.shape[1], arr.shape[2]))
            array.fill(np.nan)
            array[:arr.shape[0]] = arr
        array = array.copy() if deep else array
        
        # add units for new channel
        if units is None:
            self.units[channame] = UNITS[UNKNOWN]
        else:
            if units in self.units:
                self.units[channame] = UNITS[units]
            else:
                print(("add_channel: Unknown units '{}'. Default"
                    " units used instead.").format(units))
                self.units[channame] = UNITS[UNKNOWN]
        
        # initialize offinds and adjust them to match number of
        # layers
        if xofflist is None:
            print("add_channel: No index x-offsets specified.")
            xofflist = [0]*self.numlay
        else:
            xofflist = list(xofflist[:self.numlay])
            xofflist += [0]*(self.numlay - len(xofflist))
    
        if yofflist is None:
            print("add_channel: No index y-offsets specified.")
            yofflist = [0]*self.numlay
        else:
            yofflist = list(yofflist[:self.numlay])
            yofflist += [0]*(self.numlay - len(yofflist))

        # maximal dimensions
        xn = [xofflist[i] + lay.xnum
            for i, lay in enumerate(self.layers)]
        yn = [yofflist[i] + lay.ynum
            for i, lay in enumerate(self.layers)]

        # xoffind + xnum must be always <= x-shape of the array,
        # and similarly for y-shape 
        if (np.array(xn) > array.shape[1]).any() \
            or (np.array(yn) > array.shape[2]).any():
            print("add_channel: Index offsets specified do not"
                " match array shape.")
            return
        
        # if channame is already present in the structure,
        # delete it
        if channame in self.channels:
            self.delchan(channame)

        # update offinds for layers        
        for i, lay in enumerate(self.layers):
            lay.xoffind[channame] = xofflist[i]
            lay.yoffind[channame] = yofflist[i]
            
            # setchan routine needs to know all channels present
            # in the lay.channels, so we add it artificially
            lay.channels[channame] = None

        # set a channel to the structure
        self.setchan(channame, array, updatelayers=True)
        self.xsteps[channame] = xstep
        self.ysteps[channame] = ystep
        
        print(("add_channel: New channel '{}' successfully"
            " added.").format(channame))
        
    # PROPERTIES, AUXILIARY INQUIRING FUNCTIONS, GETTERS
        
    @property
    def channels_counts(self):
        """dictionary, where for each channel there is a list of
        layers containing the given channel
        
        self.channels_counts = {ELEV: [0, 1, 3, 4, 7, ...], ...}
        """
        
        return {chan: [i for i, lay in enumerate(self.layers)
            if chan in lay.channels.keys()]
            for chan in self.channels}
        
    @property
    def arrshape(self):
        """dictionary of (x, y) shapes of all arrays corresponding
        to channels in the structure
        self.arrshape = {ELEV:(xdim of 3D array for ELEV,
        ydim of 3D array for ELEV), ...}
        """
        
        return {chan: self.channels[chan].shape[1:]
            for chan in self.channels}
        
    def __len__(self):
        """length, i.e. number of layers
        """

        return len(self.layers)
        
    def print_channels(self):
        """print channels present in self
        """
        
        print(("print_channels: There are {} channels: {}."
            ).format(len(self.channels.keys()),
            list(self.channels.keys())))
                
    @property
    def numlay(self):
    	"""number of layers, the same output as for len()
    	"""
    
    	return len(self.layers)
        
    @property
    def xoffinds(self):
        """index offsets in x-direction for each channel
        """
        
        return {chan: [lay.xoffind[chan]
            for lay in self.layers if chan in lay.channels]
            for chan in self.channels.keys()}

    @property
    def yoffinds(self):
        """index offsets in y-direction for each channel
        """
        
        return {chan: [lay.yoffind[chan]
            for lay in self.layers if chan in lay.channels]
            for chan in self.channels.keys()}
    
    def print_offinds(self, chan):
        """print xoffind and yoffind for all layers for channel
        'chan'
        """

        if chan not in self.channels.keys():
            print(("print_offinds: Invalid channel '{}'."
                ).format(chan))
            return

        print("print_offinds: In total {} layer{}{}".format(
            self.numlay, "s" if self.numlay != 1 else "",
    		"." if self.numlay == 0 else (
    		" with offinds (xoffind, yoffind) for channel '{}':"
    		).format(chan)))
        for i, lay in enumerate(self.layers):
            if chan not in lay.channels.keys():
                xoff, yoff = "absent", "absent"
            else:
                xoff, yoff = lay.xoffind[chan], lay.yoffind[chan]
            print("\tLayer no. {:<3}: {}, {} ".format(
                i, xoff, yoff))
    
    @property
    def xoffsets(self):
        """physical offsets in x-direction for each layer
        """
        
        return [lay.xoffset for lay in self.layers]

    @property
    def yoffsets(self):
        """physical offsets in y-direction for each layer
        """

        return [lay.yoffset for lay in self.layers]
    
    def print_offsets(self):
        """print xoffset and yoffset for all layers
        """

        offs = set(zip(self.xoffsets, self.yoffsets))        
        if len(offs) == 1:
            print("print_offsets: In total {} layer{} with".format(
                self.numlay, "s" if self.numlay != 1 else "") + \
                " {}offsets: (xoffset, yoffset) = {} {}".format(
        		"identical " if self.numlay > 1 else "",
        		offs.pop(), UNITS[ELEV]))
        else:
            print("print_offsets: In total {} layers{}".format(
                self.numlay, "." if self.numlay == 0 else \
                " with offsets (xoffset, yoffset):"))
            for i, lay in enumerate(self.layers):
                print("\tLayer no. {:<3}: {}, {} {}".format(
                    i, lay.xoffset, lay.yoffset, UNITS[ELEV]))
    
    def print_steps(self):
        """print xsteps and ysteps for all layers
        """

        chanlen = len(self.channels.keys())
        steps = set(zip(self.xsteps, self.ysteps))        
        if len(steps) == 1:
            print("print_steps: In total {} channel{}".format(
                chanlen, "s" if chanlen != 1 else "") + \
                " with {}steps: (xstep, ystep) = {} {}".format(
        		"identical " if chanlen > 1 else "",
        		steps.pop(), UNITS[ELEV]))
        else:
            print("print_steps: In total {} channel{}".format(
                chanlen, "." if chanlen == 0 else \
                " with steps (xstep, ystep):"))
            for chan, val in self.channels.items():
                print("\tChannel '{:7}': {:.6}, {:.6} {}".format(
                    chan, self.xsteps[chan], self.ysteps[chan],
                    UNITS[ELEV]))

    @property
    def xranges(self):
        """physical range of valid data in x-direction for each
        layer
        """
        
        return [lay.xran for lay in self.layers]

    @property
    def yranges(self):
        """physical range of valid data in y-direction for each
        layer
        """

        return [lay.yran for lay in self.layers]
    
    def print_ranges(self):
        """print xrange and yrange for all layers
        """

        rans = set(zip(self.xranges, self.yranges))        
        if len(rans) == 1:
            print("print_ranges: In total {} layer{} with".format(
                self.numlay, "s" if self.numlay != 1 else "") + \
                " {}ranges: (xran, yran) = {} {}".format(
        		"identical " if self.numlay > 1 else "",
        		rans.pop(), UNITS[ELEV]))
        else:
            print("print_ranges: In total {} layers{}".format(
                self.numlay, "." if self.numlay == 0 else \
                " with ranges (xran, yran):"))
            for i, lay in enumerate(self.layers):
                print("\tLayer no. {:<3}: {}, {} {}".format(
                    i, lay.xran, lay.yran, UNITS[ELEV]))  

    @property
    def heights(self):
        """heigts of scans for each layer
        """
        
        return [lay.height for lay in self.layers]
        
    def print_heights(self):
    	"""print height for all layers
    	"""
    	
    	print("print_heights: In total {} layer{}{}".format(
    	    self.numlay, "s" if self.numlay != 1 else "",
    		"." if self.numlay == 0 else " with heights:"))
    	for i, lay in enumerate(self.layers):
    	    print("\tLayer no. {:<3}: {} {}".format(
    	        i, lay.height, UNITS[ELEV]))
	
    @property
    def filenames(self):
        """names of files from which data were taken for each
        layer
        """

        return [lay.filename for lay in self.layers]
	
    def print_filenames(self):   
    	"""print filename for all layers
    	"""
    	
    	print("print_filenames: In total {} layer{}{}".format(
    	    self.numlay, "s" if self.numlay != 1 else "",
    		"." if self.numlay == 0 else " with filenames:"))
    	for i, lay in enumerate(self.layers):
    	    print("\tLayer no. {:<3}: {}".format(i, lay.filename))
    	    
    @property
    def nums(self):
        """dimensions of valid data
        """
        
        return [(lay.xnum, lay.ynum) for lay in self.layers]

    @property
    def xnums(self):
        """x-dimensions of valid data
        """
        
        return [lay.xnum for lay in self.layers]

    @property
    def ynums(self):
        """y-dimensions of valid data
        """
        
        return [lay.ynum for lay in self.layers]

    def print_nums(self):
        """print dimensions of valid data
        """

        dimsl = set(self.nums)
        if len(dimsl) == 1:
            print("print_nums: In total {} layer{} with".format(
                self.numlay, "s" if self.numlay != 1 else "") + \
                " {}dimensions: (xnum, ynum) = {}".format(
        		"identical " if self.numlay > 1 else "",
        		dimsl.pop()))
        else:
            print("print_nums: In total {} layers{}".format(
                self.numlay, "." if self.numlay == 0 else \
                " with dimensions (xnum, ynum):"))
            for i, lay in enumerate(self.layers):
                print("\tLayer no. {:<3}: {}, {}".format(
                    i, lay.xnum, lay.ynum))        

    # INDEXING

    def ind_from_height(self, height):
        """take height of the scan and return its corresponding
        or most approximate layer index
        """
        
        try:
            height = float(height)
        except ValueError:
            print("ind_from_height: Invalid height specification.")
            return None

        diffs = np.array(list(map(lambda x: abs(x.height - height),
            self.layers)))
        if diffs.min() > 1: # 1 or any other comparatively big number
            print("ind_from_height: Input height quite far from "
                "stored data.")

        return np.argmin(diffs)

    def __getitem__(self, item):
        """if the item is a string, it is assumed to be a height
        of some layer, otherwise the item is assumed to be a normal
        index
        """
        
        if isinstance(item, str):
            itemloc = self.ind_from_height(item)
            return self.layers[itemloc]
        else:
            return self.layers[item]

    # DATA LOADING

    def load_data(self, paths=(), suffixes=SUFFIXES, key=None,
        finer=True):
        """loads data from external files

        paths - name of the file or tuple of filenames storing
            the data;
            if 'paths' is evaluated to False, then a graphical
            window is shown for choosing files to be loaded;
            if 'paths'='cwd', then current working directory is
            used to exctract files from;
            if 'paths' contains a string with file name, then
            this file is loaded;
            if 'paths' contains a string with a directory name,
            then all valid files from that directory are loaded;
            if 'paths' is a sequence of strings containing file
            names, then all valid files with these file names are
            loaded
        suffixes - relevant format suffixes
        key - function evaluated for each layer, according to
        	outputs of this function the loaded layers are
        	embedded into 3D arrays, if None, then sorting
        	according to heights is used
        finer - if data resolution for individual layers differs,
            two situations may happen:
            if finer=True, then data with lower resolution are
            upsampled so that they can be added to data with higher
            resolution;
            if finer=False, then data with higher resolution are
            downsampled so that they can be added to data with lower
            resolution;
            each layer contains attributes horstep and verstep,
            which represent data resolution in x- and y-direction;
            these attributes should not be modified, since they
            express resolution of the original data
        """

        spm_load_data(self, paths=paths, suffixes=suffixes,
            key=key, finer=finer)

    # DATA MANIPULATION
    
    def setchan(self, chan, val, updatelayers=False):
        """set channel to the structure
        
        chan - channel to be set
        val - value to be set to channel 'chan'
        updatelayers - iff True, then pointers to valid data
            frame are updated for each layer    
        """
        
        setattr(self, chan, val)
        self.channels[chan] = val
        
        if updatelayers:
            for i, lay in enumerate(self.layers):
                if chan not in lay.channels.keys(): continue
                lay.setchan(chan, self.channels[chan][i,
                    lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                    lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum
                    ])

        
    def delchan(self, chan):
        """delete channel from the structure
        
        chan - channel to be deleted
        """
        
        if chan not in self.channels:
            return
        
        del self.channels[chan]    
        for lay in self.layers:
            lay.delchan(chan)
            
        delattr(self, chan)
    	    	
    def average_channels(self, cl=None, fun=None, spmout=True):
        """average forward and backward scans for given channels
        
        THIS PROCEDURE SHOULD BE USED PRIOR TO ANY ALIGN PROCEDURE!!!
        
        cl - list of channels to be averaged; if None, then all
            channels are averaged
        fun - average function, if None, then arithmetic average
            is used
        spmout - if True, then averaged channels will be put into
            SPMdata structure and returned;
            if False this routine returns a dictionary of the form
            '{chan:arrchan, ...}' for 'chan' in 'cl' and
            corresponding array 'arrchan'
            
        returns:
            if 'spmout'==True, then SPMdata structure is returned,
            containing average channels as channels;
            if 'spmout'==False, then a dictionary of channels is
            returned
        """

        # choose averaging function
        if fun is None:
            factor = 0.5
            
            # already assumes that y is reverted:
            avfun = lambda x, y: factor*x + (1 - factor)*y
        else:
            avfun = fun
        
        # which channels to average:
        # if cl is None, then take all channels
        cl = self.channels.keys() if cl is None else cl
        cl = [chan for chan in cl if chan in self.channels.keys()]
        # take only such forward channels, for which there is
        # corresponding backward channel
        cl = [chan for chan in cl 
                if not is_backward_scan(chan)
                if get_backward_scan(chan) in self.channels.keys()]
        
        res = {}    
        
        # average channels...
        for chan in clf:
            forscan = self.channels[chan]
            bckscan = self.channels[get_backward_scan(chan)]
            res[chan] = avfun(forscan, bckscan)
            print(("average_channels: Channel '{}' averaged."
                ).format(chan))

        if spmout:
            spmdata = create_spm()
            spmdata.channels = res
            spmdata.determine()
            return spmdata        
        else:
            return res    

    def determine(self, hew=True):
        """this routine should be used wisely, since it rewrites
        all original values of offinds, nums and pointers to
        valid data frame for all channels in all layers
        
        this routine takes self.channels and based on data stored
        in each array for each channel determines xnums, ynums,
        xoffinds and yoffinds for each layer and updates pointer
        to valid data frames accordingly
        
        hew - if True, then xnums and ynums are determined to
            define valid data frame, which contains really only
            a valid data for each channel;
            if False, then some channels might have NaNs in a 
            valid data frame for some channel;
            if no valid data are found in a given channel for
            a given layer, then the channel is not saved into
            the layer at all
        """
        
        xnums, ynums = [], []
        xoffs, yoffs = [], []

        # for each layer...
        for i, lay in enumerate(self.layers):

            xnumslay, ynumslay = {}, {}
            xoffslay, yoffslay = {}, {}
            
            # for each channel...
            for chan, val in self.channels.items():

                # find out where are NaNs, elementwise
                mask = np.isnan(val[i])
                
                if hew:
                    # if we want no NaNs whatsoever to
                    # appear in valid data frame,
                    # then search for any NaNs in rows
                    # and columns
                    maskrow = mask.all(axis=1)
                    maskrow = mask[~maskrow]
                    maskclm = maskrow.any(axis=0)                
                    maskrow = mask[:,~maskclm]
                    maskrow = maskrow.any(axis=1)
                else:
                    # otherwise search for rows and columns
                    # containing nothing but NaNs
                    maskrow = mask.all(axis=1)
                    maskclm = mask.all(axis=0)
        
                # np.argmin returns index of the first occurence
                # of the minimum value for given array
                minrow = np.argmin(maskrow)
                minclm = np.argmin(maskclm)

                # the last occurence is taken as the first 
                # occurence in the reversed array, in which
                # case we have to add the array dimension
                maxrow, maxclm = maskrow.shape[0], maskclm.shape[0]
                maxrow -= np.argmin(maskrow[::-1])
                maxclm -= np.argmin(maskclm[::-1])
                        
                # if there are no valid data for channel 'chan'
                # in the present layer, delete the channel from
                # the layer
                if mask.all() and chan in lay.channels.keys():
                    # condition something like:
                    # if maxrow <= minrow or maxclm <= minclm: ...
                    # cannot be used, because if the array is full
                    # of NaNs, then mask is full of True's and
                    # argmin returns the first occurence of True,
                    # instead of the first occurence of False
                    print(("determine: Channel '{}' in layer no."
                        " {} deleted.").format(chan, i))
                    lay.delchan(chan)
                else:
                    # store xnum, ynum, xoffind, yoffind for each
                    # channel and given layer
                    xnumslay[chan] = (maxrow - minrow)
                    ynumslay[chan] = (maxclm - minclm)
                    xoffslay[chan] = minrow
                    yoffslay[chan] = minclm
        
            xnums.append(xnumslay)
            ynums.append(ynumslay)

            # update xoffinds, yoffind for given layer
            lay.xoffind = xoffslay
            lay.yoffind = yoffslay
        
        # calculate such xnum, ynum for each layer, so that it
        # defines valid data frame for all channels simultaneously
        xnuml = [min(lay.values()) if lay.values() else 0
            for lay in xnums]
        ynuml = [min(lay.values()) if lay.values() else 0
            for lay in ynums]
        
        # update xnum, ynum and pointers for each layer
        for i, lay in enumerate(self.layers):        
            lay.xnum, lay.xnum = xnuml[i], ynuml[i]
            
            # here no self.setchan is used since we change only
            # layer properties, not the arrays
            for chan in lay.channels.keys():
                xoff, yoff = lay.xoffind[chan], lay.yoffind[chan]
                lay.setchan(chan, self.channels[chan][i, 
                    xoff:xoff + lay.xnum,
                    yoff:yoff + lay.ynum])

    def strip_nans(self, *cl, same=True):
        """strip NaNs from arrays corresponding to channels in
        'cl', so that all valid data are preserved
        
        cl - sequence of channels, if empty, then all channels
            are stripped
        same - if True, then all channel arrays are stripped to
            have the same dimensions, this option is useful for
            already aligned arrays;
            if False, then each channel is stripped separately
        """
        
        if not cl:
            cl = list(self.channels.keys())
        else:
            cl = [chan for chan in cl
                if chan in self.channels.keys()]
            
        xmin, ymin = [None]*len(cl), [None]*len(cl)
        xmax, ymax = [None]*len(cl), [None]*len(cl)
        
        # find frame limits
        for i, chan in enumerate(cl):
            # take only layers containing channel 'chan'
            laylist = [lay for lay in self.layers
                if chan in lay.channels.keys()]

            # from where to take valid data
            xmin[i] = min([lay.xoffind[chan] for lay in laylist])    
            ymin[i] = min([lay.yoffind[chan] for lay in laylist])    
            
            # up to where to take valid data
            xmax[i] = max([lay.xoffind[chan] + lay.xnum
                for lay in laylist])
            ymax[i] = max([lay.yoffind[chan] + lay.ynum
                for lay in laylist])

        if same:            
            # if all arrays have to have the same dimensions...
            xminv, yminv = min(xmin), min(ymin)
            xmaxv, ymaxv = max(xmax), max(ymax)

            for i, chan in enumerate(cl):
                # create new array
                aux = np.empty(
                    (self.numlay, xmaxv - xminv, ymaxv - yminv),
                    dtype=float)
                aux.fill(np.nan)                
                aux[:,
                    (xmin[i] - xminv):(xmax[i] - xminv),
                    (ymin[i] - yminv):(ymax[i] - yminv)
                    ] = self.channels[chan][:,
                    xmin[i]:xmax[i], ymin[i]:ymax[i]]
                  
                # update xoffinds and yoffinds  
                for j, lay in enumerate(self.layers):
                    if chan not in lay.channels.keys(): continue
                    lay.xoffind[chan] -= xminv
                    lay.yoffind[chan] -= yminv
                    
                # update channel
                self.setchan(chan, aux.copy(), updatelayers=True)            
        else:            
            # if channel arrays do not need to have the sam
            # dimensions...
            for i, chan in enumerate(cl):
                # create new array
                aux = self.channels[chan][:,
                    xmin[i]:xmax[i], ymin[i]:ymax[i]]            

                # update xoffinds and yoffinds  
                for j, lay in enumerate(self.layers):
                    if chan not in lay.channels.keys(): continue
                    lay.xoffind[chan] -= xmin[i]
                    lay.yoffind[chan] -= ymin[i]
                    
                # update channel
                self.setchan(chan, aux.copy(), updatelayers=True)            

    def hew_nans(self, minl=None, maxl=None, same=True):
        """hew_nans NANs from 3D arrays corresponding to each
        channel, the hew_nans is performed for layers whose
        indices lie between minl and maxl, all remaining layers
        are deleted, resulting 3D arrays contain subset of purely
        valid data
        
        minl, maxl - cut off is performed for layers:
            minl, minl + 1, ..., maxl - 1, maxl
        same - if True, then all channel arrays are hewed to have
            the same dimensions, this option is useful for already
            aligned arrays; if False, then each channel is hewed
            separately
        """
        
        # range of layers for hew_nans
        lmin = 0 if minl is None else minl
        lmax = self.numlay if maxl is None else maxl
        lmin = lmin % (self.numlay + 1)
        lmax = lmax % (self.numlay + 1)
        if lmax < lmin: lmax, lmin = lmin, lmax

        print(("hew_nans: Cutoff performed for layers {} through"
            " {}.").format(lmin, lmax))
            
        # delete outside layers
        self.layers = self.layers[lmin:lmax].copy()
        chanlen = len(self.channels)
        
        xmin, ymin = [None]*chanlen, [None]*chanlen
        xmax, ymax = [None]*chanlen, [None]*chanlen

        # find frame limits
        for i, chan in enumerate(self.channels):
            # take only layers containing channel 'chan'
            laylist = [lay for lay in self.layers
                if chan in lay.channels.keys()]

            # from where to take valid data
            xmin[i] = max([lay.xoffind[chan] for lay in laylist])
            ymin[i] = max([lay.yoffind[chan] for lay in laylist])    
            
            # up to where to take valid data
            xmax[i] = min([lay.xoffind[chan] + (lay.xnum - 1)
                for lay in laylist])
            ymax[i] = min([lay.yoffind[chan] + (lay.ynum - 1)
                for lay in laylist])
            
        if same:
            # if all arrays have to have the same dimensions...
            xminv, yminv = max(xmin), max(ymin)
            xmaxv, ymaxv = min(xmax), min(ymax)

            xmin, xmax = [xminv]*chanlen, [xmaxv]*chanlen
            ymin, ymax = [yminv]*chanlen, [ymaxv]*chanlen

        # update channels
        for i, chan in enumerate(self.channels):      
            # take valid data
            aux = self.channels[chan][lmin:lmax,
                xmin[i]:xmax[i], ymin[i]:ymax[i]]            
            
            # update layer properties
            for j, lay in enumerate(self.layers):
                if chan not in lay.channels.keys(): continue     
                
                # layer offinds
                lay.xoffind[chan], lay.yoffind[chan] = 0, 0
                
                # since each layer does not have to have all
                # channels, we set xnum, ynum only when 
                # chan is present in the given layer; due to
                # being in the loop xnum and ynum are set
                # as many times as there are channels in the
                # layer, but whatever
                lay.xnum, lay.ynum = aux.shape[1:]
            
            # update channel
            self.setchan(chan, aux.copy(), updatelayers=True)


    # ALIGNMENT

    def align_correlate(self, *cl, rchan=None, rlay=0, alay=None,
        ilay=None, actll=None, fitll=None, misll=None,
        show=False, xcoefs=None, ycoefs=None, fac=1,
        downsample=True, order=1):
        """align arrays according to rchan channel, i.e. offsets
        are calculated only once for rchan and then used multiple
        times for all channels in cl, unknown values are
        represented by NANs, resulting array may have bigger
        dimensions than original arrays
        
        cl - list of channels to be aligned, if empty, then all
            channels are aligned
        rlay - the index of a fixed reference layer with which
            all layers are compared
        alay - if actll is not None or not all layers, then alay
            specifies to which layer inactive layers should adhere,
            i.e. offind of ilay layer will be set to the offind
            of alay;
            if alay is None or lies outside actll, then it is set
                to be the last layer in actll
        ilay - analog of alay for inactive layers, if None or in 
            actll, then it is set to be the first layer not in
            actll
        rchan - all channels in cl are aligned according to
            the reference channel rchan
        actll - list of 'active' layers, which are to be aligned;
            if None, then all layers are used
        fitll - list of 'fitting' layers which are used to
            determine offsets for layers in misll;
            if None or empty, then fitll is effectively put equal
                to all layers not lying in misll
        misll - list of 'missing' layers for which offsets should
            be calculated by inter-/extra-polation;
            if None or empty, then correlation calculation is
                performed for all layers and no
                inter-/extra-polation is used;
            if misll and fitll are not disjoint, then misll is
                put equal to misll minus fitll;
            if misll='all', then all layers are used for 
                inter-/extra-polation
        show - if True, then plot depicting fitll points and
            extrapolation trends is shown,
            control points - fitll points
            calculated points - misll points
            x fit, y fit - fit of x-offsets and y-offsets
            (manual) - shown if fit is set manually via xcoefs
                or ycoefs
            x offsets, y offsets - resulting x- and y-offsets
            red area - for this area no alignment is active,
                i.e. corresponding layers do not lay in actll
        xcoefs, ycoefs - linear and absolute coefficients for manual
            interpolation for x and y dimension, respectively;
            if not None, then these parameters are used instead of
                automatically calculated ones, in which case 
                xcoefs=(xlincoef, xabscoef) etc.;
            if not None, then fitll is not used
        fac - upsampling factor
        downsample - whether to perform downsampling after
            correlation
        order - order of interpolation for resampling
            
        FOR EXAMPLES SEE ALIGN_CORRELATE IN MODULE SPM_DATA_ALIGN
        """
       
        spm_align_correlate(self, *cl, rchan=rchan, rlay=rlay,
            alay=alay, ilay=ilay, actll=actll, fitll=fitll,
            misll=misll, show=show, xcoefs=xcoefs, ycoefs=ycoefs,
            fac=fac, downsample=downsample, order=order)
       
    def align_offsets(self, *cl):
        """align arrays corresponding to each channel according
        their offsets, unknown values are represented by NANs, 
        resulting array may have bigger dimensions than original
        arrays, it is assumed that xranges and yranges of all
        layers are identical
        
        cl - list of channels to be aligned, if None, then all
            channels are aligned
        """

        spm_align_offsets(self, *cl)

    def align_forback(self, *cl, rchan=None, same=True, rlay=0,
        xdir=True, ydir=True, order=0):
        """align all channels in 'cl' according to 'rchan' channel
        so that forward and backward scans are moved towards each
        other and aligned 'to the centre', this routine does not
        align layers within one channel, but channels among
        themselves instead
        
        cl - sequence of channels to be aligned; if empty, then
            all channels are aligned
        rchan - reference channel; if None, then each channel is
            aligned independently
        same - if True, then resulting 3D arrays for all channels
            in cl are of the same dimensions;
            if False, then 3D arrays for forward and backward
            direction of each physical quantity have the same
            dimension
        rlay - reference layer 
        xdir - if True, than allow shift in x-direction
        ydir - if True, than allow shift in y-direction
        order - order of interpolation during shifting
        """

        spm_align_forback(self, *cl, rchan=rchan, same=same,
            rlay=rlay, xdir=xdir, ydir=ydir, order=order)
           
    # DISPLAY PROCEDURES
    
    def show_layer(self, *cl, rlay=0, rangelayer='each'):
        """show valid data for channel 'chan' in layer 'rlay'
        
        cl - sequence of channels to be shown
        rlay - which layer to show
        rangelayer - which layer to choose for scaling; if 'each',
            then scaling is based on layer 'rlay'; if 'all', then
            scaling is based on maximum and minimum of the whole
            3D array pertaining to channels in 'cl'
        """
        
        # take only those channels, which are present in 'rlay'
        cl = [chan for chan in cl
            if chan in self.layers[rlay].channels.keys()]
        if not cl:
            print("show_layer: No valid channel to show.")
            return
        
        # create plot
        fig, axes = plt.subplots(1, len(cl), squeeze=False,
            num="Valid data for layer no. {}".format(rlay))
        
        lay = self.layers[rlay]
        
        # initialize each subplot
        for ax, chan in zip(axes.flat, cl):
            # determine scaling for imshow
            if rangelayer == 'each':
                vmin, vmax = None, None
            elif rangelayer == 'all':
                vmin = np.nanmin(self.channels[chan])
                vmax = np.nanmax(self.channels[chan])
            else:
                pass
#                vmin = np.nanmin(self.channels[chan][rangelayer])
#                vmax = np.nanmax(self.channels[chan][rangelayer])

            # determine dimensions of the image
            bottolim = self.xsteps[chan] #lay.xstep
            bottolim *= lay.channels[chan].shape[0]                
            rightlim = self.ysteps[chan] #lay.ystep
            rightlim *= lay.channels[chan].shape[1]
            extent = (0, rightlim, bottolim, 0)
            
            # create image
            ax.imshow(self.layers[rlay].channels[chan],
                vmin=vmin, vmax=vmax, origin='lower',
                extent=extent)
            ax.set_title("Channel '{}' [units: {}]".format(
                chan, self.units[chan]))
        
        plt.show()
            
    def show_channel(self, chan, laylist=None, limit=8, rowlimit=4,
        aspect='equal', uniscale='present', interpolation='none',
        origin='lower'):
        """plot channel for various layers
        
        chan - which channel to plot
        laylist - which layers to plot, sequence containing indices,
            default is first 'limit' layers
        limit - maximum number of layers plotted
        rowlimit - maximum number of plots per row
        aspect - aspect of imshow
        uniscale - how to scale plots, for 'all' plots are rescaled
            to global maxima and minima of the channel throughout
            all layers, for None each plot is rescaled individually,
        	otherwise it is assumed that uniscale is a list of
        	layers over which the maxima and minima are computed
        	and used for plot scaling, for 'present' the same list
        	as 'laylist' is used
        interpolation - how to interpolate values in imshow
        origin - either 'lower' or 'upper'; determines where imshow
            should put origin of the image
        """

        if chan not in self.channels:
            print(("show_channel: Invalid channel"
                " '{}'.").format(chan))
            return

        # only first 'limit' (or less) entries plotted, layloc
        # stores layer indices
        if laylist:
            layloc = [i for i in set(laylist) if i < self.numlay]
            layloc = sorted(layloc)[:limit]
        if not laylist or not layloc:
            layloc = range(min(limit, self.numlay))
        length = len(layloc)

        # determine number of rows and columns
        if length < rowlimit:
            nrows, ncols = 1, length
        else:
            nrows = length // rowlimit
            if length % rowlimit > 0: nrows += 1
            ncols = rowlimit

        # create figure and axes
        figtit = 'Channel: ' + chan
        figtit += " [units: " + self.units[chan] + "]"
        fig, ax = plt.subplots(nrows, ncols, squeeze=False,
            num=figtit)

	    # how to scale plots
        if uniscale is None:
            vmin, vmax = None, None	
        elif uniscale == 'all':
            vmin = np.nanmin(self.channels[chan])
            vmax = np.nanmax(self.channels[chan])
        elif uniscale == 'present':
            vmin = np.nanmin(self.channels[chan][laylist])
            vmax = np.nanmax(self.channels[chan][laylist])
        else:
            vmin = np.nanmin(self.channels[chan][uniscale])
            vmax = np.nanmax(self.channels[chan][uniscale])
            
        # create plots
        l = 0
        for k in layloc:
            ind = tuple(divmod(l, rowlimit))
            try:
            	# create image
                img = ax[ind].imshow(self.channels[chan][k],
                    aspect=aspect, vmin=vmin, vmax=vmax,
                    interpolation=interpolation, origin=origin)
                ax[ind].set_title((chan + \
                    " - layer no. " + str(k % self.numlay)))
                
                # add colorbar
                divider = make_axes_locatable(ax[ind])
                cax = divider.append_axes("right", size=0.1,
                    pad=(0.5 if aspect == 'auto' else 0.1))
                plt.colorbar(img, cax=cax)
                
                # increment counter
                l += 1
            except KeyError as err:
                print(("show_channel: Key"
                    " {} not known.").format(err))
                continue

        # improve layout and launch the interactive window
        fig.tight_layout()	
        plt.show()

    def inspect_channels(self, *cl, cutlist=None, cutzlist=None,
        layer=0, xind=None, yind=None, rangelayer='all',
        linestyle='-', linecolor='w', pointcolor='g', aspect='auto',
        size="30%", interpolation='none', scaling=True): 
        """show one or more channels together with its cuts for
        arbitrary layer

        spmdatalist - list of SPMData structures used for plotting
        chanlist - list of channels to be drawn; for each index
            'i' channel chanlist[i] is shown for SPMData structure
            spmdatalist[i]
        cutlist - list into which cut control points are stored;
            if None, then nothing is stored, otherwise after each
            selection is made by mouse and saved by ctrl+p new
            triple [zpos, xpos, ypos] is appended to cutlist so
            that these data are accessible even after the
            termination of inspect_channels routine
        cutzlist - list into which z-cut points are stored;
            if None, then nothing is stored, otherwise after each
            selection is made by mouse right-button new double
            [xpos, ypos] is appended to cutzlist so that these
            data are accesssible even after the termination of
            inspect_channels routine
        layer - which layer of all the channels in chanlist to
            draw at the beginning, layer is identical for all
            channels in chanlist
        xind, yind - initial position of the intersection point
        linestyle - format specification of the intersecting lines
        linecolor - color specification of the intersecting lines
        pointcolor - color specification of the z-cut points and
            corresponding lines
        rangelayer - which layer to choose for image min and max
            determination;
            if 'all', then global minima and maxima are used;
            if 'each', then maxima and minima are determined
            independently for each layer;
            if rangelayer is not a valid number of a layer,
            then 'layer' is used instead
        aspect - aspect of the cuts in the figure
        size - relative size of cuts in the figure
        interpolation - how to interpolate values in imshow
        scaling - if True, then image labels are chosen
            to correspond to real physical dimensions of the data;
            if False, then array indices are used instead
        """
        
        spm_inspect_channels([self]*len(cl), cl, cutlist=cutlist,
            cutzlist=cutzlist, layer=layer, xind=xind, yind=yind, 
            rangelayer=rangelayer, linestyle=linestyle,
            linecolor=linecolor, pointcolor=pointcolor,
            aspect=aspect, size=size, interpolation=interpolation,
            scaling=scaling)
        
    # IMPORT, EXPORT

    def copy(self, *cl, ll=None):
        """create a copy of spmdata containing only channels
        in 'cl' and only layers in 'll'
        """
        
        # process lists of channels and layers
        cl = self.channels if not cl else cl
        cl = [chan for chan in cl if chan in self.channels]
        ll = list(range(self.numlay)) if ll is None else ll
        ll = [lay for lay in ll if 0 <= lay < self.numlay]
        
        print("copy: Channels to be copied: {}".format(cl))
        print("copy: Layers to be copied: {}".format(ll))
        
        # allocate new SPMdata structure
        copyspm = SPMdata()
        
        # copy auxiliary attributes
        itemlist = ["ysteps", "xsteps", "units"]
        for item in itemlist:
            setattr(copyspm, item, getattr(self, item).copy())
        
        # copy channels
        for chan, val in self.channels.items():
            if chan not in cl: continue
            # set new channel to copyspm
            copyspm.setchan(chan, val[ll].copy())
        
        # copy layers
        k = 0
        for i, lay in enumerate(self.layers):
            if i not in ll: continue
            copyspm.layers.append(lay.copy(k, copyspm.channels))
            k += 1
        
        # return new SPMdata structure
        return copyspm

    def __copy__(self):
        self.copy(())

    def export_data(self, filename):
        """export data to external binary file
        """
        
        with open(filename, "wb") as f:
            # TODO
            pass
            
        print("export_data: Data exported to file: {}".format(filename))
        
    def import_data(self, filename, suffixes=SUFFIXES):
        """import data from external binary file
        """
        
        validname = self._check_filenames(filenames, suffixes)
        with open(validname, "rb") as f:
            buff = f.read().split(b'\x1a\x04') # sequence \x1a\x04 taken from nanonis file
            
        # process header
        buff[0] # header
        
        self.xnum, self.ynum = [], [] # TODO
        numchan = [] # TODO
        cl = [] # LIST of relevant channels sorted by their respective appearance in data list below
        
        # process binary data
        ENDIAN = '>'
        DATA_TYPE = 'f'

        datalist = struct.unpack(ENDIAN + DATA_TYPE*self.xnum*self.ynum*numchan, buff[-1])
        self.data = np.array(dataList)
        self.data.resize((numchan, self.numlay, self.xnum, self.ynum))

#        for i, chan in enumerate(cl):
#            setattr(self, chan, data[i])

        # 
        
        # TODO
       

