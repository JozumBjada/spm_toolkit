from spm_data_layer import *
from spm_data_inspect import inspect_channels as spm_inspect_channels

import tkinter, os, math, pickle
from mpl_toolkits.axes_grid1 import make_axes_locatable    
from scipy.ndimage import map_coordinates
from scipy import signal


# TODO: BUDE POTREBA PORESIT PRIPADY, KDY NE VSECHNY VRSTVY MAJI VSECHNY KANALY, RELEVANTNI TO JE NAPR. PRI NASKLADANI DVOU SPMData NA SEBE
# JE POTREBA TAKY PORESIT NASTAVOVANI channels_counts V PROCEDURE setchan, COZ SE DOSUD DELA DOST NA HRUBO

# PORESIT, ZDA BY SE channels_counts NEDALO UDALAT JAKO PROPERTY



def stack_spm(spm1, spm2):
    """stack two SPMdata objects into a single SPMdata object 'spmnew' and return it
    
    spm1, spm2 - input SPMdata structures
    """
     
    # channels to be in spmnew and number of layers in spmnew
    chanlist = list(set(list(spm1.channels.keys()) + list(spm2.channels.keys())))
    numlay = spm1.numlay + spm2.numlay
    
    # create spmnew
    spmnew = SPMdata()
    spmnew.layers = [None]*numlay
    
    # add all channels from chanlist into spmnew
    for chan in chanlist:
    
        # get dimensions of new array provided chan is in spm1
        if chan in spm1.channels:
            arrsh1x, arrsh1y = spm1.arrshape[chan]
        else:
            arrsh1x, arrsh1y = 0, 0
            print("stack_spm: The first SPMdata structure does not contain channel"
                " '{}'".format(chan))

        # get dimensions of new array provided chan is in spm2
        if chan in spm2.channels:
            arrsh2x, arrsh2y = spm2.arrshape[chan]
        else:
            arrsh2x, arrsh2y = 0, 0
            print("stack_spm: The second SPMdata structure does not contain channel"
                " '{}'".format(chan))
            
        # allocate new array with appropriate dimensions
        arr = np.empty((numlay, max(arrsh1x, arrsh2x), max(arrsh1y, arrsh2y)))
        arr.fill(np.nan)
        
        # fill new array with valid data from spm1
        if chan in spm1.channels:
            arr[:len(spm1), :arrsh1x, :arrsh1y] = spm1.channels[chan]
            
        # fill new array with valid data from spm2
        if chan in spm2.channels:
            arr[len(spm1):, :arrsh2x, :arrsh2y] = spm2.channels[chan]
        
        # save new array into spmnew
        spmnew.setchan(chan, arr, laylist=[])
        
    # update data for layers from spm1
    for i, lay in enumerate(spm1.layers):
        spmnew.layers[i] = lay.copy(i, spmnew.channels)
        for chan in lay.channels.keys():
            spmnew.channels_counts.setdefault(chan, [])
            spmnew.channels_counts[chan].append(i)

    # update data for layers from spm2
    for i, lay in enumerate(spm2.layers):
        spmnew.layers[i + len(spm1)] = lay.copy(i + len(spm1), spmnew.channels)
        for chan in lay.channels.keys():
            spmnew.channels_counts.setdefault(chan, [])
            spmnew.channels_counts[chan].append(i + len(spm1))
    
#    print("spmnew.channels_counts: ", spmnew.channels_counts)
    
    print("stack_spm: New SPMdata object has {} layers and channels: {}".format(
        spmnew.numlay, list(spmnew.channels.keys())))
    
    return spmnew

def create_spm(paths=None, numlay=0):
    """create SPMdata structure with all necessary attributes
    
    paths - if None, then new SPMdata structure is generated with no data; otherwise
        'paths' is passed to 'load_data' routine
    numlay - number of layers in the new structure; this parameter applies only
        when 'paths' is None
    """
    
    spmdata = SPMdata(paths)
    if paths is None:
        if numlay < 0:
            print("create_spm: Invalid number of layers: {}.".format(numlay))
            del spmdata
            return
        spmdata.layers = [None]*numlay
        for i in range(numlay):
            spmdata.layers[i] = SPMdataLayer()
        # BEWARE: spmdata.layers = [SPMdataLayer()]*numlay cannot be used
        # since there is only ONE SPMdataLayer() allocated and spmdata.layers
        # then is a list of numlay  pointer to the same SPMdataLayer() object
    
    return spmdata

def pickle_data(obj, filename):
    """export data to external file by mechanism provided by pickle
    """
    
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)


def unpickle_data(obj, filename):
    """import data from external file by mechanism provided by pickle
    
    obj - list the first element of which will be assigned data from file
    """

    with open(filename, 'rb') as f:
        obj[0] = pickle.load(f)

class SPMdata:
    """class embracing all the relevant data for SPM images
    """
    
    def __init__(self, paths=None):
    
        # self.channels = {ELEV:self.elev (any scan), ...}
        self.channels = {}
        # self.channels_fb = {ELEV:(self.elev (forward scan), self.elevb (backward scan), False (is averaged)), ...}
        self.channels_fb = {}
        # self.channels_counts = {ELEV: [0, 1, 3, 4, 7, ...] (list of layers containing the given channel), ...}
        self.channels_counts = {}
        # self.layers = [<SPMdataLayer object 1>, ...] (list of layers)
        self.layers = []
        
        # sequences of channel names and units; it is possible that each SPMdata structure has its
        # own ad hoc channels with their own names and units
        self.units = UNITS
        
        # BEWARE: self.layers[i].channels[chan] is a view of self.channels[chan][i, xoffind:xoffind + num, yoffind:yoffind + ynum];
        # in other words, self.layers[i].channels[chan] points to the area of self.channels[chan][i] containing only valid data;
        # in fact, after aligning it might be the case that valid data are padded by NANs in self.channels[chan][i]
        
        if paths is not None:
            self.load_data(paths)
        
    def __iter__(self):
        """iterator over SPMdata structure is set to be iterator over its layers
        """
        
        return iter(self.layers)
    	
    def view_channel(self, chan, view, key=None, laylist=None):
        """put all layers in array corresponding to channel 'chan' on top of each other,
        such that layers with smaller dimensions are seen in front of layers with larger
        dimensions
        
        chan - channel to be viewed
        view - Cut2D structure into which the view is saved
        key - key function for sorting
        laylist - list of layers to be viewed; if None, then all layers are used
        """
        
        if chan not in self.channels.keys():
            print("view_channel: Invalid channel '{}'.".format(chan))
            return
        
        if laylist is None: laylist = list(range(self.numlay))
        key = (lambda x: x.xnum) if key is None else key
        
        # allocate view array
        val = np.empty_like(self.channels[chan][0])
        val.fill(np.nan)

        # store data into view array
        for lay in sorted(self.layers, key=key):
            if lay not in self.channels_counts[chan]: continue
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
        
        # list of layers to be preserved and list of channels to be deleted
        complist = [i for i in range(self.numlay) if i not in laylist]
        chandellist = []
        
        # leave only relevant layers
        self.layers = [lay for i, lay in enumerate(self.layers) if i in complist] 
               
        for chan, val in self.channels.items():
            # for each channel preserve only relevant layers
            arr = val[complist].copy()            
            
            # update data for each layer
            nonmissing = []
            for i, lay in enumerate(self.layers):
                # if chan is not present in the given layer, skip the layer
                if complist[i] not in self.channels_counts[chan]:
                    continue
                    
                # otherwise update the data pointer                
                xoff, yoff = lay.xoffind[chan], lay.yoffind[chan]
                xnum, ynum = lay.xnum, lay.ynum
                lay.setchan(chan, arr[i, xoff:xoff + xnum, yoff:yoff + ynum])
                nonmissing.append(i)
                
            # if chan is missing in all preserved layers, delete chan
            if len(nonmissing) == 0:
                print("delete_layers: Channel '{}' ".format(chan) + \
                    "was present only in deleted layers and will be deleted as well.")
                # since direct deletion would alter self.channels.items() over which
                # the iteration is done, we store the channel name in chandellist
                chandellist.append(chan)
            else:
                self.setchan(chan, arr, laylist=nonmissing)
                
        # delete all useless channels
        for chan in chandellist:
            self.delchan(chan)
            
        print("delete_layers: There remain {} layers and channels {}.".format(
            self.numlay, list(self.channels.keys())))

    def delete_channels(self, chanlist):
        """delete channels specified by 'chanlist'; if there remain no channels
        in the SPMdata structure after deletion, layers are NOT deleted, i.e.
        they preserve their metadata
        """
        
        for chan in chanlist:
            self.delchan(chan)
        print("delete_channels: Channels '{}' deleted.".format(chanlist))
        
    def sort_spm(self, key=None):
        """sort layers
        
        key - sorting function; if None, then sorting according to heights is used
        """
        
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
            
            # permute channels_counts
            chancounts = [(i in self.channels_counts[chan]) for i in range(self.numlay)]
            chancounts = chancounts[indsort]
            chancounts = [i for i in range(self.numlay) if chancounts[i]]
            
            # update all relevant layers
            for i, lay in enumerate(self.layers):
                if i not in chancounts: continue
                lay.setchan(chan, arr[i,
                    lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                    lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum
                    ])
            
            # update channel
            self.setchan(chan, arr, laylist=chancounts)
        
    def rewrite_channel(self, channame, arr):
        """replace array corresponding to 'channame' by 'arr'; all other metadata
        are preserved and layer pointer are automatically updated
        """
        
        if channame not in self.channels.keys():
            print("rewrite_channel: Invalid channel name '{}'.".format(channame))
            return
            
        # check whether in the valid data frames for each layer there are
        # NaNs; if so, halt with no change in SPMdata structure
        for i, lay in enumerate(self.layers):
            if i not in self.channels_counts[channame]: continue
            arrslice = arr[i,
                lay.xoffind[channame]:lay.xoffind[channame] + lay.xnum,
                lay.yoffind[channame]:lay.yoffind[channame] + lay.ynum]
                
            if np.any(np.isnan(arrslice.flat)):
                print("rewrite_channel: Nuexpected NaNs in a valid data frame. Halt.")
                return

        # update pointers in each layer
        for i, lay in enumerate(self.layers):
            if i not in self.channels_counts[channame]: continue
            arrslice = arr[i,
                lay.xoffind[channame]:lay.xoffind[channame] + lay.xnum,
                lay.yoffind[channame]:lay.yoffind[channame] + lay.ynum]                
            lay.setchan(channame, arrslice)

        # save new data            
        self.setchan(channame, arr, laylist=self.channels_counts[channame])
        
    def check_valid_frames(self):
        """maintenance tool for analysing whether frames defined by lay.xoffind, 
        lay.yoffind, lay.xnum, lay.ynum for each layer lay contain only valid data
        """
        
        # THIS PROCEDURE SHOULD BE AUGMENTED BY ADDITIONAL FEATURES, NOW IT ONLY
        # CHECKS WHETHER THERE ARE SOME NANS IN VALID DATA FRAMES, BUT NOT CHECK
        # WHETHER THERE ARE VALID DATA OUTSIDE THE FRAME ETC.
        
        nonans = True
        for i, lay in enumerate(self.layers):
            for chan, val in lay.channels.items():
                if np.any(np.isnan(val.flat)):
                    print("check_valid_frames: In layer no. {:3}, channel '{}'".format(
                        i, chan) + " there are redundant NaNs.")
                    nonans = False
        if nonans:
            print("check_valid_frames: There are no redundant NaNs.")
        
    def add_layers(self, laylist, heilist=None):
        """insert additional layers full of NaNs into SPMdata structure
        
        laylist - list of layer indices to be inserted, e.g. if there are 5 layers
            initially and one wants to add one layer after the first layer in SPMdata
            and two layers at the end of SPMdata, one sets laylist=[2,6,7], since the 
            final SPMdata looks like:
            [1(old), 2(new), 3(old), 4(old), 5(old), 6(new), 7(new)];
            it is assumed that laylist does not contain two identical indices, not even
            of the type: let number of layers of the final SPMdata structure be 6, then
            ind1=2 denotes the same position as ind2 = -4;
        heilist -  list of heights for new layers; if None, then lay.height is left
            as None
        """
        
        # process list of layers
        laylist = set(laylist)
        numl = len(laylist) + self.numlay   
        laylist = [ind % numl for ind in laylist]
        laylist.sort()
        numl = len(laylist) + self.numlay        
        
        if max(laylist) >= numl:
            print("add_layers: Invalid list of layer indices.")

        # process list of heights
        heilist = [None]*numl if heilist is None else heilist[:numl]
        heilist += [None]*(numl - len(heilist))

        # reallocate channel arrays
        for chan, val in self.channels.items():
            arr = np.empty((numl, val.shape[1], val.shape[2]))
            arr.fill(np.nan)
            
            complist = [i for i in range(numl) if i not in laylist]
            arr[complist] = val
            
            listcount = [(i in self.channels_counts[chan]) for i in range(self.numlay)]
            for i, ind in enumerate(laylist):
                listcount[ind:ind] = [False]
            listcount = [i for i in range(numl) if listcount[i]]

            self.setchan(chan, arr, laylist=listcount)

        # allocate new layers
        for i, ind in enumerate(laylist):
            self.layers[ind:ind] = [SPMdataLayer()]
            lay = self.layers[ind]
            lay.height = heilist[i]
        
    def add_channel(self, channame, arr, direction='forward', units=None, deep=True, 
        xofflist=None, yofflist=None, xnum=None, ynum=None, laylist=None):
        """add a new channel to the structure
        
        arr - array to represent the channel, if arr corresponds to the number of layers
            greater than self.numlay, then arr is reallocated and padded with NaNs to
            match the self.numlay, if arr corresponds to the number of layer less than
            self.numlay, then only first self.numlay layers of arr is used
        channame - name of the channel, if channame is already present in the structure,
            the original array is overriden
        direction - which scan direction the arr represents, it can be either
            'forward' or 'backward', if the latter is used and channame does not ends with
            'b', this 'b' is appended to the channame
        units - in which units the data in arr should be interpreted; if None, then 
            default UNKNOWN value is used
        deep - if True, then deep copy of arr is used
        xofflist, yofflist - x- and y-offinds for each layer; if None, then zero offinds
            are used
        xnum, ynum - lists of xnums and ynums for each layer; if None, then xnum and ynum
            of each layer is set to arr.shape[1] and arr.shape[2], respectively; these
            parameters are used only when the SPMdata structure is empty before adding
            channame; if there are already some channels in SPMdata structure, then
            these parameters ARE NOT USED
        laylist - list to be stored in channels_counts[channame]; if None, then
            all layers are put in channels_counts[channame]
        """
        
        # process the channel name
        if direction == 'backward':
            if not is_backward_scan(channame):
                print("add_channel: Invalid channel name '{}' for backward scanning.".format(
                    channame))
                return
        elif direction == 'forward':
            if is_backward_scan(channame):
                print("add_channel: Invalid channel name '{}' for forward scanning.".format(
                    channame))
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
            
            for i, lay in enumerate(self.layers):
                lay.xnum, lay.ynum = xnum[i], ynum[i]
        
        # new channel array must be at least as big as valid data frame
        xnum, ynum = max(self.xnums), max(self.ynums)
        if arr.shape[1] < xnum or arr.shape[2] < ynum:
            print(("add_channel: Input array (x, y)-shape must be at least "
                "({}, {}). Provided shape is ({}, {}).").format(xnum, ynum,
                arr.shape[1], arr.shape[2]))          
            return      
        
        # process data array
        if arr.shape[0] >= self.numlay:
            array = arr[:self.numlay]
        else:
            array = np.empty((self.numlay, arr.shape[1], arr.shape[2]))
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
                print("add_channel: Unknown units '{}'.".format(units) + \
                " Default units used instead.")
                self.units[channame] = UNITS[UNKNOWN]

        # initialize offinds and adjust them to match number of layers
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

        xn = [xofflist[i] + lay.xnum for i, lay in enumerate(self.layers)]
        yn = [yofflist[i] + lay.ynum for i, lay in enumerate(self.layers)]

        # xoffind + xnum must be always <= x-shape of the array, and similarly for y-shape 
        if (np.array(xn) > array.shape[1]).any() \
            or (np.array(yn) > array.shape[2]).any():
            print("add_channel: Index offsets specified do not match array shape.")
            return
        
        # initialize laylist and adjust it
        if laylist is None:
            print("add_channel: No laylist specified.")
            laylist = list(range(self.numlay))
        else:
            laylist = list(laylist[:self.numlay])
        
        # if channame is already present in the structure, delete it
        if channame in self.channels:
            self.delchan(channame)

        # set a channel to the structure
        self.setchan(channame, arr, laylist=laylist)

        # update data for layers        
        for i, lay in enumerate(self.layers):
            arr = self.channels[channame][i]
            lay.setchan(channame,
                arr[xofflist[i]:xofflist[i] + lay.xnum,
                    yofflist[i]:yofflist[i] + lay.ynum])
            lay.xoffind[channame] = xofflist[i]
            lay.yoffind[channame] = yofflist[i]
        
        print("add_channel: New channel '{}' successfully added.".format(channame))
        
    # PROPERTIES, AUXILIARY INQUIRING FUNCTIONS, GETTERS
        
    @property
    def arrshape(self):
        """self.arrshape = {ELEV:(xdim of 3D array corresponding to ELEV, ydim of 3D array corresponding to ELEV), ...}
        """
        return {chan: self.channels[chan].shape[1:] for chan in self.channels}
        
    def __len__(self):
    	return len(self.layers)
        
    def print_channels(self):
        """print channels present in self
        """
        
        chans = self.channels.keys()
        print("print_channels: There are {} channels: {}.".format(len(chans), list(chans)))
        
    @property
    def numlay(self):
    	"""number of layers, the same output as for len()
    	"""
    
    	return len(self.layers)
    
    @property
    def xoffinds(self):
        return {chan: [lay.xoffind[chan]
                    for lay in self.layers if chan in lay.channels]
                    for chan in self.channels.keys()}

    @property
    def yoffinds(self):
        return {chan: [lay.yoffind[chan]
                    for lay in self.layers if chan in lay.channels]
                    for chan in self.channels.keys()}
    
    @property
    def xoffsets(self):
        return [lay.xoffset for lay in self.layers]

    @property
    def yoffsets(self):
        return [lay.yoffset for lay in self.layers]
    
    def print_offsets(self):
        """print xoffset and yoffset for all layers
        """

        offs = set(zip(self.xoffsets, self.yoffsets))        
        if len(offs) == 1:
            print("print_offsets: In total {} layer{}{}".format(
                self.numlay, "s" if self.numlay != 1 else "",
        		"." if self.numlay == 0 else " with {}offsets: (xoffset, yoffset) = {} {}".format(
        		"identical " if self.numlay > 1 else "", offs.pop(), self.units[ELEV])))
        else:
            print("print_offsets: In total {} layer{}{}".format(
                self.numlay, "s" if self.numlay != 1 else "",
        		"." if self.numlay == 0 else " with offsets (xoffset, yoffset):"))
            for i, lay in enumerate(self.layers):
                print("\tLayer no. {:<3}: {}, {} {}".format(i, lay.xoffset, lay.yoffset, UNITS[ELEV]))

    @property
    def xranges(self):
        return [lay.xran for lay in self.layers]

    @property
    def yranges(self):
        return [lay.yran for lay in self.layers]
    
    def print_ranges(self):
        """print xrange and yrange for all layers
        """

        rans = set(zip(self.xranges, self.yranges))        
        if len(rans) == 1:
            print("print_ranges: In total {} layer{}{}".format(
                self.numlay, "s" if self.numlay != 1 else "",
        		"." if self.numlay == 0 else " with {}ranges: (xran, yran) = {} {}".format(
        		"identical " if self.numlay > 1 else "", rans.pop(), UNITS[ELEV])))
        else:
            print("print_ranges: In total {} layer{}{}".format(
                self.numlay, "s" if self.numlay != 1 else "",
        		"." if self.numlay == 0 else " with ranges (xran, yran):"))
            for i, lay in enumerate(self.layers):
                print("\tLayer no. {:<3}: {}, {} {}".format(i, lay.xran, lay.yran, UNITS[ELEV]))  

    @property
    def heights(self):
    	return [lay.height for lay in self.layers]
        
    def print_heights(self):
    	"""print height for all layers
    	"""
    	
    	print("print_heights: In total {} layer{}{}".format(
    	    self.numlay, "s" if self.numlay != 1 else "",
    		"." if self.numlay == 0 else " with heights:"))
    	for i, lay in enumerate(self.layers):
    	    print("\tLayer no. {:<3}: {} {}".format(i, lay.height, UNITS[ELEV]))
	
    @property
    def filenames(self):
    	return [lay.filename for lay in self.layers]
	
    def print_filenames(self):   
    	"""print filename for all layers
    	"""
    	
    	print("print_filenames: In total {} layer{}{}".format(self.numlay, "s" if self.numlay != 1 else "",
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
            print("print_nums: In total {} layer{}{}".format(
                self.numlay, "s" if self.numlay != 1 else "",
        		"." if self.numlay == 0 else " with {}dimensions: dims = {}".format(
        		"identical " if self.numlay > 1 else "", dimsl.pop())))
        else:
            print("print_nums: In total {} layer{}{}".format(
                self.numlay, "s" if self.numlay != 1 else "",
        		"." if self.numlay == 0 else " with dimensions:"))
            for i, lay in enumerate(self.layers):
                print("\tLayer no. {:<3}: dims = {}, {}".format(i, lay.xnum, lay.ynum))        

    # INDEXING

    def ind_from_height(self, height):
        """take height of the scan and return its corresponding or most approximate layer index
        """
        
        try:
            height = float(height)
        except ValueError:
            print("ind_from_height: Invalid height specification.")
            return None

        diffs = np.array(list(map(lambda x: abs(x.height - height), self.layers)))
        if diffs.min() > 1: # 1 or any other comparatively big number
            print("ind_from_height: Input height quite far from stored data.")

        return np.argmin(diffs)

    def __getitem__(self, item):
        """if the item is a string, it is assumed to be a height of some layer,
        otherwise the item is assumed to be a normal index
        """
        
        if isinstance(item, str):
            itemloc = self.ind_from_height(item)
            return self.layers[itemloc]
        else:
            return self.layers[item]

    # DATA LOADING

    def load_data(self, paths=(), suffixes=SUFFIXES, key=None):
        """loads data from external files

        paths - name of the file or tuple of filenames storing the data;
            if 'paths' is evaluated to False, then open-graphical-window
            is shown for choosing files to be loaded;
            if 'paths'='cwd', then current working directory is used to exctract
            files from;
            if 'paths' contains a string with file name, then this file is loaded;
            if 'paths' contains a string with a directory name, then all valid
            files from that directory are loaded;
            if 'paths' is a sequence of strings containing file names, then all
            valid files with these file names are loaded
        suffixes - relevant format suffixes
        key - function evaluated for each layer, according to
        	outputs of this function the loaded layers are embedded
        	into 3D arrays, if None, then sorting according to heights
        	is used
        """

	    # if 'paths' is None, than loading window is opened,
	    # otherwise 'paths' stores filenames
        if not paths:
            root = tkinter.Tk()
            root.withdraw()
            validnames = tkinter.filedialog.askopenfilenames(parent=root)
        else:
            validnames = self._process_filenames(suffixes, paths)
            
        # check whether filenames are valid
        validnames = self._check_filenames(validnames, suffixes)
        
        # if a lot of files to be loaded, warn user that it might take long time
        if len(validnames) > 15:
            print("load_data: {} files to be loaded.".format(len(validnames)))
            
        # if no file to be loaded, halt
        if len(validnames) == 0:
            print("load_data: No file to be loaded.")
            return

    	# load individual files
        self.layers = [None]*len(validnames)
        for i, name in enumerate(validnames):
            self._load_layer(i, name)            

	    # create 3D arrays in SPMdata
        self._embed_layers(key)

	    # statistics
        namelen = len(validnames)
        print("load_data: In total {} file{} loaded{}".format(namelen,
                "s" if namelen != 1 else "", "." if namelen == 0 else ":"
                ))
        for i, name in enumerate(validnames):
            print("\tFile {:2}: {}".format(i + 1, name))

    def _check_filenames(self, filenames, suffixes):
        """check whether filenames are valid
        """
        
        sufx = tuple(suffixes)

        valfilenames = []
        for name in filenames:
            if os.path.isfile(name) and name.endswith(sufx):
                valfilenames.append(name)
            else:
                print("_check_filenames: Name '{}' is not a valid file name.".format(name))
        return valfilenames 

    def _process_filenames(self, suffixes, pathnames='cwd'):
        """process input file names

        pathnames - variable containing file names, if left unspecified, then the
            current working directory is searched through
        suffixes - format suffixes which should be loaded, if None then files with all
            appropriate suffixes are loaded
        """

        try:
            if pathnames == 'cwd':
                # get files from the current working directory
                pathnamesloc = {os.path.join(os.getcwd(), item) for item in os.listdir(os.getcwd())}
                print("_process_filenames: Data read from the current working directory: {}".format(os.getcwd()))
            elif isinstance(pathnames, str):
                # if pathnames contains a name of a directory, load every file in that directory
                if os.path.isdir(pathnames):
                    pathnamesloc = {os.path.join(pathnames, item) for item in os.listdir(pathnames)}
                    print("_process_filenames: Data read from directory: {}".format(pathnames))
                # if pathnames contains a name of a file, load that file
                else:
                    pathnamesloc = set([pathnames])
            elif all(map(lambda x: isinstance(x, str), pathnames)):
                # if all items in pathnames are strings, get those
                pathnamesloc = set(pathnames)
            else:
                raise TypeError
        except TypeError:
            print("_process_filenames: Invalid data file names.")
            pathnamesloc = []

        return pathnamesloc

    def _load_layer(self, ind, filename):
        """load a single layer
        """

        self.layers[ind] = SPMdataLayer(filename)

        # update number of layers per each channel in channels_counts
        for key in self.layers[ind].channels.keys():
            self.channels_counts.setdefault(key, [])
            self.channels_counts[key].append(ind)

    def _embed_layers(self, key=None):
        """create 3D-arrays corresponding to each channel by retrieving data from all layers,
        unknown values are represented by NANs
        
        key - key for sorting the layers
        """

	    # implicit sorting of layers is according to heights
        self.layers.sort(key=((lambda x: x.height) if key is None else key))
#        self.layers.sort(key=((lambda x: x.filename) if key is None else key))
        xdims, ydims = self.xnums, self.ynums
        xnum, ynum = max(xdims), max(ydims)

        # for each channel...
        for chan in self.channels_counts.keys():
            # create empty 3D array
            chan_arr = np.empty((self.numlay, xnum, ynum))
            chan_arr.fill(np.nan)

            # fill the array with valid data and update layers
            for i, lay in enumerate(self.layers):
                chan_arr[i, :xdims[i], :ydims[i]] = lay.channels[chan].copy()
                lay.setchan(chan, chan_arr[i])
                
            # save the array to the structure
            self.setchan(chan, chan_arr)


    # DATA MANIPULATION
    	
    def setchan(self, chan, val, aver=None, laylist=None):
        """set channel to the structure
        
        chan - channel to be set
        val - value to be set to channel 'chan'
        aver - flags stating whether the channel has been averaged already
        laylist - which layers contain channel 'chan'; if None, then all 
            layers are considered
        """
        
#        if chan not in self.chan_ident:
#            print("setchan: Invalid channel '{}'.".format(chan))
#            return
        setattr(self, chan, val)
        self.channels[chan] = val
        
        if is_backward_scan(chan):
            chanl = get_forward_scan(chan)
        else:
            chanl = chan
        
        if chanl not in self.channels_fb:
            self.channels_fb[chanl] = [None, None, None]
                
        if is_backward_scan(chan):
            self.channels_fb[chanl][1] = val
            self.channels_fb[chanl][2] = False
        else:
            self.channels_fb[chan][0] = val
            
        if aver is not None:
            self.channels_fb[chanl][2] = aver

        self.channels_counts[chan] = list(range(self.numlay)) if laylist is None else laylist

    def delchan(self, chan):
        """delete channel from the structure
        
        chan - channel to be deleted
        """
        
        if chan not in self.channels:
            return
        
        del self.channels[chan]
        del self.channels_counts[chan]
        
        if is_backward_scan(chan):
            chanl = get_forward_scan(chan)
            self.channels_fb[chanl][1] = None
        else:
            chanl = chan
            self.channels_fb[chanl][0] = None
        if self.channels_fb[chanl][:2] == [None, None]:
            del self.channels_fb[chanl]
            
        for lay in self.layers:
            lay.delchan(chan)
            
        delattr(self, chan)
    	
    def is_averaged(self, chan):
        """returns True, if forward and backward scans were averaged already,
           returns None, if backward scan has not been loaded
           returns False, if forward and backward scans have not been averaged yet
        """
    
        return self.channels_fb[chan][2]
    	
    def average_channels(self, cl=None, fun=None):
        """average forward and backward scans for given channels
        
        THIS PROCEDURE SHOULD BE USED PRIOR TO ANY ALIGN PROCEDURE!!!
        
        cl - list of channels to be averaged, if None, then all channels are averaged
        fun - average function, if None, then arithmetic average is used
        """

        if fun is None:
            factor = 0.5
            avfun = lambda x, y: factor*x + (1 - factor)*y # already assumes that y is reverted
        else:
            avfun = fun
        
        clf = cl if cl is not None else self.channels.keys()
        clf = [chan for chan in clf if not is_backward_scan(chan)]
        for chan in clf:
            forscan = getattr(self, chan)
            bckscan = getattr(self, get_backward_scan(chan), None)
            if bckscan is None:
                print("average_channels: Channel '{}' has no backward version.".format(chan))
            else:                
                self.setchan(chan, avfun(forscan, bckscan))
                self.channels_fb[chan][2] = True
                
                for i, lay in enumerate(self.layers):
                    lay.setchan(chan, self.channels[chan][i])
                
                print("average_channels: Channel '{}' averaged.".format(chan))
                self.delchan(get_backward_scan(chan))

    def strip_nans(self, *cl, same=True):
        """strip NaNs from arrays corresponding to channels in 'cl', so that
        all valid data are preserved
        
        cl - sequence of channels, if empty, then all channels are stripped
        same - if True, then all channel arrays are stripped to have the
            same dimensions, this option is useful for already aligned
            arrays; if False, then each channel is stripped separately
        """
        
        if not cl:
            cl = list(self.channels.keys())
        else:
            cl = [chan for chan in cl if chan in self.channels.keys()]
            
        xmin, ymin = [None]*len(cl), [None]*len(cl)
        xmax, ymax = [None]*len(cl), [None]*len(cl)
        
        # find frame limits
        for i, chan in enumerate(cl):
            # from where to take valid data
            xmin[i] = min([lay.xoffind[chan] for lay in self.layers])    
            ymin[i] = min([lay.yoffind[chan] for lay in self.layers])    
            
            # up to where to take valid data
            xmax[i] = max([lay.xoffind[chan] + lay.xnum for lay in self.layers])
            ymax[i] = max([lay.yoffind[chan] + lay.ynum for lay in self.layers])
            
        if same:            
            # if all arrays have to have the same dimensions...
            xminv, yminv = min(xmin), min(ymin)
            xmaxv, ymaxv = max(xmax), max(ymax)

            for i, chan in enumerate(cl):
                # in the following lines there is always subtracted 
                # 'xminv' for both x and y axes!
                aux = np.empty(
                    (self.numlay, xmaxv - xminv, ymaxv - yminv), dtype=float)
                aux.fill(np.nan)                
                aux[:,
                    (xmin[i] - xminv):(xmax[i] - xminv),
                    (ymin[i] - yminv):(ymax[i] - yminv)
                    ] = self.channels[chan][:, xmin[i]:xmax[i], ymin[i]:ymax[i]]
                self.setchan(chan, aux.copy())            
        else:            
            # if channel arrays do not need to have the same dimensions...
            for i, chan in enumerate(cl):
                # take valid data
                aux = self.channels[chan][:, xmin[i]:xmax[i], ymin[i]:ymax[i]]            
                self.setchan(chan, aux.copy())
        
        # update layer pointers to new arrays as well as xoffinds and yoffinds
        for i, chan in enumerate(cl):
            for j, lay in enumerate(self.layers):
                lay.xoffind[chan] -= xminv if same else xmin[i]
                lay.yoffind[chan] -= yminv if same else ymin[i]

                lay.channels[chan] = self.channels[chan][j, 
                        lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                        lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum]

    def hew_nans(self, minl=None, maxl=None, same=True):
        """hew_nans NANs from 3D arrays corresponding to each channel, the hew_nans is performed
        for layers whose indices lie between minl and maxl, all remaining layers are
        deleted, resulting 3D arrays contain subset of purely valid data
        
        minl, maxl - cut off is performed for layers: minl, minl + 1, ..., maxl - 1, maxl
        same - if True, then all channel arrays are hewed to have the
            same dimensions, this option is useful for already aligned
            arrays; if False, then each channel is hewed separately
        """
        
        # range of layers for hew_nans
        lmin = 0 if minl is None else minl
        lmax = self.numlay if maxl is None else maxl
        lmin, lmax = lmin % (self.numlay + 1), lmax % (self.numlay + 1)
        if lmax < lmin: lmax, lmin = lmin, lmax

        print("hew_nans: Cutoff performed for layers {} through {}.".format(lmin, lmax))
            
        # delete outside layers
        self.layers = self.layers[lmin:lmax].copy()
        
        xmin, ymin = [None]*len(self.channels), [None]*len(self.channels)
        xmax, ymax = [None]*len(self.channels), [None]*len(self.channels)

        # find frame limits
        for i, chan in enumerate(self.channels):
            # from where to take valid data
            xmin[i] = max([lay.xoffind[chan] for lay in self.layers])    
            ymin[i] = max([lay.yoffind[chan] for lay in self.layers])    
            
            # up to where to take valid data
            xmax[i] = min([lay.xoffind[chan] + (lay.xnum - 1) for lay in self.layers])
            ymax[i] = min([lay.yoffind[chan] + (lay.ynum - 1) for lay in self.layers])
            
        if same:
            # if all arrays have to have the same dimensions...
            xminv, yminv = max(xmin), max(ymin)
            xmaxv, ymaxv = min(xmax), min(ymax)

            for chan in self.channels:
                # take valid data
                aux = self.channels[chan][lmin:lmax, xminv:xmaxv, yminv:ymaxv]                           
                self.setchan(chan, aux.copy())
        else:        
            # if channel arrays do not need to have the same dimensions...
            for i, chan in enumerate(self.channels):      
                # take valid data
                aux = self.channels[chan][lmin:lmax, xmin[i]:xmax[i], ymin[i]:ymax[i]]            
                self.setchan(chan, aux.copy())

        # all channels should now have the same shape
        chan = list(self.channels.keys())[0]
        
        # update layer properties
        for i, lay in enumerate(self.layers):
            # layer pointers to new arrays and offinds
            for chan in self.channels:            
                lay.channels[chan] = self.channels[chan][i]
                lay.xoffind[chan], lay.yoffind[chan] = 0, 0
                
            # xnum, ynum
            lay.xnum, lay.ynum = lay.channels[chan].shape


    # ALIGNMENT
 
    # function determining total offsets
    def _get_shifts(self, lay, chanf, chanb, xdir, ydir):
        """get total x-shifts and y-shifts
        """

        # check presence of channels
        if self.channels_fb[chanf][0] is None:
            print("align_forback: Forward channel '{}' is missing.".format(chanf))
            return None, None

        if self.channels_fb[chanf][1] is None:
            print("align_forback: Backward channel '{}' is missing.".format(chanb))
            return None, None

        # evaluate correlation and get indices of the best correlation match
        layf = self.layers[lay].channels[chanf]
        layb = self.layers[lay].channels[chanb]
        valf, valb = np.nanmean(layf), np.nanmean(layb) 
        layf[np.isnan(layf)] = valf
        layb[np.isnan(layb)] = valb
            
        relx, rely = get_offinds_corr(layf, layb)
        
        # OPET DOST BORDEL S TIM, CO JE TEDA X-OVA A Y-OVA OSA...
            
        # calculate offsets            
        if not xdir:
            relx = 0
        else:
            relx -= int(self.layers[lay].xnum / 2)
            relx += self.layers[lay].xoffind[chanf]
            relx -= self.layers[lay].xoffind[chanb]

        if not ydir:
            rely = 0
        else:
            rely -= int(self.layers[lay].ynum / 2)
            rely += self.layers[lay].yoffind[chanf]
            rely -= self.layers[lay].yoffind[chanb]

        return relx, rely

    def align_forback(self, *cl, same=True, ilay=0, xdir=True, ydir=True, order=0):
        """align all channels in 'cl' according to 'ichan' channel so that forward
        and backward scans are moved towards each other and aligned 'to
        the centre', this routine does not align layers within one channel, but
        channels among themselves instead
        
        cl - sequence of channels to be aligned, if empty, than all
            channels are aligned
        same - if True, then resulting 3D arrays for all channels in cl are of the
            same dimensions; if False, then 3D arrays for forward and backward
            direction of each physical quantity have the same dimension only
        ilay - reference layer 
        xdir - if True, than allow shift in x-direction
        ydir - if True, than allow shift in y-direction
        order - order of interpolation during shifting
        """

        # preprocess the list of channels to align, let cl contain only forward scans
        if not cl:
            cl = [chan for chan in self.channels.keys() if not is_backward_scan(chan)]
        else:
            clf = [chan for chan in cl if not is_backward_scan(chan)]
            clb = [get_forward_scan(chan) for chan in cl if is_backward_scan(chan)]
            claux, cl = list(set(clf + clb)), []
            for chan in claux:
                if chan in self.channels:
                    cl.append(chan)
                else:
                    print("aling_forback: Channel {} not present.".format(chan))
        print("align_forback: Channels to be aligned: {}.".format(cl))
        
        vecx,  vecy  = [None]*len(cl), [None]*len(cl)
        xlimf, xlimb = [None]*len(cl), [None]*len(cl)
        ylimf, ylimb = [None]*len(cl), [None]*len(cl)
        
        # align each pair of forw./backw. channels
        for i, chanf in enumerate(cl):
            chanb = get_backward_scan(chanf)
            vecx[i], vecy[i] = self._get_shifts(ilay, chanf, chanb, xdir, ydir)            

            # if xlimf != xlimb or ylimf != ylimb, then
            # resulting 3D arrays are not ensured to be
            # properly aligned w.r.t. forward/backward!
            xlimf[i], ylimf[i] = self.channels[chanf].shape[1:]
            xlimb[i], ylimb[i] = self.channels[chanb].shape[1:]
        
        # get dimensions of new arrays
        xlim = [int(np.abs(vi).max()) + max(lfi, lbi) for vi, lfi, lbi in zip(vecx, xlimf, xlimb)]
        ylim = [int(np.abs(vi).max()) + max(lfi, lbi) for vi, lfi, lbi in zip(vecy, ylimf, ylimb)]
         
        # if same == True, then modify vecx* so that all channels
        # have the same array dimension
        if same:
            xlv, ylv = max(xlim), max(ylim)
            xlim, ylim = [xlv]*len(cl), [ylv]*len(cl) 

            vecxf = [int((xlv - f - v)/2) for v, f in zip(vecx, xlimf)] 
            vecxb = [int((xlv - b + v)/2) for v, b in zip(vecx, xlimb)] 
            vecyf = [int((ylv - f - v)/2) for v, f in zip(vecy, ylimf)] 
            vecyb = [int((ylv - b + v)/2) for v, b in zip(vecy, ylimb)]
        else:
            vecxf = [0 if v > 0 else int(-v) for v in vecx]
            vecxb = [0 if v < 0 else int(+v) for v in vecx]
            vecyf = [0 if v > 0 else int(-v) for v in vecy]
            vecyb = [0 if v < 0 else int(+v) for v in vecy]
            
        # update 3D arrays
        for i, chanf in enumerate(cl):
            chanb = get_backward_scan(chanf)
            
            # create auxiliary arrays
            auxf = np.empty((self.numlay, xlim[i], ylim[i]))
            auxb = np.empty((self.numlay, xlim[i], ylim[i]))
            auxf.fill(np.nan)
            auxb.fill(np.nan)
            
            # forward channel - fill array with data
            xoff, yoff = vecxf[i], vecyf[i]
            auxf[:, xoff:xoff + xlimf[i], yoff:yoff + ylimf[i]] = self.channels[chanf]

            # backward channel - fill array with data
            xoff, yoff = vecxb[i], vecyb[i]
            auxb[:, xoff:xoff + xlimb[i], yoff:yoff + ylimb[i]] = self.channels[chanb]

            # update channels
            self.setchan(chanf, auxf)
            self.setchan(chanb, auxb)
            
            # update layer pointers and xoffinds, yoffinds
            for j, lay in enumerate(self.layers):            
                xnum, ynum = lay.xnum, lay.ynum
                
                # forward channel
                xoff, yoff = lay.xoffind[chanf] + vecxf[i], lay.yoffind[chanf] + vecyf[i]
                lay.xoffind[chanf], lay.yoffind[chanf] = xoff, yoff                
                lay.channels[chanf] = self.channels[chanf][j][xoff:xoff + xnum, yoff:yoff + ynum]
#                print("for: xoff, yoff: ", xoff, yoff)

                # backward channel
                xoff, yoff = lay.xoffind[chanb] + vecxb[i], lay.yoffind[chanb] + vecyb[i]
                lay.xoffind[chanb], lay.yoffind[chanb] = xoff, yoff
                lay.channels[chanb] = self.channels[chanb][j][xoff:xoff + xnum, yoff:yoff + ynum]
#                print("bck: xoff, yoff: ", xoff, yoff)

    def _showfit_align(self, coefs, fac, xn, yn, xoffmin, yoffmin):
        """calculate fit points to be shown in procedure align_corr_chan_fac
        """
    
        acoefx, bcoefx, acoefy, bcoefy = coefs
    
        xfit, yfit = None, None
    
        if acoefx is not None:
            xfit = acoefx*np.arange(self.numlay) + bcoefx        
            xfit = np.ceil(xfit - fac * xn / 2)
            xfit -= xoffmin
            xfit[xfit < 0] = 0

        if acoefy is not None:
            yfit = acoefy*np.arange(self.numlay) + bcoefy
            yfit = np.ceil(yfit - fac * yn / 2)
            yfit -= yoffmin
            yfit[yfit < 0] = 0

        return xfit, yfit
        
    def _show_align(self, rellay, revactlay, mislay,
        xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs):
        """show alignments in procedures align_corr_chan and align_corr_chan_fac
        """

        # create plot
        fig, ax = plt.subplots(1, 1, num="Extrapolation trends and control points")
        ax.set_title("Extrapolation trends and control points")
        ax.set_xlim([0, self.numlay])

        # show offsets for all layers
        ax.step(range(self.numlay), xoffind, "b-", label="x offsets",
            where="mid", linewidth=2.5)
        ax.step(range(self.numlay), yoffind, "g-", label="y offsets",
            where="mid", linewidth=2.5)
        
        # show areas with inactive layers
        pending = False
        for i in revactlay:
            if not pending: i0 = i
            if i + 1 in revactlay:
                pending = True
                continue
            else:
                plt.axvspan(i0 - 0.5, i + 0.5, color='r', alpha=0.5)
                pending = False

        # if there are layers whose offsets are to be inter-/extra-polated
        # plot control points, calculated points and fits
        if mislay:
            # manual setting?
            manstrx = "" if xcoefs is None else " (manual)"
            manstry = "" if ycoefs is None else " (manual)"
                    
            # plot reference points used for fitting
            if xcoefs is None:
                ax.plot(rellay, xfit[rellay], "bo", label="x control points", markersize=10)
            if ycoefs is None:        
                ax.plot(rellay, yfit[rellay], "go", label="y control points", markersize=10)

            # plot inter-/extra-polated points
            ax.plot(mislay, xoffind[mislay], "bv", label="x calculated points", markersize=10)
            ax.plot(mislay, yoffind[mislay], "gv", label="y calculated points", markersize=10)

            # plot fits
            ax.plot(range(self.numlay), xfit, "b--",
                label="x fit" + manstrx + "\nslope = {:.4}".format(coefs[0]))        
            ax.plot(range(self.numlay), yfit, "g--",
                label="y fit" + manstry + "\nslope = {:.4}".format(coefs[2]))        
                    
        ax.legend(loc='best')
        plt.show()
        
    def _preprocess_align(self, ichan, actlayi, mislayi, rellayi, xcoefs, ycoefs, nlay=None, rnlay=None):
        """preprocessing stage of procedure align_corr_chan_fac
        """
    
        # reference channel
        if ichan is None:
            print("align_corr_chan_fac: No reference channel specified.")
            return

        # active layers
        if actlayi is None:
            actlay = range(self.numlay)
        else:
            actlay = {lay % self.numlay for lay in actlayi}
        actlay = sorted(list(actlay))                  
                    
        if not actlay:
            print("align_corr_chan_fac: No layer to be aligned.")
            return
        
        # nonactive layers
        revactlay = [lay for lay in range(self.numlay) if lay not in actlay]
        
        # reference layer for (non)active layers
        nlay = actlay[-1] if nlay is None or nlay not in actlay else nlay
        if revactlay:
            rnlay = revactlay[0] if rnlay is None or rnlay in actlay else rnlay
        else:
            rnlay = None
        print("align_corr_chan_fac: Active layers: {}".format(actlay))
        print("align_corr_chan_fac: nlay, rnlay: {}, {}".format(nlay, rnlay))

        # if there are layers to fit, preprocess mislay and rellay
        if mislayi == 'all':
            mislay = list(range(self.numlay))
        else:
            mislay = mislayi
        
        rellay = ()
        if mislay:
            mislay = list({lay % self.numlay for lay in mislay})
            
            # if at least one manual alignment coefficient is not given
            # use fitting by rellay
            if xcoefs is None or ycoefs is None:
                if not rellayi:
                    rellay = [lay for lay in range(self.numlay) if lay not in mislay]            
                else:
                    rellay = list({lay % self.numlay for lay in rellayi})
                    mislay = [lay for lay in mislay if lay not in rellay]
                rellay = sorted(rellay)
    
            mislay = sorted(mislay)
            
            print("align_corr_chan_fac: Refer. layers: {}".format(rellay))
            print("align_corr_chan_fac: Missi. layers: {}".format(mislay))

        return actlay, revactlay, mislay, rellay, nlay, rnlay

    def _corr_align(self, ichan, actlay, rellay, mislay, ilay, arr):
        """cross-correlation stage of procedure align_corr_chan_fac
        """
        
        print("align_corr_chan_fac: Alignment according to channel '{}'.".format(ichan))
        relx, rely = np.zeros(self.numlay), np.zeros(self.numlay)
        arr1, k = arr[ilay], 1        
        arr1[np.isnan(arr1)] = np.nanmean(arr1)
        
        # for each layer in actlay or rellay find out the most matching offsets
        arlist = sorted([lay for lay in range(self.numlay) if lay in actlay or lay in rellay])
        
        print("align_corr_chan_fac: Cross-correlation:")
        for i in arlist:
            print("\tLayer {:2} out of {} being processed...".format(k, len(arlist)))
            k += 1
            if mislay and i in mislay: continue

            # evaluate correlation and get indices of the best correlation match
            arrtemp, arr1tem = arr[i], arr1
            arrtemp[np.isnan(arrtemp)] = np.nanmean(arrtemp)
            arr1tem[np.isnan(arrtemp)] = np.nanmean(arr1tem)

            relx[i], rely[i] = get_offinds_corr(arr1tem, arrtemp)

            # OPET DOST BORDEL S TIM, CO JE TEDA X-OVA A Y-OVA OSA...
        
        return relx, rely
  
    def _fit_align(self, xcoefs, ycoefs, mislay, rellay, relx, rely):
        """fitting stage of procedure align_corr_chan_fac
        """

        acoefx, bcoefx = None, None
        acoefy, bcoefy = None, None
        xfit, yfit = None, None

        # inter-/extra-polation
        if mislay:
            # find linear fit
            from scipy.stats import linregress
                
            if xcoefs is None:
                acoefx, bcoefx = linregress(rellay, relx[rellay])[:2]                
            else:
                acoefx, bcoefx = xcoefs

            if ycoefs is None:
                acoefy, bcoefy = linregress(rellay, rely[rellay])[:2]                
            else:
                acoefy, bcoefy = ycoefs
       
            # calculate offsets for layers in mislay
            for i in mislay:
                relx[i] = acoefx*i + bcoefx
                rely[i] = acoefy*i + bcoefy

        return relx, rely, acoefx, bcoefx, acoefy, bcoefy
        
    def _postprocess_align(self, ichan, relx, rely, fac, revactlay, actlay, nlay, rnlay):
        """postprocessing stage of procedure align_corr_chan_fac
        """
        
        xn, yn = np.array(self.xnums), np.array(self.ynums)
        xo = np.array([lay.xoffind[ichan] for lay in self.layers])
        yo = np.array([lay.yoffind[ichan] for lay in self.layers])

        # offsets        
        relx = np.round(relx - fac * xn / 2)
        rely = np.round(rely - fac * yn / 2)

        # reference offsets with respect to the array origin
        xoffind = relx
        yoffind = rely
        
        # for nonactive layers use the original offsets
        if revactlay:
            print("align_corr_chan_fac: Nonactive layer offsets may be affected.")
            # now we set xoffind and yoffind of all channels to values of
            # xoffind and yoffind of the ichan channel
            xoffind[revactlay] = fac * xo[revactlay] + xoffind[nlay] - fac * xo[rnlay]
            yoffind[revactlay] = fac * yo[revactlay] + yoffind[nlay] - fac * yo[rnlay]
            
        # all offsets with respect to the array origin
        xoffmin, yoffmin = xoffind.min(), yoffind.min()
        xoffind -= xoffmin
        yoffind -= yoffmin

        # xoffmin and yoffmin are used in demonstrating plot
        xoffmin += relx.min()
        yoffmin += rely.min()

        return xn, yn, xoffind, yoffind, xoffmin, yoffmin

    def _upsample_align(self, ichan, fac, chanlist, order):
        """upsample stage of procedure align_corr_chan_fac
        """
        
        # upsampling of reference channel
        arr, arrdict = [None]*self.numlay, {}

        print("align_corr_chan_fac: Channel '{}' being upsampled.".format(ichan))
        for i, lay in enumerate(self.layers):        
            xno, yno = lay.channels[ichan].shape                
            x = np.linspace(0, xno - 1, fac * xno)
            y = np.linspace(0, yno - 1, fac * yno)
            coords = np.meshgrid(x, y, indexing='ij')
            arr[i] = map_coordinates(lay.channels[ichan], coords, order=order, cval=np.nan) # cval!!!

        # upsampling of channels in chanlist
        # assume that xnum and ynum is equal for all channels to xnum
        # and ynum of the reference channel for each layer
        for chan in chanlist:
            print("align_corr_chan_fac: Channel '{}' being upsampled.".format(chan))
            arrtmp = [None]*self.numlay
                        
            for i, lay in enumerate(self.layers):
                xno, yno = lay.channels[chan].shape                
                x = np.linspace(0, xno - 1, fac * xno)
                y = np.linspace(0, yno - 1, fac * yno)
#                x = np.arange(fac * xno)
#                y = np.arange(fac * yno)

                coords = np.meshgrid(x, y, indexing='ij')
                arrtmp[i] = map_coordinates(lay.channels[chan], coords, order=order, cval=np.nan)                
                
            arrdict[chan] = arrtmp

        return arr, arrdict
       
    def _reallocate_align(self, chanlist, xoffind, yoffind, fac, xn, yn, arrdict, downsample, order=0):
        """reallocating stage of procedure align_corr_chan_fac
        """
    
        # rescale valid data dimensions
        xn, yn = fac * xn, fac * yn
        auxlist = {}
    
        # dimensions of the new array
        xlim = max(xoffind + xn) - min(xoffind)
        ylim = max(yoffind + yn) - min(yoffind)
    
        # align each channel...
        for chan in chanlist:
            
            # allocate new array for the channel 'chan'
            aux = np.empty((self.numlay, xlim, ylim), dtype=float)    
            aux.fill(np.nan)

            # fill new array with valid data and appropriate offsets
            for i in range(self.numlay):      
                aux[i, 
                    xoffind[i]:xoffind[i] + xn[i],
                    yoffind[i]:yoffind[i] + yn[i]
                    ] = arrdict[chan][i].copy()

            # store channel
            auxlist[chan] = aux.copy()

        # downsampling and updating data attributes...
        if not downsample: fac = 1

        for chan in chanlist:
            aux = auxlist[chan]
                    
            # save nums and offinds for each layer
            for i, lay in enumerate(self.layers):                
                lay.xnum = int(xn[i] / fac)
                lay.ynum = int(yn[i] / fac)
                lay.xoffind[chan] = math.floor(xoffind[i] / fac)
                lay.yoffind[chan] = math.floor(yoffind[i] / fac)

            # downsampling                    
            if downsample:            
                print("align_corr_chan_fac: Channel '{}' being downsampled.".format(chan))
                x = np.linspace(0, aux.shape[1], math.ceil(aux.shape[1] / fac), endpoint=False)
                y = np.linspace(0, aux.shape[2], math.ceil(aux.shape[2] / fac), endpoint=False)
                z = range(self.numlay)
                coords = np.meshgrid(z, x, y, indexing='ij')

                # ZDE JE DALSI PROBLEM V TOM, ZE NAN S CIMKOLI DA OPET NAN, TJ. PRI INTERPOLACI
                # SE DO AUX NA MISTA, KDE BY MOHLA BYT PLATNA DATA, DOSTANE NAN A ROZMERY VYSLEDNE
                # OBLASTI PLATNYCH DAT JSOU PO DOWNSAMPLINGU OBECNE MENSI NEZ PUVODNI ROZMERY
                # ULOZENE V self.layers[i].xnum A self.layers[i].ynum
                
                # during interpolation in map_coordinates anything with NaN is again NaN, so
                # it may be the case that valid data are truncated due to padding NaNs, to
                # circumvent this we replace NaNs with an average over aux
                aux[np.isnan(aux)] = np.nanmean(aux)
                
                # aux now stores downsampled 3D array for channel 'chan', but instead of NaNs,
                # np.nanmean(aux) is used as a padding
                aux = map_coordinates(aux, coords, order=order, mode='nearest')#cval=np.nan)
            
                # store tentative aux in auxlist
                auxlist[chan] = aux.copy()
                    
                # create new aux filled with NaNs
                aux = np.empty_like(auxlist[chan])
                aux.fill(np.nan)

                # to aux insert valid data from auxlist
                for i, lay in enumerate(self.layers):
                    aux[i,
                        lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                        lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum
                        ] = auxlist[chan][i,
                        lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                        lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum                    
                        ].copy()
                    
            # update channel
            self.setchan(chan, aux.copy())                    

            # for each layer update pointer to valid data
            for i, lay in enumerate(self.layers):
                lay.channels[chan] = self.channels[chan][i,
                    lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                    lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum]                

#        print("PO:")
#        self.show_layer(FREQ, ilay=4)            

    def align_corr_chan_fac(self, *cl, ichan=None, ilay=0, nlay=None, rnlay=None, actlay=None, rellay=None, mislay=None, show=False,
        xcoefs=None, ycoefs=None, fac=1, downsample=True, order=1):
        """align arrays according to ichan channel, i.e. offsets are calculated only once for ichan and then used
        multiple times for all channels in cl, unknown values are represented by NANs,
        resulting array may have bigger dimensions
        than original arrays, this routine cannot work with NANs (!) (due to correlate2d)
        
        cl - list of channels to be aligned, if empty, then all channels are aligned
        ilay - all layers are compared to fixed layer with index 'lay'
        nlay - if actlay is not None or not all layers, then nlay specifies to which layer inactive
            layers should adhere, i.e. offind of rnlay layer will be set to the offind of nlay;
            if nlay is None or lies outside actlay, then it is set to be the last layer in actlay
        rnlay - analog of nlay for inactive layers, if None or in actlay, then it is set to be the
            first layer not in actlay
        ichan - all channels in cl are aligned according to channel ichan
        actlay - list of 'active' layers, which are to be aligned, if None, then all layers are used
        rellay - list of 'relevant' layers which are used to determine offsets for layers in mislay,
            if None or empty, then rellay is effectively put equal to all layers not lying in mislay
        mislay - list of 'missing' layers for which offsets should be calculated by inter-/extra-polation,
            if None or empty, then correlation calculation is performed for all layers and no 
            inter-/extra-polation is used, if mislay and rellay are not disjoint, then mislay is
            put equal to mislay minus rellay, if mislay='all', then all layers are used for
            inter-/extra-polation
        show - if True, then plot depicting rellay points and extrapolation trends is shown,
            control points - rellay points
            calculated points - mislay points
            x fit, y fit - fit of x-offsets and y-offsets, respectively
            (manual) - shown if fit is set manually via xcoefs or ycoefs
            x offsets, y offsets - resulting x-offsets and y-offsets, respectively
            red area - for this area no alignment is active, i.e. corresponding layers do not
                lay in actlay
        xcoefs, ycoefs - linear and absolute coefficients for manual interpolation for x and y
            dimension, respectively; if not None, then these parameters are used instead of
            automatically calculated ones, in which case xcoefs=(xlincoef, xabscoef) etc.;
            if not None, then rellay is not used
        fac - upsampling factor
        downsample - whether to perform downsampling after correlation
        order - order of interpolation for resampling
        
        BEWARE - there might be problems for interpolation order > 1; due to round off errors
            the self.layers[i][chan] arrays may contain NaNs with which map_coordinates cannot
            work for order > 1; this issue takes place in '_reallocate_align' routine in
            downsampling stage, specifically in command np.linspace in its third parameter
     
            
        Example 1: I have given SPMdata structure 'dta' with 30 layers. I want to align
            channels CURR and PHAS. As a reference channel I want to choose FREQ,
            as a reference layer I choose layer no. 3. That is, I write:
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, ilay=3)'.
            But, from the result I see that cross-correlation does not work well for layers
            no. 10, 11, 12, 15, 16. I would like to use linear fit for offsets of these 
            layers. So I append 'mislay=[10, 11, 12, 15, 16]' to argument list.
            New result I obtain nevertheless does not fit offsets for these layers well.
            Reason for this may be that offsets for layers outside 'mislay' vary a lot, so I
            would like to calculate a linear fit only from offsets pertaining to layers
            no. 5, 6, 7. Thus I append 'rellay=[5, 6, 7]' to argument list and also set 
            'show' to True to see how offsets are calculated. At the end I have
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, ilay=3, mislay=[10, 11, 12, 15, 16],
                rellay=[5, 6, 7])'.
                
        Example 2: Suppose I have given SPMdata structure as in Example 1. But now 
            I want to align according to layer no. 0, which is default choice. Now I would
            like to align only layers no. 10 through no. 29 and leave layers no. 0 through
            no. 5 unaffected (i.e. layers no. 0 through no. 5 of channel 'ichan'! unaffected,
            since routine 'align_corr_chan' aligns channels according to channel 'ichan', all
            channels in 'cl' different from 'ichan' may get different offsets even
            for layers no. 0 through no. 5). So I write:
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, actlay=range(6, 30), show=True,
                mislay=range(10, 30))'.
            But, the resulting alignment is not good, so I read actual offsets used by
            routine from the plot shown and set fit coefficients for x-offsets by hand,
            so I write
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, actlay=range(6, 30), show=True,
                mislay=range(10, 30), xcoefs=(5., 10.))'.        
        """

        chanlist = list(self.channels if not cl else cl)
        chanlist = [chan for chan in chanlist if chan in self.channels]
        
        # since each layer contains attributes xnum and ynum, i.e. dimensions of valid
        # data area, which are identical for all channels, let us exclude possibility,
        # that 'align_corr_chan_fac' changes dimension of valid data area for different
        # channels
        if not downsample and fac > 1 and len(chanlist) < len(self.channels):
            print("align_corr_chan_fac: Downsampling on, but not all channels included.")
            return

        # layer lists initialization
        actlay, revactlay, mislay, rellay, nlay, rnlay = self._preprocess_align(
                                    ichan, actlay, mislay, rellay, xcoefs, ycoefs, nlay, rnlay)

        # upsampling        
        arr, arrdict = self._upsample_align(ichan, fac, chanlist, order)
       
        # cross-correlation         
        relx, rely = self._corr_align(ichan, actlay, rellay, mislay, ilay, arr)
    
        # inter-/extra-polation
        relx, rely, *coefs = self._fit_align(xcoefs, ycoefs, mislay, rellay, relx, rely)

        # postprocessing
        xn, yn, xoffind, yoffind, xoffmin, yoffmin = self._postprocess_align(
                                    ichan, relx, rely, fac, revactlay, actlay, nlay, rnlay)

        # plotting - show trends and control points
        if show:
            xfit, yfit = self._showfit_align(coefs, fac, xn, yn, xoffmin, yoffmin)
            self._show_align(rellay, revactlay, mislay,
                        xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs)

        # reallocation of arrays according to offsets and downsampling
        self._reallocate_align(chanlist, xoffind, yoffind, fac,
            xn, yn, arrdict, downsample, order)

        print("align_corr_chan_fac: Channels {} aligned.".format(chanlist))            

    def _showfit_align_trend(self, coefs, fac, xn, yn, xoffmin, yoffmin):
        """calculate fit points to be shown in procedure align_trend
        """
    
        acoefx, bcoefx, acoefy, bcoefy = coefs
        xfit, yfit = None, None
        
        if acoefx is not None:
            xfit = acoefx*np.arange(self.numlay) + bcoefx
            xfit = xfit - fac * xn / 2
            xfit -= xoffmin
            xfit = np.round(xfit / fac)
            xfit[xfit < 0] = 0
        
        if acoefy is not None:
            yfit = acoefy*np.arange(self.numlay) + bcoefy
            yfit = yfit - fac * yn / 2
            yfit -= yoffmin
            yfit = np.round(yfit / fac)
            yfit[yfit < 0] = 0

        return xfit, yfit

    def _reallocate_align_trend(self, chanlist, xoffind, yoffind, fac, xn, yn, arrdict, order=0):
        """reallocating stage of procedure align_trend
        """
    
        # exact offinds
        xoflt, yoflt = xoffind, yoffind
        
        # integral offinds
        xoffind, yoffind = np.round(xoflt).astype(int), np.round(yoflt).astype(int)
#        xoffind, yoffind = np.floor(xoflt).astype(int), np.floor(yoflt).astype(int)
    
        # dimensions of the new array
        xlim = max(xoffind + xn) - min(xoffind)
        ylim = max(yoffind + yn) - min(yoffind)
    
        print("xoffind - xoflt: ", xoffind - xoflt)
    
        # align each channel...
        for chan in chanlist:
            # allocate new array for the channel 'chan'
            aux = np.empty((self.numlay, xlim, ylim), dtype=float)
            aux.fill(np.nan)

            # fill new array with valid data and appropriate offsets            
            for i, lay in enumerate(self.layers):
                x = (xoffind[i] - xoflt[i]) + np.arange(lay.xnum)
                y = (yoffind[i] - yoflt[i]) + np.arange(lay.ynum)
                coords = np.meshgrid(x, y, indexing='ij')
                aux[i, 
                    xoffind[i]:xoffind[i] + xn[i],
                    yoffind[i]:yoffind[i] + yn[i]
                    ] = map_coordinates(
                        lay.channels[chan], coords,
                        order=order, mode='nearest')

            # store channel
            self.setchan(chan, aux.copy())
            
            # update offinds and channels for each layer
            for i, lay in enumerate(self.layers):
                lay.xoffind[chan] = xoffind[i]
                lay.yoffind[chan] = yoffind[i]

                lay.channels[chan] = self.channels[chan][i, 
                    xoffind[i]:xoffind[i] + xn[i],
                    yoffind[i]:yoffind[i] + yn[i]]

    def _upsample_align_trend(self, ichan, fac, order):
        """upsample stage of procedure align_trend
        """
        
        # upsampling of reference channel
        arr = [None]*self.numlay

        print("align_corr_chan_fac: Channel '{}' being upsampled.".format(ichan))
        for i, lay in enumerate(self.layers):        
            xno, yno = lay.channels[ichan].shape                
            x = np.linspace(0, xno - 1, fac * xno)
            y = np.linspace(0, yno - 1, fac * yno)
            coords = np.meshgrid(x, y, indexing='ij')
            arr[i] = map_coordinates(lay.channels[ichan], coords, order=order)#, cval=np.nan) # cval!!!

        return arr

    def _postprocess_align_trend(self, ichan, relx, rely, fac, revactlay, actlay):
        """postprocessing stage of procedure align_trend
        """
        
        xn, yn = np.array(self.xnums), np.array(self.ynums)
        xo = np.array([lay.xoffind[ichan] for lay in self.layers])
        yo = np.array([lay.yoffind[ichan] for lay in self.layers])

        # offsets        
        relx = relx - fac * xn / 2
        rely = rely - fac * yn / 2

        relx[revactlay] = relx[actlay].min()
        rely[revactlay] = rely[actlay].min()

        # reference offsets with respect to the array origin
        xoffind = relx - relx.min()
        yoffind = rely - rely.min()
        
        # for nonactive layers use the original offsets
        if revactlay:
            print("align_corr_chan_fac: Nonactive layer offsets may be affected.")
            # now we set xoffind and yoffind of all channels to values of
            # xoffind and yoffind of the ichan channel
            xoffind[revactlay] = fac * xo[revactlay]
            yoffind[revactlay] = fac * yo[revactlay]
            
        # all offsets with respect to the array origin
        xoffmin, yoffmin = xoffind.min(), yoffind.min()
        xoffind = (xoffind - xoffmin).astype(float) / fac
        yoffind = (yoffind - yoffmin).astype(float) / fac

        xoffmin += relx.min()
        yoffmin += rely.min()

        return xn, yn, xoffind, yoffind, xoffmin, yoffmin
                    
    def align_trend(self, *cl, ichan=None, ilay=0, actlay=None, rellay=None, mislay=None, show=False,
        xcoefs=None, ycoefs=None, fac=1, order=1):
        """align arrays according to ichan channel, i.e. offsets are calculated only once for ichan and then used
        multiple times for all channels in cl, unknown values are represented by NANs,
        resulting array may have bigger dimensions than original arrays
        
        cl - list of channels to be aligned, if empty, then all channels are aligned
        ilay - all layers are compared to fixed layer with index 'lay'
        ichan - all channels in cl are aligned according to channel ichan
        actlay - list of 'active' layers, which are to be aligned, if None, then all layers are used
        rellay - list of 'relevant' layers which are used to determine offsets for layers in mislay,
            if None or empty, then rellay is effectively put equal to all layers not lying in mislay
        mislay - list of 'missing' layers for which offsets should be calculated by inter-/extra-polation,
            if None or empty, then correlation calculation is performed for all layers and no 
            inter-/extra-polation is used, if mislay and rellay are not disjoint, then mislay is
            put equal to mislay minus rellay, if mislay='all', then all layers are used for
            inter-/extra-polation
        show - if True, then plot depicting rellay points and extrapolation trends is shown,
            control points - rellay points
            calculated points - mislay points
            x fit, y fit - fit of x-offsets and y-offsets, respectively
            (manual) - shown if fit is set manually via xcoefs or ycoefs
            x offsets, y offsets - resulting x-offsets and y-offsets, respectively
            red area - for this area no alignment is active, i.e. corresponding layers do not
                lay in actlay
        xcoefs, ycoefs - linear and absolute coefficients for manual interpolation for x and y
            dimension, respectively; if not None, then these parameters are used instead of
            automatically calculated ones, in which case xcoefs=(xlincoef, xabscoef) etc.;
            if not None, then rellay is not used
        fac - upsampling factor
        order - order of interpolation for resampling
        
        BEWARE - there might be problems for interpolation order > 1; due to round off errors
            the self.layers[i][chan] arrays may contain NaNs with which map_coordinates cannot
            work for order > 1; this issue takes place in '_reallocate_align' routine in
            downsampling stage, specifically in command np.linspace in its third parameter
     
            
        Example 1: I have given SPMdata structure 'dta' with 30 layers. I want to align
            channels CURR and PHAS. As a reference channel I want to choose FREQ,
            as a reference layer I choose layer no. 3. That is, I write:
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, ilay=3)'.
            But, from the result I see that cross-correlation does not work well for layers
            no. 10, 11, 12, 15, 16. I would like to use linear fit for offsets of these 
            layers. So I append 'mislay=[10, 11, 12, 15, 16]' to argument list.
            New result I obtain nevertheless does not fit offsets for these layers well.
            Reason for this may be that offsets for layers outside 'mislay' vary a lot, so I
            would like to calculate a linear fit only from offsets pertaining to layers
            no. 5, 6, 7. Thus I append 'rellay=[5, 6, 7]' to argument list and also set 
            'show' to True to see how offsets are calculated. At the end I have
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, ilay=3, mislay=[10, 11, 12, 15, 16],
                rellay=[5, 6, 7])'.
                
        Example 2: Suppose I have given SPMdata structure as in Example 1. But now 
            I want to align according to layer no. 0, which is default choice. Now I would
            like to align only layers no. 10 through no. 29 and leave layers no. 0 through
            no. 5 unaffected (i.e. layers no. 0 through no. 5 of channel 'ichan'! unaffected,
            since routine 'align_corr_chan' aligns channels according to channel 'ichan', all
            channels in 'cl' different from 'ichan' may get different offsets even
            for layers no. 0 through no. 5). So I write:
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, actlay=range(6, 30), show=True,
                mislay=range(10, 30))'.
            But, the resulting alignment is not good, so I read actual offsets used by
            routine from the plot shown and set fit coefficients for x-offsets by hand,
            so I write
            'dta.align_corr_chan(CURR, PHAS, ichan=FREQ, actlay=range(6, 30), show=True,
                mislay=range(10, 30), xcoefs=(5., 10.))'.        
        """

        chanlist = list(self.channels if not cl else cl)
        chanlist = [chan for chan in chanlist if chan in self.channels]
        
        # layer lists initialization
        actlay, revactlay, mislay, rellay = self._preprocess_align(
                                    ichan, actlay, mislay, rellay, xcoefs, ycoefs)

        # upsampling        
        arr = self._upsample_align_trend(ichan, fac, order)
       
        # cross-correlation         
        relx, rely = self._corr_align(ichan, actlay, rellay, mislay, ilay, arr)
    
        # inter-/extra-polation
        relx, rely, *coefs = self._fit_align(xcoefs, ycoefs, mislay, rellay, relx, rely)

        # postprocessing
        xn, yn, xoffind, yoffind, xoffmin, yoffmin = self._postprocess_align_trend(
                                    ichan, relx, rely, fac, revactlay, actlay)

        # plotting - show trends and control points
        if show:
            xfit, yfit = self._showfit_align_trend(coefs, fac, xn, yn, xoffmin, yoffmin)
            self._show_align(rellay, revactlay, mislay,
                        xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs)

        # reallocation of arrays according to offsets
        self._reallocate_align_trend(chanlist, xoffind, yoffind, fac,
            xn, yn, order)
    
        print("align_trend: Channels {} aligned.".format(chanlist))  
        
    def align_offsets(self, *cl):
        """align arrays corresponding to each channel according their offsets,
        unknown values are represented by NANs, resulting array may have bigger dimensions than original arrays,
        it is assumed that xranges and yranges of all layers are identical
        
        cl - list of channels to be aligned, if None, then all channels are aligned
        """

        # POZOR!!! obracene x a y ve dvou prikazech nize, obr. takto konecne vypada dobre, ale je zde nekonzistence
        # v tom, jake parametry nahravam do x a y slozek ze vstupniho souboru nanonisem, patrne
        # PROBLEM JE V TOM, ZE PRVNI INDEX POLE SE VYKRESLUJE DO RADKU OBRAZKU, ZATIMCO DRUHY INDEX POLE SE VYKRESLUJE
        # DO SLOUPCU OBRAZKU, TJ. INDEX x V NPARRAY JE CHAPAN JAKO OSA y V OBRAZKU
        # KTEROU Z TECHTO MOZNOSTI POUZIVA NANONIS ??????

        yoffsets, xoffsets = self.xoffsets, self.yoffsets
#        xoffsets, yoffsets = self.xoffsets, self.yoffsets

        xn, yn = self.xnums, self.ynums
        xdimnew, ydimnew = max(xn), max(yn)
        
        # it is assumed that xranges and yranges of all layers are identical, xstep and ystep
        # are thus calculated according to ranges of the first layer
        xstep, ystep = self.xranges[0] / xdimnew, self.yranges[0] / ydimnew
        xdimnew += math.ceil((max(xoffsets) - min(xoffsets)) / xstep)
        ydimnew += math.ceil((max(yoffsets) - min(yoffsets)) / ystep)

        # xoffind = [math.floor((offset - min(xoffsets)) / xstep) for offset in xoffsets]
        # yoffind = [math.floor((offset - min(yoffsets)) / ystep) for offset in yoffsets]
        xoffind = [math.ceil((offset - min(xoffsets)) / xstep) for offset in xoffsets]
        yoffind = [math.ceil((offset - min(yoffsets)) / ystep) for offset in yoffsets]

        for chan in (self.channels if not cl else cl):
            # create new array
            aux = np.empty((self.numlay, xdimnew, ydimnew), dtype=float)
            aux.fill(np.nan)

            # update lay.xoffind and lay.yoffind
            for i, lay in enumerate(self.layers):
                lay.xoffind[chan] = int(xoffind[i])
                lay.yoffind[chan] = int(yoffind[i])

            for i, lay in enumerate(self.layers):
                aux[i, 
                    xoffind[i]:xoffind[i] + xn[i],
                    yoffind[i]:yoffind[i] + yn[i]
                    ] = lay.channels[chan].copy()
                lay.channels[chan] = aux[i, xoffind[i]:xoffind[i] + xn[i],
                    yoffind[i]:yoffind[i] + yn[i]]
            
            self.setchan(chan, aux)
            
    # DISPLAY PROCEDURES
    
    def show_layer(self, *cl, ilay=0, rangelayer='each'):
        """show valid data for channel 'chan' in layer 'ilay'
        
        cl - sequence of channels to be shown
        ilay - which layer to show
        rangelayer - which layer to choose for scaling; if 'each',
            then scaling is based on layer 'ilay'; if 'all', then
            scaling is based on maximum and minimum of the whole
            3D array pertaining to channels in 'cl'
        """
        
        cl = [chan for chan in cl if chan in self.layers[ilay].channels]
        if len(cl) == 0:
            print("show_layer: No valid channel to show.")
            return
        
        fig, axes = plt.subplots(1, len(cl), squeeze=False,
            num="Valid data for layer no. {}".format(ilay))
        
#        print(type(axes), axes.shape, axes)
        
        for ax, chan in zip(axes.flat, cl):
            if rangelayer == 'each':
                vmin, vmax = None, None
            elif rangelayer == 'all':
                vmin = np.nanmin(self.channels[chan])
                vmax = np.nanmax(self.channels[chan])
            else:
                pass
#                vmin = np.nanmin(self.channels[chan][rangelayer])
#                vmax = np.nanmax(self.channels[chan][rangelayer])
                
            ax.imshow(self.layers[ilay].channels[chan], vmin=vmin, vmax=vmax)#, origin='lower')
            ax.set_title("Channel '{}' [units: {}]".format(chan, self.units[chan]))
        
#        print("dims: ", self.layers[ilay].channels[chan].shape)
#        print("dims: ", self.layers[ilay].xnum, self.layers[ilay].ynum)
#        print("offi: ", self.layers[ilay].xoffind[chan], self.layers[ilay].yoffind[chan])
#        xoff = find_offind(0, self.channels[chan][ilay].T, 100)
#        yoff = find_offind(0, self.channels[chan][ilay]  , 100)
#        print("offi: ", xoff, yoff)
#        count_nans(self.channels[chan][ilay][xoff + 2])
#        count_nans(self.channels[chan][ilay][:,yoff + 2])
        
        plt.show()
            
    def show_channel(self, chan, laylist=None, limit=8, rowlimit=4, aspect='equal',
        uniscale='present', interpolation='none', origin='lower'):
        """plot channel for various layers
        
        chan - which channel to plot
        laylist - which layers to plot, sequence containing indices, default is first 'limit' layers
        limit - maximum number of layers plotted
        rowlimit - maximum number of plots per row
        aspect - aspect of imshow
        uniscale - how to scale plots, for 'all' plots are rescaled to global maxima and minima of
        	the channel throughout all layers, for None each plot is rescaled individually,
        	otherwise it is assumed that uniscale is a list of layers over which the maxima and
        	minima are computed and used for plot scaling, for 'present' the same list as 'laylist'
        	is used
        interpolation - how to interpolate values in imshow
        origin - either 'lower' or 'upper'; determines where imshow should put origin of the image
        """

        if chan not in self.channels:
            print("show_channel: Invalid channel '{}'.".format(chan))
            return

        # only first 'limit' (or less) entries plotted, layloc stores layer indices
        if laylist: layloc = sorted([i for i in set(laylist) if i < self.numlay])[:limit]
        if not laylist or not layloc: layloc = range(min(limit, self.numlay))
        length = len(layloc)

        # determine number of rows and columns
        if length < rowlimit:
            nrows, ncols = 1, length
        else:
            nrows = length // rowlimit
            if length % rowlimit > 0: nrows += 1
            ncols = rowlimit

	    # create figure and axes
        fig, ax = plt.subplots(nrows, ncols, squeeze=False,
                        num='Channel: ' + chan + " [units: " + self.units[chan] + "]")

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
                img = ax[ind].imshow(self.channels[chan][k], aspect=aspect, vmin=vmin, vmax=vmax,
                                interpolation=interpolation, origin=origin)
                ax[ind].set_title(chan + " - layer no. " + str(k % self.numlay))
                
                # add colorbar
                divider = make_axes_locatable(ax[ind])
                cax = divider.append_axes("right", size=0.1, pad=(0.5 if aspect == 'auto' else 0.1))
                plt.colorbar(img, cax=cax)
                
                # increment counter
                l += 1
            except KeyError as err:
                print("show_channel: Key {} not known.".format(err))
                continue

        # improve layout and launch the interactive window
        fig.tight_layout()	
        plt.show()

    def inspect_channels(self, *cl, cutlist=None, cutzlist=None, layer=0, xind=None, yind=None,
        rangelayer='all', linestyle='-', linecolor='w', pointcolor='g', aspect='auto',
        size="30%", interpolation='none'): 
        """show one or more channels together with its cuts for arbitrary layer

        spmdatalist - list of SPMData structures used for plotting
        chanlist - list of channels to be drawn; for each index 'i' channel chanlist[i]
            is shown for SPMData structure spmdatalist[i]
        cutlist - list into which cut control points are stored; if None, then nothing
            is stored, otherwise after each selection is made by mouse and saved by ctrl+p
            new triple [zpos, xpos, ypos] is appended to cutlist so that these data are
            accessible even after the termination of inspect_channels routine
        cutzlist - list into which z-cut points are stored, if None, then nothing is stored,
            otherwise after each selection is made by mouse right-button new double
            [xpos, ypos] is appended to cutzlist so that these data are accesssible
            even after the termination of inspect_channels routine
        layer - which layer of all the channels in chanlist to draw at the beginning, layer
            is identical for all channels in chanlist
        xind, yind - initial position of the intersection point
        linestyle - format specification of the intersecting lines
        linecolor - color specification of the intersecting lines
        pointcolor - color specification of the z-cut points and corresponding lines
        rangelayer - which layer to choose for image min and max determination;
            if 'all', then global minima and maxima are used;
            if 'each', then maxima and minima are determined independently for each layer;
            if rangelayer is not a valid number of a layer, then 'layer' is used instead
        aspect - aspect of the cuts in the figure
        size - relative size of cuts in the figure
        interpolation - how to interpolate values in imshow
        """
        
        spm_inspect_channels([self], cl, cutlist=cutlist, cutzlist=cutzlist, layer=layer, xind=xind, yind=yind,
            rangelayer=rangelayer, linestyle=linestyle, linecolor=linecolor, pointcolor=pointcolor, aspect=aspect,
            size=size, interpolation=interpolation)
        
    # IMPORT, EXPORT

    def copy(self, *cl, ll=None):
        """create a copy of spmdata containing only channels in 'cl'
        and only layers in 'll'
        """
        
        # process lists of channels and layers
        cl = self.channels if not cl else cl
        cl = [chan for chan in cl if chan in self.channels]
        ll = list(range(self.numlay)) if ll is None else ll
        
        # allocate new SPMdata structure
        copyspm = SPMdata()
        
        # copy channels
        for chan, val in self.channels.items():
            if chan not in cl: continue
         
            # laylist will be copyspm.channels_counts[chan]   
            laylist = self.channels_counts[chan]
            laylist = [i for i in laylist if i in ll]
            
            # set new channel to copyspm
            copyspm.setchan(chan, val[ll].copy(), laylist=laylist)
        
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
        
        self.channels_counts = {} # TODO
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
       

