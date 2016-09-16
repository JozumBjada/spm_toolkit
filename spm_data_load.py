from .spm_data_layer import *

"""module containing loading routines for SPMdata structure
"""

def load_data(spmdata, paths=(), suffixes=SUFFIXES, key=None,
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
        three situations may happen:
        if finer=None, then no resampling is performed;
        if finer=True, then data with lower resolution are
        upsampled so that they can be added to data with higher
        resolution;
        if finer=False, then data with higher resolution are
        downsampled so that they can be added to data with lower
        resolution;
        each layer contains attributes xstep and ystep, which
        represent data resolution in x- and y-direction; these
        attributes should not be modified, since they express
        resolution of the original data
    """

    # if 'paths' is None, than loading window is opened,
    # otherwise 'paths' stores filenames
    if not paths:
        root = tkinter.Tk()
        root.withdraw()
        validnames = tkinter.filedialog.askopenfilenames(
            parent=root)
    else:
        validnames = _process_filenames(suffixes, paths)
        
    # check whether filenames are valid
    validnames = _check_filenames(validnames, suffixes)
    
    # if a lot of files to be loaded, warn user that it might
    # take a long time
    if len(validnames) > 15:
        print("load_data: {} files to be loaded.".format(
            len(validnames)))
        
    # if no file to be loaded, halt
    if len(validnames) == 0:
        print("load_data: No file to be loaded.")
        return

	# load individual files
    spmdata.layers = [None]*len(validnames)
    for i, name in enumerate(validnames):
        spmdata.layers[i] = SPMdataLayer(name)

    # create 3D arrays in SPMdata
    print("load_data: Embedding layers...")
    _embed_layers(spmdata, finer, key)

    # set units for channels
    if spmdata.layers:
        for chan in spmdata.channels.keys():

            # get the unit from the first layer
            unit = spmdata.layers[0].units[chan]
            
            # if units from other layers and the same channel
            # differ, let the user know, but do nothing
            for lay in spmdata.layers:
                if unit != lay.units[chan]:
                    print(("load_data: Different units '{}' and "
                        "'{}' for the same channel '{}'."
                        ).format(unit, lay.units[chan], chan))
            spmdata.units[chan] = unit

    # statistics
    namelen = len(validnames)
    print("load_data: In total {} file{} loaded{}".format(
            namelen, "s" if namelen != 1 else "",
            "." if namelen == 0 else ":"))
    for i, name in enumerate(validnames):
        print("\tFile {:2}: {}".format(i + 1, name))

def _embed_layers(spmdata, finer=None, key=None):
    """create 3D-arrays corresponding to each channel by
    retrieving data from all layers, unknown values are
    represented by NANs
    
    finer - if True, then layer data with lower resolution
        are upsampled; if False layer data with higher
        resolution are downsampled; if None, then no
        resampling is performed
    key - key for sorting the layers
    """

    # implicit sorting of layers is according to heights
    key = (lambda x: x.height) if key is None else key

    # sort layers
    spmdata.layers.sort(key=key)

    # resolution in x- and y-direction
    xsteps = [lay.xstep for lay in spmdata.layers]
    ysteps = [lay.ystep for lay in spmdata.layers]

    # original dimensions of data
    xnumsold = [lay.xnum for lay in spmdata.layers]
    ynumsold = [lay.ynum for lay in spmdata.layers]
    
    # do not resample?
    if finer is None:
        # do not care about possibly different data
        # resolution
        print("load_data: Data are not resampled to the same "
            "resolution, so inconsistencies may arise. "
            "Especially xsteps and ysteps now do not correspoind "
            "to real resolution of individual channels.")
        
        # when array does not contain consistent data
        # (in terms of data resolution), then xsteps
        # and ysteps for spmdata are chosen as minima
        # of xsteps and ysteps over all layers
        xstep, ystep = min(xsteps), min(ysteps)
        xnums, ynums = xnumsold.copy(), ynumsold.copy()
    
    # upsample or downsample?
    else:
        if finer:
            xstep, ystep = min(xsteps), min(ysteps)
            print("load_data: Lower-resolution data will be "
                "upsampled.")
        else:
            xstep, ystep = max(xsteps), max(ysteps)        
            print("load_data: Higher-resolution data will be "
                "downsampled.")

        print(("load_data: Original dimensions:\n\txnums = {}"
            "\n\tynums = {}").format(xnumsold, ynumsold))
        
        # adjust array dimensions to contain data with 
        # possibly different resolution        
        xnums = [int(lay.xran / xstep) for lay in spmdata.layers]
        ynums = [int(lay.yran / ystep) for lay in spmdata.layers]
        print(("load_data: New dimensions:\n\txnums = {}"
            "\n\tynums = {}").format(xnums, ynums))
        
    # dimensions of new arrays
    xnum, ynum = max(xnums), max(ynums)

    # for which channels to allocate new arrays
    chanlist = []
    for lay in spmdata.layers:
        chanlist.extend([chan for chan in lay.channels.keys()])

    # save horizontal and vertical resolution for each channel
    spmdata.xsteps = {chan:xstep for chan in set(chanlist)}
    spmdata.ysteps = {chan:ystep for chan in set(chanlist)}

    # update dimensions in layers; this has to be done before the
    # spmdata.set_chan is called and before the loop:
    # for chan in ... etc.
    for i, lay in enumerate(spmdata.layers):
        lay.xnum, lay.ynum = xnums[i], ynums[i]

    # for each channel...
    for chan in set(chanlist):
        # create empty 3D array
        arr = np.empty((spmdata.numlay, xnum, ynum))
        arr.fill(np.nan)

        # fill the array with valid data and update layers
        for i, lay in enumerate(spmdata.layers):
            if chan not in lay.channels.keys(): continue
            
            if  xnums[i] == xnumsold[i] and \
                ynums[i] == ynumsold[i]:
                # if no resampling is necessary, use directly
                # data from the layer
                arrlay = lay.channels[chan]
            else:
                # otherwise resample the data
#                print(("_embed_layers: Data from layer no. {}"
#                    " resampled").format(i))

                xno, yno = lay.channels[chan].shape                               
                arrlay = gen_resample_data_2d(
                    lay.channels[chan].copy(),
                    xno, yno, xnums[i], ynums[i])

            # store data from layers to the array
            arr[i, :xnums[i], :ynums[i]] = arrlay.copy()
                            
        # save the array to the structure and update pointers in
        # layers
        spmdata.set_chan(chan, arr, updatelayers=True)

def _check_filenames(filenames, suffixes):
    """check whether filenames are valid
    """
    
    sufx = tuple(suffixes)

    valfilenames = []
    for name in filenames:
        if os.path.isfile(name) and name.endswith(sufx):
            valfilenames.append(name)
        else:
            print(("load_data: Name '{}' is not a "
                "valid file name.").format(name))
    return valfilenames 

def _process_filenames(suffixes, pathnames='cwd'):
    """process input file names

    pathnames - variable containing file names, if left
        unspecified, then the current working directory
        is searched through
    suffixes - format suffixes which should be loaded;
        if None then files with all appropriate suffixes
        are loaded
    """

    try:
        if pathnames == 'cwd':
            # get files from the current working directory
            pathnamesloc = {os.path.join(os.getcwd(), item)
                for item in os.listdir(os.getcwd())}
            print(("load_data: Data read from the "
               "current working directory: {}").format(
               os.getcwd()))
        elif isinstance(pathnames, str):
            # if pathnames contains a name of a directory,
            # load every file in that directory
            if os.path.isdir(pathnames):
                pathnamesloc = {os.path.join(pathnames, item)
                    for item in os.listdir(pathnames)}
                print(("load_data: Data read from "
                    "directory: {}").format(pathnames))
            # if pathnames contains a name of a file, load that
            # file
            else:
                pathnamesloc = set([pathnames])
        elif all(map(lambda x: isinstance(x, str), pathnames)):
            # if all items in pathnames are strings, get those
            pathnamesloc = set(pathnames)
        else:
            raise TypeError
    except TypeError:
        print("load_data: Invalid data file names.")
        pathnamesloc = []

    return pathnamesloc

