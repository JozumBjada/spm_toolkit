from .spm_data_struct import *

def create_spm(paths=None, numlay=1, finer=True):
    """create SPMdata structure with all necessary attributes
    
    paths - if None, then new SPMdata structure is generated with
        no data; otherwise 'paths' is passed to 'load_data'
        routine
    numlay - number of layers in the new structure; this 
        parameter applies only when 'paths' is None
    finer - if True, then lower-resolution data are upsampled to
        have the same "resolution" as higher-resolution data
    """
    
    spmdata = SPMdata(paths, finer=finer)
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
        # and all elements of spmdata.layers are just pointer to
        # the same SPMdataLayer() object
    
    return spmdata

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
            # for all channels they share, so we copy units from
            # spm2
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
                print("stack_spm: The second SPMdata structure "
                    "does not contain channel '{}'.".format(chan))    
                spmaux = spm1
                ranaux = slice(0, len(spm1))
            else:
                print("stack_spm: The first SPMdata structure "
                    "does not contain channel '{}'.".format(chan))
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
        spmnew.set_chan(chan, arr, updatelayers=False)
        
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

