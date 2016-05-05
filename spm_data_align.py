from spm_data_base import *
from scipy.ndimage import map_coordinates

"""this module implements procedures:
     "align_correlate", 
     "align_forback",
     "align_offsets"
     
     all the other routines in this file are auxiliary routines
"""

class InvalidInputError(Exception): pass

def _showfit_align(numlay, coefs, fac, xn, yn, xoffmin, yoffmin):
    """calculate fit points to be shown in procedure align_correlate
    """

    acoefx, bcoefx, acoefy, bcoefy = coefs

    xfit, yfit = None, None

    if acoefx is not None:
        xfit = acoefx*np.arange(numlay) + bcoefx        
        xfit = np.ceil(xfit - fac * xn / 2)
        xfit -= xoffmin
        xfit[xfit < 0] = 0

    if acoefy is not None:
        yfit = acoefy*np.arange(numlay) + bcoefy
        yfit = np.ceil(yfit - fac * yn / 2)
        yfit -= yoffmin
        yfit[yfit < 0] = 0

    return xfit, yfit
    
def _show_align(numlay, fitll, revactlay, misll,
    xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs):
    """show alignments in procedures align_corr_chan and align_correlate
    """

    # create plot
    fig, ax = plt.subplots(1, 1, num="Extrapolation trends and control points")
    ax.set_title("Extrapolation trends and control points")
    ax.set_xlim([0, numlay])

    # show offsets for all layers
    ax.step(range(numlay), xoffind, "b-", label="x offsets",
        where="mid", linewidth=2.5)
    ax.step(range(numlay), yoffind, "g-", label="y offsets",
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
    if misll:
        # manual setting?
        manstrx = "" if xcoefs is None else " (manual)"
        manstry = "" if ycoefs is None else " (manual)"
                
        # plot reference points used for fitting
        if xcoefs is None:
            ax.plot(fitll, xfit[fitll], "bo",
                label="x control points", markersize=10)
        if ycoefs is None:        
            ax.plot(fitll, yfit[fitll], "go",
                label="y control points", markersize=10)

        # plot inter-/extra-polated points
        ax.plot(misll, xoffind[misll], "bv",
            label="x calculated points", markersize=10)
        ax.plot(misll, yoffind[misll], "gv",
            label="y calculated points", markersize=10)

        # plot fits
        ax.plot(range(numlay), xfit, "b--",
            label="x fit" + manstrx + "\nslope = {:.4}".format(coefs[0]))        
        ax.plot(range(numlay), yfit, "g--",
            label="y fit" + manstry + "\nslope = {:.4}".format(coefs[2]))        
                
    ax.legend(loc='best')
    plt.show()
    

def _preprocess_align(numlay, rchan, actlayi, mislayi, rellayi,
    xcoefs, ycoefs, alay=None, ilay=None):
    """preprocessing stage of procedure align_correlate
    """

    # reference channel
    if rchan is None:
        print("align_correlate: No reference channel specified.")
        raise InvalidInputError

    # active layers
    if actlayi is None:
        actll = range(numlay)
    else:
        actll = {lay % numlay for lay in actlayi}
    actll = sorted(list(actll))                  
                
    if not actll:
        print("align_correlate: No layer to be aligned.")
        raise InvalidInputError
    
    # nonactive layers
    revactlay = [lay for lay in range(numlay) if lay not in actll]
    
    # reference layer for (non)active layers
    alay = actll[-1] if alay is None or alay not in actll else alay
    if revactlay:
        ilay = revactlay[0] if ilay is None or ilay in actll else ilay
    else:
        ilay = None
    print("align_correlate: Active layers: {}".format(actll))
    print("align_correlate: alay, ilay: {}, {}".format(alay, ilay))

    # if there are layers to fit, preprocess misll and fitll
    if mislayi == 'all':
        misll = list(range(numlay))
    else:
        misll = mislayi
    
    fitll = ()
    if misll:
        misll = list({lay % numlay for lay in misll})
        
        # if at least one manual alignment coefficient is not given
        # use fitting by fitll
        if xcoefs is None or ycoefs is None:
            if not rellayi:
                fitll = [lay for lay in range(numlay) if lay not in misll]            
            else:
                fitll = list({lay % numlay for lay in rellayi})
                misll = [lay for lay in misll if lay not in fitll]
            fitll = sorted(fitll)

        misll = sorted(misll)
        
        print("align_correlate: Refer. layers: {}".format(fitll))
        print("align_correlate: Missi. layers: {}".format(misll))

    return actll, revactlay, misll, fitll, alay, ilay


def _corr_align(numlay, rchan, actll, fitll, misll, rlay, arr):
    """cross-correlation stage of procedure align_correlate
    """
    
    print("align_correlate: Alignment according to channel '{}'.".format(rchan))
    relx, rely = np.zeros(numlay), np.zeros(numlay)
    arr1, k = arr[rlay], 1        
    
    if arr1 is None:
        print("align_correlate: Reference layer in reference channel full of NaNs. Halt.")
        raise InvalidInputError

    arr1[np.isnan(arr1)] = np.nanmean(arr1)
    
    # for each layer in actll or fitll find out the most matching offsets
    arlist = sorted([lay for lay in range(numlay) if lay in actll or lay in fitll])
    
    print("align_correlate: Cross-correlation:")
    for i in arlist:
        print("\tLayer {:2} out of {} being processed...".format(k, len(arlist)))
        k += 1
        if misll and i in misll: continue

        # evaluate correlation and get indices of the best correlation match
        arrtemp, arr1tem = arr[i], arr1
        
        if arrtemp is None:
            relx[i], rely[i] = 0, 0
        else:
            arrtemp[np.isnan(arrtemp)] = np.nanmean(arrtemp)
            arr1tem[np.isnan(arrtemp)] = np.nanmean(arr1tem)

            relx[i], rely[i] = get_offinds_corr(arr1tem, arrtemp)

        # OPET DOST BORDEL S TIM, CO JE TEDA X-OVA A Y-OVA OSA...
    
    return relx, rely

def _fit_align(xcoefs, ycoefs, misll, fitll, relx, rely):
    """fitting stage of procedure align_correlate
    """

    acoefx, bcoefx = None, None
    acoefy, bcoefy = None, None
    xfit, yfit = None, None

    # inter-/extra-polation
    if misll:
        # find linear fit
        from scipy.stats import linregress
            
        if xcoefs is None:
            acoefx, bcoefx = linregress(fitll, relx[fitll])[:2]                
        else:
            acoefx, bcoefx = xcoefs

        if ycoefs is None:
            acoefy, bcoefy = linregress(fitll, rely[fitll])[:2]                
        else:
            acoefy, bcoefy = ycoefs
   
        # calculate offsets for layers in misll
        for i in misll:
            relx[i] = acoefx*i + bcoefx
            rely[i] = acoefy*i + bcoefy

    return relx, rely, acoefx, bcoefx, acoefy, bcoefy
    
def _postprocess_align(spmdata, rchan, relx, rely, fac, revactlay, actll, alay, ilay):
    """postprocessing stage of procedure align_correlate
    """
    
    xn, yn = np.array(spmdata.xnums), np.array(spmdata.ynums)
    xo = np.array([lay.xoffind[rchan] if rchan in lay.channels.keys() else 0
        for lay in spmdata.layers])
    yo = np.array([lay.yoffind[rchan] if rchan in lay.channels.keys() else 0
        for lay in spmdata.layers])

    # offsets        
    relx = np.round(relx - fac * xn / 2)
    rely = np.round(rely - fac * yn / 2)

    # reference offsets with respect to the array origin
    xoffind = relx
    yoffind = rely
    
    # for nonactive layers use the original offsets
    if revactlay:
        print("align_correlate: Nonactive layer offsets may be affected.")
        # now we set xoffind and yoffind of all channels to values of
        # xoffind and yoffind of the rchan channel
        xoffind[revactlay] = fac * xo[revactlay] + xoffind[alay] - fac * xo[ilay]
        yoffind[revactlay] = fac * yo[revactlay] + yoffind[alay] - fac * yo[ilay]
        
    # all offsets with respect to the array origin
    xoffmin, yoffmin = xoffind.min(), yoffind.min()
    xoffind -= xoffmin
    yoffind -= yoffmin

    # xoffmin and yoffmin are used in demonstrating plot
    xoffmin += relx.min()
    yoffmin += rely.min()

    return xn, yn, xoffind, yoffind, xoffmin, yoffmin


def _upsample_align(numlay, layers, rchan, fac, chanlist, order):
    """upsample stage of procedure align_correlate
    """
    
    arrdict = {}
    if rchan not in chanlist: chanlist.append(rchan)

    # upsampling of channels in chanlist
    # assume that xnum and ynum is equal for all channels to xnum
    # and ynum of the reference channel for each layer
    for chan in chanlist:
        arrtmp = [None]*numlay
        
        # for fac = 1 no upsampling is done
        if fac > 1:
            print("align_correlate: Channel '{}' being upsampled.".format(chan)) 
            for i, lay in enumerate(layers):
                if chan not in lay.channels.keys(): continue
                
                xno, yno = lay.channels[chan].shape                
                x = np.linspace(0, xno - 1, fac * xno)
                y = np.linspace(0, yno - 1, fac * yno)
                coords = np.meshgrid(x, y, indexing='ij')
                arrtmp[i] = map_coordinates(lay.channels[chan],
                    coords, order=order, cval=np.nan)                
        else:
            for i, lay in enumerate(layers):
                if chan not in lay.channels.keys(): continue
                arrtmp[i] = lay.channels[chan]
        
        arrdict[chan] = arrtmp

    return arrdict[rchan], arrdict
    
   
def _reallocate_align(spmdata, chanlist, xoffind, yoffind,
    fac, xn, yn, arrdict, downsample, order=0):
    """reallocating stage of procedure align_correlate
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
        aux = np.empty((spmdata.numlay, xlim, ylim), dtype=float)    
        aux.fill(np.nan)

        # fill new array with valid data and appropriate offsets
        for i, lay in enumerate(spmdata.layers):      
            if arrdict[chan][i] is None:
                aux2 = np.empty((lay.xnum, lay.ynum))
                aux2.fill(np.nan)
            else:
                aux2 = arrdict[chan][i]
                
            aux[i, 
                xoffind[i]:xoffind[i] + xn[i],
                yoffind[i]:yoffind[i] + yn[i]
                ] = aux2.copy()

        # save channel
        auxlist[chan] = aux.copy()

    # downsampling and updating data attributes...
    if not downsample: fac = 1

    # save nums and offinds for each layer
    for i, lay in enumerate(spmdata.layers):                
        lay.xnum = int(xn[i] / fac)
        lay.ynum = int(yn[i] / fac)
            
    for chan in chanlist:
        aux = auxlist[chan]
                
        for i, lay in enumerate(spmdata.layers):                
            if chan not in lay.channels.keys(): continue
            lay.xoffind[chan] = math.floor(xoffind[i] / fac)
            lay.yoffind[chan] = math.floor(yoffind[i] / fac)

        # downsampling                    
        if downsample and fac > 1:            
            print(("align_correlate: Channel '{}' "
                "being downsampled.").format(chan))
            x = np.linspace(0, aux.shape[1],
                math.ceil(aux.shape[1] / fac), endpoint=False)
            y = np.linspace(0, aux.shape[2],
                math.ceil(aux.shape[2] / fac), endpoint=False)
            z = range(spmdata.numlay)
            coords = np.meshgrid(z, x, y, indexing='ij')

            # ZDE JE DALSI PROBLEM V TOM, ZE NAN S CIMKOLI DA OPET
            # NAN, TJ. PRI INTERPOLACI
            # SE DO AUX NA MISTA, KDE BY MOHLA BYT PLATNA DATA, 
            # DOSTANE NAN A ROZMERY VYSLEDNE
            # OBLASTI PLATNYCH DAT JSOU PO DOWNSAMPLINGU OBECNE 
            # MENSI NEZ PUVODNI ROZMERY
            # ULOZENE V spmdata.layers[i].xnum A spmdata.layers[i].ynum
            
            # during interpolation in map_coordinates anything
            # with NaN is again NaN, so it may be the case that
            # valid data are truncated due to padding NaNs, to
            # circumvent this we replace NaNs with an average
            # over aux
            aux[np.isnan(aux)] = np.nanmean(aux)
            
            # aux now stores downsampled 3D array for channel
            # 'chan', but instead of NaNs,
            # np.nanmean(aux) is used as a padding
            aux = map_coordinates(aux, coords,
                order=order, mode='nearest')#cval=np.nan)
        
            # store tentative aux in auxlist
            auxlist[chan] = aux.copy()
                
            # create new aux filled with NaNs
            aux = np.empty_like(auxlist[chan])
            aux.fill(np.nan)

            # to aux insert valid data from auxlist
            for i, lay in enumerate(spmdata.layers):
                if chan not in lay.channels.keys(): continue
                aux[i,
                    lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                    lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum
                    ] = auxlist[chan][i,
                    lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                    lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum                    
                    ].copy()
                
        # update channel
        spmdata.setchan(chan, aux.copy())                    

        # for each layer update pointer to valid data
        for i, lay in enumerate(spmdata.layers):
            if chan not in lay.channels.keys(): continue
            lay.setchan(chan, spmdata.channels[chan][i,
                lay.xoffind[chan]:lay.xoffind[chan] + lay.xnum,
                lay.yoffind[chan]:lay.yoffind[chan] + lay.ynum])

def align_correlate(spmdata, *cl, rchan=None, rlay=0, alay=None,
    ilay=None, actll=None, fitll=None, misll=None, show=False,
    xcoefs=None, ycoefs=None, fac=1, downsample=True, order=1):
    """align arrays according to rchan channel, i.e. offsets are
    calculated only once for rchan and then used multiple times
    for all channels in cl, unknown values are represented by NANs,
    resulting array may have bigger dimensions than original arrays
    
    cl - list of channels to be aligned, if empty, then all
        channels are aligned
    rlay - index of a fixed reference layer with which all layers
        are compared
    alay - if actll is not None or not all layers, then alay
        specifies to which layer inactive layers should adhere,
        i.e. offind of ilay layer will be set to the offind of alay;
        if alay is None or lies outside actll, then it is set
        to be the last layer in actll
    ilay - analog of alay for inactive layers, if None or in 
        actll, then it is set to be the first layer not in actll
    rchan - reference channel; all channels in cl are aligned
        according to channel 'rchan'
    actll - list of 'active' layers, which are to be aligned;
        if None, then all layers are used
    fitll - list of layers which are used to determine by fitting
        offsets for layers in misll, if None or empty, then
        fitll is effectively put equal to all layers not lying in misll
    misll - list of 'missing' layers for which offsets should
        be calculated by inter-/extra-polation;
        if None or empty, then correlation calculation is performed
        for all layers and no inter-/extra-polation is used;
        if misll and fitll are not disjoint, then misll is
        put equal to misll minus fitll;
        if misll='all', then all layers are used for 
        inter-/extra-polation
    show - if True, then plot depicting fitll points and
        extrapolation trends is shown,
        control points - fitll points
        calculated points - misll points
        x fit, y fit - fit of x-offsets and y-offsets, respectively
        (manual) - shown if fit is set manually via xcoefs or ycoefs
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
    downsample - whether to perform downsampling after correlation
    order - order of interpolation for resampling
    
    BEWARE - there might be problems for interpolation order > 1;
        due to round off errors the spmdata.layers[i][chan] arrays
        may contain NaNs with which map_coordinates cannot
        work for order > 1; this issue takes place in '_reallocate_align'
        routine in downsampling stage, specifically in command
        np.linspace in its third parameter
 
        
    Example 1: I have given SPMdata structure 'dta' with 30 layers. I want to align
        channels CURR and PHAS. As a reference channel I want to choose FREQ,
        as a reference layer I choose layer no. 3. That is, I write:
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, rlay=3)'.
        But, from the result I see that cross-correlation does not work well for layers
        no. 10, 11, 12, 15, 16. I would like to use linear fit for offsets of these 
        layers. So I append 'misll=[10, 11, 12, 15, 16]' to argument list.
        New result I obtain nevertheless does not fit offsets for these layers well.
        Reason for this may be that offsets for layers outside 'misll' vary a lot, so I
        would like to calculate a linear fit only from offsets pertaining to layers
        no. 5, 6, 7. Thus I append 'fitll=[5, 6, 7]' to argument list and also set 
        'show' to True to see how offsets are calculated. At the end I have
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, rlay=3, misll=[10, 11, 12, 15, 16],
            fitll=[5, 6, 7])'.
            
    Example 2: Suppose I have given SPMdata structure as in Example 1. But now 
        I want to align according to layer no. 0, which is default choice. Now I would
        like to align only layers no. 10 through no. 29 and leave layers no. 0 through
        no. 5 unaffected (i.e. layers no. 0 through no. 5 of channel 'rchan'! unaffected,
        since routine 'align_corr_chan' aligns channels according to channel 'rchan', all
        channels in 'cl' different from 'rchan' may get different offsets even
        for layers no. 0 through no. 5). So I write:
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, actll=range(6, 30), show=True,
            misll=range(10, 30))'.
        But, the resulting alignment is not good, so I read actual offsets used by
        routine from the plot shown and set fit coefficients for x-offsets by hand,
        so I write
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, actll=range(6, 30), show=True,
            misll=range(10, 30), xcoefs=(5., 10.))'.        
    """

    chanlist = list(spmdata.channels if not cl else cl)
    chanlist = [chan for chan in chanlist if chan in spmdata.channels]
    
    # since each layer contains attributes xnum and ynum, i.e. dimensions of valid
    # data area, which are identical for all channels, let us exclude possibility,
    # that 'align_correlate' changes dimension of valid data area for different
    # channels
    if not downsample and fac > 1 and len(chanlist) < len(spmdata.channels):
        print("align_correlate: Downsampling on, but not all channels included.")
        return

    numlay, layers = spmdata.numlay, spmdata.layers

    # layer lists initialization
    try:
        actll, revactlay, misll, fitll, alay, ilay = _preprocess_align(
            numlay, rchan, actll, misll, fitll, xcoefs, ycoefs, alay, ilay)
    except InvalidInputError:
        return

    # upsampling        
    arr, arrdict = _upsample_align(numlay, layers, rchan, fac, chanlist, order)
   
    # cross-correlation     
    try:    
        relx, rely = _corr_align(numlay, rchan, actll, fitll, misll, rlay, arr)
    except InvalidInputError:
        return

    # inter-/extra-polation
    relx, rely, *coefs = _fit_align(xcoefs, ycoefs, misll, fitll, relx, rely)

    # postprocessing
    xn, yn, xoffind, yoffind, xoffmin, yoffmin = _postprocess_align(
        spmdata, rchan, relx, rely, fac, revactlay, actll, alay, ilay)

    # plotting - show trends and control points
    if show:
        xfit, yfit = _showfit_align(numlay, coefs, fac, xn, yn, xoffmin, yoffmin)
        _show_align(numlay, fitll, revactlay, misll,
                    xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs)

    # reallocation of arrays according to offsets and downsampling
    _reallocate_align(spmdata, chanlist, xoffind, yoffind, fac,
        xn, yn, arrdict, downsample, order)

    print("align_correlate: Channels {} aligned.".format(chanlist))            

















# JESTLI TO SPRAVNE CHAPU, TAK align_trend SE SNAZI OBEJIT RESAMPLOVANI KANALU, KTERE NEJSOU REFERENCNI KVULI USPORE VYPOCETNIHO CASU


#def _showfit_align_trend(self, coefs, fac, xn, yn, xoffmin, yoffmin):
#    """calculate fit points to be shown in procedure align_trend
#    """

#    acoefx, bcoefx, acoefy, bcoefy = coefs
#    xfit, yfit = None, None
#    
#    if acoefx is not None:
#        xfit = acoefx*np.arange(self.numlay) + bcoefx
#        xfit = xfit - fac * xn / 2
#        xfit -= xoffmin
#        xfit = np.round(xfit / fac)
#        xfit[xfit < 0] = 0
#    
#    if acoefy is not None:
#        yfit = acoefy*np.arange(self.numlay) + bcoefy
#        yfit = yfit - fac * yn / 2
#        yfit -= yoffmin
#        yfit = np.round(yfit / fac)
#        yfit[yfit < 0] = 0

#    return xfit, yfit

#def _reallocate_align_trend(self, chanlist, xoffind, yoffind, fac, xn, yn, arrdict, order=0):
#    """reallocating stage of procedure align_trend
#    """

#    # exact offinds
#    xoflt, yoflt = xoffind, yoffind
#    
#    # integral offinds
#    xoffind, yoffind = np.round(xoflt).astype(int), np.round(yoflt).astype(int)
##        xoffind, yoffind = np.floor(xoflt).astype(int), np.floor(yoflt).astype(int)

#    # dimensions of the new array
#    xlim = max(xoffind + xn) - min(xoffind)
#    ylim = max(yoffind + yn) - min(yoffind)

#    print("xoffind - xoflt: ", xoffind - xoflt)

#    # align each channel...
#    for chan in chanlist:
#        # allocate new array for the channel 'chan'
#        aux = np.empty((self.numlay, xlim, ylim), dtype=float)
#        aux.fill(np.nan)

#        # fill new array with valid data and appropriate offsets            
#        for i, lay in enumerate(self.layers):
#            x = (xoffind[i] - xoflt[i]) + np.arange(lay.xnum)
#            y = (yoffind[i] - yoflt[i]) + np.arange(lay.ynum)
#            coords = np.meshgrid(x, y, indexing='ij')
#            aux[i, 
#                xoffind[i]:xoffind[i] + xn[i],
#                yoffind[i]:yoffind[i] + yn[i]
#                ] = map_coordinates(
#                    lay.channels[chan], coords,
#                    order=order, mode='nearest')

#        # store channel
#        self.setchan(chan, aux.copy())
#        
#        # update offinds and channels for each layer
#        for i, lay in enumerate(self.layers):
#            lay.xoffind[chan] = xoffind[i]
#            lay.yoffind[chan] = yoffind[i]

#            lay.channels[chan] = self.channels[chan][i, 
#                xoffind[i]:xoffind[i] + xn[i],
#                yoffind[i]:yoffind[i] + yn[i]]

#def _upsample_align_trend(self, rchan, fac, order):
#    """upsample stage of procedure align_trend
#    """
#    
#    # upsampling of reference channel
#    arr = [None]*self.numlay

#    print("align_correlate: Channel '{}' being upsampled.".format(rchan))
#    for i, lay in enumerate(self.layers):        
#        xno, yno = lay.channels[rchan].shape                
#        x = np.linspace(0, xno - 1, fac * xno)
#        y = np.linspace(0, yno - 1, fac * yno)
#        coords = np.meshgrid(x, y, indexing='ij')
#        arr[i] = map_coordinates(lay.channels[rchan], coords, order=order)#, cval=np.nan) # cval!!!

#    return arr

#def _postprocess_align_trend(self, rchan, relx, rely, fac, revactlay, actll):
#    """postprocessing stage of procedure align_trend
#    """
#    
#    xn, yn = np.array(self.xnums), np.array(self.ynums)
#    xo = np.array([lay.xoffind[rchan] for lay in self.layers])
#    yo = np.array([lay.yoffind[rchan] for lay in self.layers])

#    # offsets        
#    relx = relx - fac * xn / 2
#    rely = rely - fac * yn / 2

#    relx[revactlay] = relx[actll].min()
#    rely[revactlay] = rely[actll].min()

#    # reference offsets with respect to the array origin
#    xoffind = relx - relx.min()
#    yoffind = rely - rely.min()
#    
#    # for nonactive layers use the original offsets
#    if revactlay:
#        print("align_correlate: Nonactive layer offsets may be affected.")
#        # now we set xoffind and yoffind of all channels to values of
#        # xoffind and yoffind of the rchan channel
#        xoffind[revactlay] = fac * xo[revactlay]
#        yoffind[revactlay] = fac * yo[revactlay]
#        
#    # all offsets with respect to the array origin
#    xoffmin, yoffmin = xoffind.min(), yoffind.min()
#    xoffind = (xoffind - xoffmin).astype(float) / fac
#    yoffind = (yoffind - yoffmin).astype(float) / fac

#    xoffmin += relx.min()
#    yoffmin += rely.min()

#    return xn, yn, xoffind, yoffind, xoffmin, yoffmin
#                
#def align_trend(self, *cl, rchan=None, rlay=0, actll=None, fitll=None, misll=None, show=False,
#    xcoefs=None, ycoefs=None, fac=1, order=1):
#    """align arrays according to rchan channel, i.e. offsets are calculated only once for rchan and then used
#    multiple times for all channels in cl, unknown values are represented by NANs,
#    resulting array may have bigger dimensions than original arrays
#    
#    cl - list of channels to be aligned, if empty, then all channels are aligned
#    rlay - all layers are compared to fixed layer with index 'lay'
#    rchan - all channels in cl are aligned according to channel rchan
#    actll - list of 'active' layers, which are to be aligned, if None, then all layers are used
#    fitll - list of 'relevant' layers which are used to determine offsets for layers in misll,
#        if None or empty, then fitll is effectively put equal to all layers not lying in misll
#    misll - list of 'missing' layers for which offsets should be calculated by inter-/extra-polation,
#        if None or empty, then correlation calculation is performed for all layers and no 
#        inter-/extra-polation is used, if misll and fitll are not disjoint, then misll is
#        put equal to misll minus fitll, if misll='all', then all layers are used for
#        inter-/extra-polation
#    show - if True, then plot depicting fitll points and extrapolation trends is shown,
#        control points - fitll points
#        calculated points - misll points
#        x fit, y fit - fit of x-offsets and y-offsets, respectively
#        (manual) - shown if fit is set manually via xcoefs or ycoefs
#        x offsets, y offsets - resulting x-offsets and y-offsets, respectively
#        red area - for this area no alignment is active, i.e. corresponding layers do not
#            lay in actll
#    xcoefs, ycoefs - linear and absolute coefficients for manual interpolation for x and y
#        dimension, respectively; if not None, then these parameters are used instead of
#        automatically calculated ones, in which case xcoefs=(xlincoef, xabscoef) etc.;
#        if not None, then fitll is not used
#    fac - upsampling factor
#    order - order of interpolation for resampling
#    
#    BEWARE - there might be problems for interpolation order > 1; due to round off errors
#        the self.layers[i][chan] arrays may contain NaNs with which map_coordinates cannot
#        work for order > 1; this issue takes place in '_reallocate_align' routine in
#        downsampling stage, specifically in command np.linspace in its third parameter
# 
#        
#    Example 1: I have given SPMdata structure 'dta' with 30 layers. I want to align
#        channels CURR and PHAS. As a reference channel I want to choose FREQ,
#        as a reference layer I choose layer no. 3. That is, I write:
#        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, rlay=3)'.
#        But, from the result I see that cross-correlation does not work well for layers
#        no. 10, 11, 12, 15, 16. I would like to use linear fit for offsets of these 
#        layers. So I append 'misll=[10, 11, 12, 15, 16]' to argument list.
#        New result I obtain nevertheless does not fit offsets for these layers well.
#        Reason for this may be that offsets for layers outside 'misll' vary a lot, so I
#        would like to calculate a linear fit only from offsets pertaining to layers
#        no. 5, 6, 7. Thus I append 'fitll=[5, 6, 7]' to argument list and also set 
#        'show' to True to see how offsets are calculated. At the end I have
#        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, rlay=3, misll=[10, 11, 12, 15, 16],
#            fitll=[5, 6, 7])'.
#            
#    Example 2: Suppose I have given SPMdata structure as in Example 1. But now 
#        I want to align according to layer no. 0, which is default choice. Now I would
#        like to align only layers no. 10 through no. 29 and leave layers no. 0 through
#        no. 5 unaffected (i.e. layers no. 0 through no. 5 of channel 'rchan'! unaffected,
#        since routine 'align_corr_chan' aligns channels according to channel 'rchan', all
#        channels in 'cl' different from 'rchan' may get different offsets even
#        for layers no. 0 through no. 5). So I write:
#        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, actll=range(6, 30), show=True,
#            misll=range(10, 30))'.
#        But, the resulting alignment is not good, so I read actual offsets used by
#        routine from the plot shown and set fit coefficients for x-offsets by hand,
#        so I write
#        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, actll=range(6, 30), show=True,
#            misll=range(10, 30), xcoefs=(5., 10.))'.        
#    """

#    chanlist = list(self.channels if not cl else cl)
#    chanlist = [chan for chan in chanlist if chan in self.channels]
#    
#    # layer lists initialization
#    actll, revactlay, misll, fitll = self._preprocess_align(
#                                rchan, actll, misll, fitll, xcoefs, ycoefs)

#    # upsampling        
#    arr = self._upsample_align_trend(rchan, fac, order)
#   
#    # cross-correlation         
#    relx, rely = self._corr_align(rchan, actll, fitll, misll, rlay, arr)

#    # inter-/extra-polation
#    relx, rely, *coefs = self._fit_align(xcoefs, ycoefs, misll, fitll, relx, rely)

#    # postprocessing
#    xn, yn, xoffind, yoffind, xoffmin, yoffmin = self._postprocess_align_trend(
#                                rchan, relx, rely, fac, revactlay, actll)

#    # plotting - show trends and control points
#    if show:
#        xfit, yfit = self._showfit_align_trend(coefs, fac, xn, yn, xoffmin, yoffmin)
#        self._show_align(fitll, revactlay, misll,
#                    xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs)

#    # reallocation of arrays according to offsets
#    self._reallocate_align_trend(chanlist, xoffind, yoffind, fac,
#        xn, yn, order)

#    print("align_trend: Channels {} aligned.".format(chanlist))  



# function determining total offsets
def _get_shifts(layers, lay, chanf, chanb, xdir, ydir, same):
    """get total x-shifts and y-shifts
    """

#        # check presence of channels
#        if self.channels_fb[chanf][0] is None:
#            print("align_forback: Forward channel '{}' is missing.".format(chanf))
#            return None, None

#        if self.channels_fb[chanf][1] is None:
#            print("align_forback: Backward channel '{}' is missing.".format(chanb))
#            return None, None

    # evaluate correlation and get indices of the best correlation match
    layf = layers[lay].channels[chanf]
    layb = layers[lay].channels[chanb]
    valf, valb = np.nanmean(layf), np.nanmean(layb) 
    layf[np.isnan(layf)] = valf
    layb[np.isnan(layb)] = valb
        
    relx, rely = get_offinds_corr(layf, layb)
    
    # OPET DOST BORDEL S TIM, CO JE TEDA X-OVA A Y-OVA OSA...
    xdir, ydir = ydir, xdir # aby to vychazelo, tak takhle, ale WTF???!!!
        
    # calculate offsets            
    if not xdir:
        relx = 0
    else:
        relx -= int(layers[lay].xnum / 2)
        if same:
            relx += layers[lay].xoffind[chanf]
            relx -= layers[lay].xoffind[chanb]

    if not ydir:
        rely = 0
    else:
        rely -= int(layers[lay].ynum / 2)
        if same:
            rely += layers[lay].yoffind[chanf]
            rely -= layers[lay].yoffind[chanb]

    return relx, rely

def align_forback(spmdata, *cl, rchan=None, same=True, rlay=0, xdir=True, ydir=True, order=0):
    """align all channels in 'cl' according to 'rchan' channel so that forward
    and backward scans are moved towards each other and aligned 'to
    the centre', this routine does not align layers within one channel, but
    channels among themselves instead
    
    cl - sequence of channels to be aligned, if empty, than all
        channels are aligned
    rchan - reference channel; if None, then each channel is aligned
        independently
    same - if True, then resulting 3D arrays for all channels in cl are of the
        same dimensions; if False, then 3D arrays for forward and backward
        direction of each physical quantity have the same dimension
    rlay - reference layer 
    xdir - if True, than allow shift in x-direction
    ydir - if True, than allow shift in y-direction
    order - order of interpolation during shifting
    """

    # preprocess the list of channels to align, let cl contain only forward scans
    claux = [chan for chan in spmdata.channels.keys() if not is_backward_scan(chan)]
    claux = [chan for chan in claux if get_backward_scan(chan) in spmdata.channels.keys()]
    if cl:
        claux = [chan for chan in claux if chan in cl]
    cl = claux
    print("align_forback: Channels to be aligned: {}.".format(cl))
    
    vecx,  vecy  = [None]*len(cl), [None]*len(cl)
    xlimf, xlimb = [None]*len(cl), [None]*len(cl)
    ylimf, ylimb = [None]*len(cl), [None]*len(cl)
#    xnum,  ynum  = [None]*len(cl), [None]*len(cl)
    fxoffs, bxoffs = [None]*len(cl), [None]*len(cl)
    fyoffs, byoffs = [None]*len(cl), [None]*len(cl)

    # process reference channel
    if rchan is None:
        pass
    elif rchan in spmdata.channels.keys() and not is_backward_scan(rchan):
        ichanb = get_backward_scan(rchan)
        if ichanb not in spmdata.channels.keys():
            print("align_forback: Reference channel has no backward counterpart.")
            return
        else:                
            vecxi, vecyi = _get_shifts(spmdata.layers, rlay, rchan, ichanb, xdir, ydir, same)
    else:
        print("align_forback: Invalid reference channel.")
        return

    # find alignments for each pair of forw./backw. channels
    for i, chanf in enumerate(cl):
        chanb = get_backward_scan(chanf)
        
        if rchan is None:
            vecx[i], vecy[i] = _get_shifts(spmdata.layers, rlay, chanf, chanb, xdir, ydir, same)            
        else:
            vecx[i], vecy[i] = vecxi, vecyi

        # if xlimf != xlimb or ylimf != ylimb, then
        # resulting 3D arrays are not ensured to be
        # properly aligned w.r.t. forward/backward!
        xlimf[i], ylimf[i] = spmdata.channels[chanf].shape[1:]
        xlimb[i], ylimb[i] = spmdata.channels[chanb].shape[1:]
        fxoffs[i] = spmdata.layers[rlay].xoffind[chanf] 
        bxoffs[i] = spmdata.layers[rlay].xoffind[chanb] 
        fyoffs[i] = spmdata.layers[rlay].yoffind[chanf] 
        byoffs[i] = spmdata.layers[rlay].yoffind[chanb] 
#        xnum[i],  ynum[i]  = self.layers[rlay].xnum, self.layers[rlay].ynum
    
    print("ylimf: ", ylimf)
    
    # get dimensions of new arrays
    xlim = [int(np.abs(vi).max()) + max(lfi, lbi) for vi, lfi, lbi in zip(vecx, xlimf, xlimb)]
    ylim = [int(np.abs(vi).max()) + max(lfi, lbi) for vi, lfi, lbi in zip(vecy, ylimf, ylimb)]
     
    # if same == True, then modify vecxf, vecxb, ... so that all channels
    # have the same array dimension
    if same:
        xlv, ylv = max(xlim), max(ylim)
        xlim, ylim = [xlv]*len(cl), [ylv]*len(cl) 
        xnum,  ynum  = spmdata.layers[rlay].xnum, spmdata.layers[rlay].ynum

        #  TODO

        vecxf = [int((xlv - xnum - v)/2) - off for v, f, off in zip(vecx, xlimf, fxoffs)] 
        vecxb = [int((xlv - xnum + v)/2) - off for v, b, off in zip(vecx, xlimb, bxoffs)] 
        vecyf = [int((ylv - ynum - v)/2) - off for v, f, off in zip(vecy, ylimf, fyoffs)] 
        vecyb = [int((ylv - ynum + v)/2) - off for v, b, off in zip(vecy, ylimb, byoffs)] 

#        vecyf = [int((ylv - f - v)/2) for v, f in zip(vecy, ylimf)] 
#        vecyb = [int((ylv - b + v)/2) for v, b in zip(vecy, ylimb)]

#        vecxf = [int((xlv - f - v)/2) for v, f in zip(vecx, xlimf)] 
#        vecxb = [int((xlv - b + v)/2) for v, b in zip(vecx, xlimb)] 
#        vecyf = [int((ylv - f - v)/2) for v, f in zip(vecy, ylimf)] 
#        vecyb = [int((ylv - b + v)/2) for v, b in zip(vecy, ylimb)]
    else:
        vecxf = [0 if v > 0 else int(-v) for v in vecx]
        vecxb = [0 if v < 0 else int(+v) for v in vecx]
        vecyf = [0 if v > 0 else int(-v) for v in vecy]
        vecyb = [0 if v < 0 else int(+v) for v in vecy]
        
    print("vecyf: ", vecyf)
        
    # update 3D arrays
    for i, chanf in enumerate(cl):
        chanb = get_backward_scan(chanf)
        
        # create auxiliary arrays
        auxf = np.empty((spmdata.numlay, xlim[i], ylim[i]))
        auxb = np.empty((spmdata.numlay, xlim[i], ylim[i]))
        auxf.fill(np.nan)
        auxb.fill(np.nan)
        
        # forward channel - fill array with data
        xoff, yoff = vecxf[i], vecyf[i]
        auxf[:, xoff:xoff + xlimf[i], yoff:yoff + ylimf[i]] = spmdata.channels[chanf]

        # backward channel - fill array with data
        xoff, yoff = vecxb[i], vecyb[i]
        auxb[:, xoff:xoff + xlimb[i], yoff:yoff + ylimb[i]] = spmdata.channels[chanb]

        # update channels
        spmdata.setchan(chanf, auxf)
        spmdata.setchan(chanb, auxb)
        
        # update layer pointers and xoffinds, yoffinds
        for j, lay in enumerate(spmdata.layers):            
            xnum, ynum = lay.xnum, lay.ynum
            
            # forward channel
            if chanf in lay.channels.keys():
                xoff, yoff = lay.xoffind[chanf] + vecxf[i], lay.yoffind[chanf] + vecyf[i]
                lay.xoffind[chanf], lay.yoffind[chanf] = xoff, yoff                
                lay.channels[chanf] = spmdata.channels[chanf][j][xoff:xoff + xnum, yoff:yoff + ynum]
#                print("for: xoff, yoff: ", xoff, yoff)

            # backward channel
            if chanb in lay.channels.keys():
                xoff, yoff = lay.xoffind[chanb] + vecxb[i], lay.yoffind[chanb] + vecyb[i]
                lay.xoffind[chanb], lay.yoffind[chanb] = xoff, yoff
                lay.channels[chanb] = spmdata.channels[chanb][j][xoff:xoff + xnum, yoff:yoff + ynum]
#                print("bck: xoff, yoff: ", xoff, yoff)



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

    # process list of channels to align
    cl = self.channels.keys() if not cl else cl
    cl = [chan for chan in self.channels.keys() if chan in cl]

    yoffsets, xoffsets = self.xoffsets, self.yoffsets
#    xoffsets, yoffsets = self.xoffsets, self.yoffsets
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

    for chan in cl:
        # create new array
        aux = np.empty((self.numlay, xdimnew, ydimnew), dtype=float)
        aux.fill(np.nan)

        # update lay.xoffind and lay.yoffind
        for i, lay in enumerate(self.layers):
            if chan not in lay.channels.keys(): continue
            lay.xoffind[chan] = int(xoffind[i])
            lay.yoffind[chan] = int(yoffind[i])

            aux[i, 
                xoffind[i]:xoffind[i] + xn[i],
                yoffind[i]:yoffind[i] + yn[i]
                ] = lay.channels[chan].copy()

        self.setchan(chan, aux.copy())

        for i, lay in enumerate(self.layers):     
            if chan not in lay.channels.keys(): continue       
            lay.setchan(chan,
                self.channels[chan][i,
                    xoffind[i]:xoffind[i] + xn[i],
                    yoffind[i]:yoffind[i] + yn[i]])
        

