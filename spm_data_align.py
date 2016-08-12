from .spm_data_base import *
from scipy.ndimage import map_coordinates

"""this module implements procedures:
     "align_correlate", 
     "align_forback",
     "align_offsets"
     
     all the other routines in this file are auxiliary routines
"""

# TODO: Upsampling and downsampling might be done easily by
# already implemented routine 'zoom'
# in scipy.ndimage.interpolation package!!!

class InvalidInputError(Exception): pass


def _get_offinds_corr(arr1, arr2):            
    """return argmax indices of cross-correlation array, which is
    calculated by fft
    """
    corrarr = np.fft.ifftn(
        np.fft.fftn(arr1) * np.conj(np.fft.fftn(arr2))
        )        
    corrarr = np.roll(corrarr,
        int(np.floor(corrarr.shape[0] / 2)), axis=0)
    corrarr = np.roll(corrarr,
        int(np.floor(corrarr.shape[1] / 2)), axis=1)
    relx, rely = np.unravel_index(
        np.argmax(corrarr), corrarr.shape)

    return relx, rely

def _get_offinds_corr_3d(arr1, arr2):            
    """return argmax indices of cross-correlation array, which is
    calculated by fft
    """
    corrarr = np.fft.ifftn(
        np.fft.fftn(arr1) * np.conj(np.fft.fftn(arr2))
        )        

    corrarr = np.roll(corrarr,
        int(np.floor(corrarr.shape[0] / 2)), axis=0)
    corrarr = np.roll(corrarr,
        int(np.floor(corrarr.shape[1] / 2)), axis=1)
    corrarr = np.roll(corrarr,
        int(np.floor(corrarr.shape[2] / 2)), axis=2)

    corrarr = corrarr.real

#    TOLERANCY = 1e-9
#    corrarr.real[np.abs(corrarr.real) < TOLERANCY] = 0
#    corrarr.imag[np.abs(corrarr.imag) < TOLERANCY] = 0

    relx, rely, relz = np.unravel_index(
        np.argmax(corrarr), corrarr.shape)

    return relx, rely, relz, corrarr

def _showfit_align(numlay, coefs, fac, xn, yn, xoffmin, yoffmin):
    """calculate fit points to be shown in procedure
    align_correlate
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
    """show alignments in procedures align_corr_chan and
    align_correlate
    """

    # create plot
    fig, ax = plt.subplots(1, 1,
        num="Extrapolation trends and control points")
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

    # if there are layers whose offsets are to be
    # inter-/extra-polated, plot control points,
    # calculated points and fits
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
        ax.plot(range(numlay), xfit, "b--", label=("x fit" + \
            manstrx + "\nslope = {:.4}".format(coefs[0])))        
        ax.plot(range(numlay), yfit, "g--", label=("y fit" + \
            manstry + "\nslope = {:.4}".format(coefs[2])))
                
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
    revactlay = [lay
        for lay in range(numlay) if lay not in actll]
    
    # reference layer for (non)active layers
    alay = actll[-1] if alay is None or alay not in actll \
        else alay
    if revactlay:
        ilay = revactlay[0] if ilay is None or ilay in actll \
            else ilay
    else:
        ilay = None
    print("align_correlate: Active layers: {}".format(actll))
    print(("align_correlate: alay, ilay: {}, {}"
        ).format(alay, ilay))

    # if there are layers to fit, preprocess misll and fitll
    if mislayi == 'all':
        misll = list(range(numlay))
    else:
        misll = mislayi
    
    fitll = ()
    if misll:
        misll = list({lay % numlay for lay in misll})
        
        # if at least one manual alignment coefficient is not
        # given use fitting by fitll
        if xcoefs is None or ycoefs is None:
            if not rellayi:
                fitll = [lay
                    for lay in range(numlay) if lay not in misll]            
            else:
                fitll = list({lay % numlay for lay in rellayi})
                misll = [lay for lay in misll if lay not in fitll]
            fitll = sorted(fitll)

        misll = sorted(misll)
        
        print("align_correlate: Refer. layers: {}".format(fitll))
        print("align_correlate: Missi. layers: {}".format(misll))

    return actll, revactlay, misll, fitll, alay, ilay


def _corr_align(numlay, rchan, actll, fitll, misll, rlay, arr,
    verbose):
    """cross-correlation stage of procedure align_correlate
    """
    
    print("align_correlate: Alignment according to channel"
        " '{}'.".format(rchan))
    relx, rely = np.zeros(numlay), np.zeros(numlay)
    arr1, k = arr[rlay], 1        
    
    if arr1 is None:
        print("align_correlate: Reference layer in reference "
            "channel full of NaNs. Halt.")
        raise InvalidInputError

    arr1[np.isnan(arr1)] = np.nanmean(arr1)
    
    # for each layer in actll or fitll find out the most
    # matching offsets
    arlist = sorted([lay
        for lay in range(numlay) if lay in actll or lay in fitll])
    
    if verbose:
        print("align_correlate: Cross-correlation:")
    for i in arlist:
        if verbose:
            print(("\tLayer {:2} out of {} being processed..."
                ).format(k, len(arlist)))
        k += 1
        if misll and i in misll: continue

        # evaluate correlation and get indices of the best
        # correlation match
        arrtemp, arr1tem = arr[i], arr1
        
        if arrtemp is None:
            relx[i], rely[i] = 0, 0
        else:
            arrtemp[np.isnan(arrtemp)] = np.nanmean(arrtemp)
            arr1tem[np.isnan(arrtemp)] = np.nanmean(arr1tem)

            relx[i], rely[i] = _get_offinds_corr(arr1tem, arrtemp)

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
    
def _postprocess_align(spmdata, rchan, relx, rely, fac,
    revactlay, actll, alay, ilay):
    """postprocessing stage of procedure align_correlate
    """
    
    xn, yn = np.array(spmdata.xnums), np.array(spmdata.ynums)
    xo = np.array([lay.xoffind[rchan] \
        if rchan in lay.channels.keys() else 0
        for lay in spmdata.layers])
    yo = np.array([lay.yoffind[rchan] \
        if rchan in lay.channels.keys() else 0
        for lay in spmdata.layers])

    # offsets        
    relx = np.round(relx - fac * xn / 2)
    rely = np.round(rely - fac * yn / 2)

    # reference offsets with respect to the array origin
    xoffind = relx
    yoffind = rely
    
    # for nonactive layers use the original offsets
    if revactlay:
        print("align_correlate: Nonactive layer offsets may"
            "be affected.")
        # now we set xoffind and yoffind of all channels to
        # values of xoffind and yoffind of the rchan channel
        xoffind[revactlay] = fac * xo[revactlay] + xoffind[alay]
        xoffind[revactlay] -= fac * xo[ilay]
        yoffind[revactlay] = fac * yo[revactlay] + yoffind[alay]
        yoffind[revactlay] -= fac * yo[ilay]
        
    # all offsets with respect to the array origin
    xoffmin, yoffmin = xoffind.min(), yoffind.min()
    xoffind -= xoffmin
    yoffind -= yoffmin

    # xoffmin and yoffmin are used in demonstrating plot
    xoffmin += relx.min()
    yoffmin += rely.min()

    return xn, yn, xoffind, yoffind, xoffmin, yoffmin


def _upsample_align(numlay, layers, rchan, fac, chanlist, order,
    verbose):
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
            if verbose:
                print(("align_correlate: Channel '{}' being"
                    " upsampled.").format(chan))
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
    fac, xn, yn, arrdict, downsample, order, verbose):
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
        aux = np.empty((spmdata.numlay, xlim, ylim),
            dtype=float)    
        aux.fill(np.nan)

        # fill new array with valid data and appropriate
        # offsets
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
            if verbose:          
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
            # ULOZENE V spmdata.layers[i].xnum A V
            # spmdata.layers[i].ynum
            
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
        spmdata.set_chan(chan, aux.copy(), updatelayers=True)                    

def align_correlate(spmdata, *cl, rchan=None, rlay=0, alay=None,
    ilay=None, actll=None, fitll=None, misll=None, show=False,
    xcoefs=None, ycoefs=None, fac=1, downsample=True, order=1,
    verbose=False):
    """align arrays according to rchan channel, i.e. offsets are
    calculated only once for rchan and then used multiple times
    for all channels in cl, unknown values are represented by
    NANs, resulting array may have bigger dimensions than original
    arrays
    
    cl - list of channels to be aligned, if empty, then all
        channels are aligned
    rlay - index of a fixed reference layer with which all layers
        are compared
    alay - if actll is not None or not all layers, then alay
        specifies to which layer inactive layers should adhere,
        i.e. offind of ilay layer will be set to the offind of
        alay; if alay is None or lies outside actll, then it is
        set to be the last layer in actll
    ilay - analog of alay for inactive layers, if None or in 
        actll, then it is set to be the first layer not in actll
    rchan - reference channel; all channels in cl are aligned
        according to channel 'rchan'
    actll - list of 'active' layers, which are to be aligned;
        if None, then all layers are used
    fitll - list of layers which are used to determine by fitting
        offsets for layers in misll, if None or empty, then
        fitll is effectively put equal to all layers not lying in
        misll
    misll - list of 'missing' layers for which offsets should
        be calculated by inter-/extra-polation;
        if None or empty, then correlation calculation is
        performed for all layers and no inter-/extra-polation is
        used; if misll and fitll are not disjoint, then misll is
        put equal to misll minus fitll;
        if misll='all', then all layers are used for 
        inter-/extra-polation
    show - if True, then plot depicting fitll points and
        extrapolation trends is shown,
        control points - fitll points
        calculated points - misll points
        x fit, y fit - fit of x-offsets and y-offsets,
            respectively
        (manual) - shown if fit is set manually via xcoefs or
            ycoefs
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
    verbose - if False, then almost all text documenting the
        routine progress is suppressed
    
    BEWARE - there might be problems for interpolation order > 1;
        due to round off errors the spmdata.layers[i][chan] arrays
        may contain NaNs with which map_coordinates cannot
        work for order > 1; this issue takes place in 
        '_reallocate_align' routine in downsampling stage,
        specifically in command np.linspace in its third parameter
 
        
    Example 1: I have given SPMdata structure 'dta' with 30
        layers. I want to align channels CURR and PHAS. As a
        reference channel I want to choose FREQ, as a reference
        layer I choose layer no. 3. That is, I write: 
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, rlay=3)'.
        But, from the result I see that cross-correlation does not
        work well for layers no. 10, 11, 12, 15, 16. I would like
        to use linear fit for offsets of these layers. So I append
        'misll=[10, 11, 12, 15, 16]' to argument list. New result
        I obtain nevertheless does not fit offsets for these
        layers well. Reason for this may be that offsets for
        layers outside 'misll' vary a lot, so I would like to
        calculate a linear fit only from offsets pertaining to
        layers no. 5, 6, 7. Thus I append 'fitll=[5, 6, 7]' to
        argument list and also set 'show' to True to see how
        offsets are calculated. At the end I have:
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ, rlay=3,
            misll=[10, 11, 12, 15, 16], fitll=[5, 6, 7])'.
            
    Example 2: Suppose I have given SPMdata structure as in
        Example 1. But now I want to align according to layer
        no. 0, which is default choice. Now I would like to align
        only layers no. 10 through no. 29 and leave layers no. 0
        through no. 5 unaffected (i.e. layers no. 0 through no. 5
        of channel 'rchan'! unaffected, since routine
        'align_corr_chan' aligns channels according to channel
        'rchan', all channels in 'cl' different from 'rchan' may
        get different offsets even for layers no. 0 through
        no. 5). So I write:
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ,
            actll=range(6, 30), show=True, misll=range(10, 30))'.
        But, the resulting alignment is not good, so I read actual
        offsets used by routine from the plot shown and set fit
        coefficients for x-offsets by hand, so I write
        'dta.align_corr_chan(CURR, PHAS, rchan=FREQ,
            actll=range(6, 30), show=True,
            misll=range(10, 30), xcoefs=(5., 10.))'.        
    """

    chanlist = list(spmdata.channels if not cl else cl)
    chanlist = [chan for chan in chanlist
                    if chan in spmdata.channels]
    
    # since each layer contains attributes xnum and ynum, i.e.
    # dimensions of valid data area, which are identical for
    # all channels, let us exclude possibility, that
    # 'align_correlate' changes dimension of valid data area for 
    # different channels
    if not downsample and fac > 1 and \
        len(chanlist) < len(spmdata.channels):
        print("align_correlate: Downsampling on, but not all"
            " channels included.")
        return

    numlay, layers = spmdata.numlay, spmdata.layers

    # layer lists initialization
    try:
        actll, revactlay, misll, fitll, alay, ilay = \
            _preprocess_align(numlay, rchan, actll, misll, fitll,
            xcoefs, ycoefs, alay, ilay)
    except InvalidInputError:
        return

    # upsampling        
    arr, arrdict = _upsample_align(numlay, layers,
        rchan, fac, chanlist, order, verbose)
   
    # cross-correlation     
    try:    
        relx, rely = _corr_align(numlay, rchan,
            actll, fitll, misll, rlay, arr, verbose)
    except InvalidInputError:
        return

    # inter-/extra-polation
    relx, rely, *coefs = _fit_align(xcoefs, ycoefs,
        misll, fitll, relx, rely)

    # postprocessing
    xn, yn, xoffind, yoffind, xoffmin, yoffmin = \
        _postprocess_align(spmdata, rchan, relx, rely, fac,
        revactlay, actll, alay, ilay)

    # plotting - show trends and control points
    if show:
        xfit, yfit = _showfit_align(numlay, coefs, fac,
            xn, yn, xoffmin, yoffmin)
        _show_align(numlay, fitll, revactlay, misll,
            xoffind, yoffind, xfit, yfit, xcoefs, ycoefs, coefs)

    # reallocation of arrays according to offsets and downsampling
    _reallocate_align(spmdata, chanlist, xoffind, yoffind, fac,
        xn, yn, arrdict, downsample, order, verbose)

    print(("align_correlate: Channels {} aligned."
        ).format(chanlist))            

# function determining total offsets
def _get_shifts(layers, lay, chanf, chanb, xdir, ydir, same):
    """get total x-shifts and y-shifts
    """

    # evaluate correlation and get indices of the best
    # correlation match
    layf = layers[lay].channels[chanf]
    layb = layers[lay].channels[chanb]
    valf, valb = np.nanmean(layf), np.nanmean(layb) 
    layf[np.isnan(layf)] = valf
    layb[np.isnan(layb)] = valb
        
    relx, rely = _get_offinds_corr(layf, layb)
    
    # OPET DOST BORDEL S TIM, CO JE TEDA X-OVA A Y-OVA OSA...
    # aby to vychazelo, tak takhle, ale WTF???!!!
    xdir, ydir = ydir, xdir 
        
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

def align_forback(spmdata, *cl, rchan=None, same=True, rlay=0,
    xdir=True, ydir=True, order=0):
    """align all channels in 'cl' according to 'rchan' channel
    so that forward and backward scans are moved towards each
    other and aligned 'to the centre', this routine does not
    align layers within one channel, but channels among
    themselves instead
    
    cl - sequence of channels to be aligned, if empty, than all
        channels are aligned
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

    # preprocess the list of channels to align, let cl contain
    # only forward scans
    claux = [chan
        for chan in spmdata.channels.keys()
        if not is_backward_scan(chan) and \
        get_backward_scan(chan) in spmdata.channels.keys()]
    if cl:
        claux = [chan for chan in claux if chan in cl]
    cl = claux
    print("align_forback: Channels to be aligned: {}.".format(cl))
    
    xshift, yshift = [None]*len(cl), [None]*len(cl)
    fxnum,  bxnum  = [None]*len(cl), [None]*len(cl)
    fynum,  bynum  = [None]*len(cl), [None]*len(cl)
    fxoffs, bxoffs = [None]*len(cl), [None]*len(cl)
    fyoffs, byoffs = [None]*len(cl), [None]*len(cl)

    # process reference channel and get corresponding
    # x-/y-offsets
    if rchan is None:
        pass
    elif rchan in spmdata.channels.keys() and \
        not is_backward_scan(rchan):
        rchanb = get_backward_scan(rchan)
        if rchanb not in spmdata.channels.keys():
            print("align_forback: Reference channel has no"
                " backward counterpart.")
            return
        else:                
            rxshift, ryshift = _get_shifts(spmdata.layers, rlay,
                rchan, rchanb, xdir, ydir, same)
    else:
        print("align_forback: Invalid reference channel.")
        return

    # find alignments for each pair of forw./backw. channels
    for i, chanf in enumerate(cl):
        chanb = get_backward_scan(chanf)
        
        # if rchan is None, then calculate offsets for each
        # channel separately...
        if rchan is None:
            xshift[i], yshift[i] = _get_shifts(spmdata.layers,
                rlay, chanf, chanb, xdir, ydir, same)            
        # ...otherwise use the same offsets as that for reference
        # channel
        else:
            xshift[i], yshift[i] = rxshift, ryshift

        # if fxnum != bxnum or fynum != bynum, then
        # resulting 3D arrays are not ensured to be
        # properly aligned w.r.t. forward/backward!
        fxnum[i], fynum[i] = spmdata.channels[chanf].shape[1:]
        bxnum[i], bynum[i] = spmdata.channels[chanb].shape[1:]
        
        fxoffs[i] = spmdata.layers[rlay].xoffind[chanf] 
        bxoffs[i] = spmdata.layers[rlay].xoffind[chanb] 
        
        fyoffs[i] = spmdata.layers[rlay].yoffind[chanf] 
        byoffs[i] = spmdata.layers[rlay].yoffind[chanb]
    
    # get dimensions of new arrays
    xlim = [int(abs(vi)) + max(lfi, lbi)
        for vi, lfi, lbi in zip(xshift, fxnum, bxnum)]
    ylim = [int(abs(vi)) + max(lfi, lbi)
        for vi, lfi, lbi in zip(yshift, fynum, bynum)]
     
    # if same == True, then modify fxshift, bxshift, ... so that
    # all channels have the same array dimension
    if same:
        # let all channels have the same dimension
        xlv, ylv = max(xlim), max(ylim)
        xlim, ylim = [xlv]*len(cl), [ylv]*len(cl)
        
        # dimensions of reference layer
        xnum = spmdata.layers[rlay].xnum
        ynum = spmdata.layers[rlay].ynum

        # calculate offsets for forward and backward channels

#        fxshift = [int((xlv - xnum - v)/2) - off
#            for v, f, off in zip(xshift, fxnum, fxoffs)] 
#        bxshift = [int((xlv - xnum + v)/2) - off
#            for v, b, off in zip(xshift, bxnum, bxoffs)] 
#        fyshift = [int((ylv - ynum - v)/2) - off
#            for v, f, off in zip(yshift, fynum, fyoffs)] 
#        byshift = [int((ylv - ynum + v)/2) - off
#            for v, b, off in zip(yshift, bynum, byoffs)] 

        fxshift = [int((xlv - f - v)/2)
            for v, f in zip(xshift, fxnum)] 
        bxshift = [int((xlv - b + v)/2)
            for v, b in zip(xshift, bxnum)] 
        fyshift = [int((ylv - f - v)/2)
            for v, f in zip(yshift, fynum)] 
        byshift = [int((ylv - b + v)/2)
            for v, b in zip(yshift, bynum)]
    else:
        # calculate offsets for forward and backward channels
        fxshift = [0 if v > 0 else int(-v) for v in xshift]
        bxshift = [0 if v < 0 else int(+v) for v in xshift]
        fyshift = [0 if v > 0 else int(-v) for v in yshift]
        byshift = [0 if v < 0 else int(+v) for v in yshift]

#    print("xshift: ", xshift)        
#    print("yshift: ", yshift)

#    print("fxshift: ", fxshift)        
#    print("fyshift: ", fyshift)
        
    # update 3D arrays
    for i, chanf in enumerate(cl):
        chanb = get_backward_scan(chanf)
        
        # create auxiliary arrays
        auxf = np.empty((spmdata.numlay, xlim[i], ylim[i]))
        auxb = np.empty((spmdata.numlay, xlim[i], ylim[i]))
        auxf.fill(np.nan)
        auxb.fill(np.nan)
        
#        print("hax: ", fxshift[i], fxnum[i])
#        print("hay: ", fyshift[i], fynum[i])
        
        # forward channel - fill array with data
        xoff, yoff = fxshift[i], fyshift[i]
        auxf[:, xoff:xoff + fxnum[i], yoff:yoff + fynum[i]
            ] = spmdata.channels[chanf]

        # backward channel - fill array with data
        xoff, yoff = bxshift[i], byshift[i]
        auxb[:, xoff:xoff + bxnum[i], yoff:yoff + bynum[i]
            ] = spmdata.channels[chanb]

        # update channels
        spmdata.set_chan(chanf, auxf, updatelayers=False)
        spmdata.set_chan(chanb, auxb, updatelayers=False)
        
        # update layer pointers and xoffinds, yoffinds
        for j, lay in enumerate(spmdata.layers):            
            xnum, ynum = lay.xnum, lay.ynum
            
            # forward channel
            if chanf in lay.channels.keys():
                xoff = lay.xoffind[chanf] + fxshift[i]
                yoff = lay.yoffind[chanf] + fyshift[i]
                lay.xoffind[chanf], lay.yoffind[chanf] = xoff, yoff                
                lay.channels[chanf] = spmdata.channels[chanf][j,
                    xoff:xoff + xnum, yoff:yoff + ynum]
#                print("for: xoff, yoff: ", xoff, yoff)

            # backward channel
            if chanb in lay.channels.keys():
                xoff = lay.xoffind[chanb] + bxshift[i]
                yoff = lay.yoffind[chanb] + byshift[i]
                lay.xoffind[chanb], lay.yoffind[chanb] = xoff, yoff
                lay.channels[chanb] = spmdata.channels[chanb][j,
                    xoff:xoff + xnum, yoff:yoff + ynum]
#                print("bck: xoff, yoff: ", xoff, yoff)

def align_offsets(self, *cl):
    """align arrays corresponding to each channel according their
    offsets, unknown values are represented by NANs, resulting
    array may have bigger dimensions than original arrays, it is
    assumed that xranges and yranges of all layers are identical
    
    cl - list of channels to be aligned; if None, then all
        channels are aligned
    """

    # POZOR!!! obracene x a y ve dvou prikazech nize, obr. takto
    # konecne vypada dobre, ale je zde nekonzistence
    # v tom, jake parametry nahravam do x a y slozek ze vstupniho
    # souboru nanonisem, patrne
    # PROBLEM JE V TOM, ZE PRVNI INDEX POLE SE VYKRESLUJE DO RADKU
    # OBRAZKU, ZATIMCO DRUHY INDEX POLE SE VYKRESLUJE
    # DO SLOUPCU OBRAZKU, TJ. INDEX x V NPARRAY JE CHAPAN JAKO
    # OSA y V OBRAZKU
    # KTEROU Z TECHTO MOZNOSTI POUZIVA NANONIS ??????

    # process list of channels to align
    cl = self.channels.keys() if not cl else cl
    cl = [chan for chan in self.channels.keys() if chan in cl]

#    xoffsets, yoffsets = self.xoffsets, self.yoffsets
    xn, yn = self.xnums, self.ynums
    
    for chan in cl:
        yoffsets, xoffsets = self.xoffsets, self.yoffsets
        xdimnew, ydimnew = max(xn), max(yn)
            
        # it is assumed that xranges and yranges of all layers are
        # identical, xstep and ystep are thus calculated according
        # to ranges of the first layer
    #    xstep = self.xranges[0] / xdimnew
    #    ystep = self.yranges[0] / ydimnew
    #    xsteps = [ran / xdimnew for ran in self.xranges]
    #    ysteps = [ran / ydimnew for ran in self.yranges]
        xstep, ystep = self.xsteps[chan], self.ysteps[chan]

        xdimnew += math.ceil((max(xoffsets) - min(xoffsets)) / xstep)
        ydimnew += math.ceil((max(yoffsets) - min(yoffsets)) / ystep)

    #    xdimnew += max([math.ceil((xoffset - min(xoffsets)) / xstep)
    #        for xoffset, xstep in zip(xoffsets, xsteps)])
    #    ydimnew += max([math.ceil((yoffset - min(yoffsets)) / ystep)
    #        for yoffset, ystep in zip(yoffsets, ysteps)])


        # xoffind = [math.floor((offset - min(xoffsets)) / xstep) for offset in xoffsets]
        # yoffind = [math.floor((offset - min(yoffsets)) / ystep) for offset in yoffsets]

        xoffind = [math.ceil((offset - min(xoffsets)) / xstep)
            for offset in xoffsets]
        yoffind = [math.ceil((offset - min(yoffsets)) / ystep)
            for offset in yoffsets]

#        xoffind = [math.ceil((offset - min(xoffsets)) / xstep)
#            for offset, xstep in zip(xoffsets, xsteps)]
#        yoffind = [math.ceil((offset - min(yoffsets)) / ystep)
#            for offset, ystep in zip(yoffsets, ysteps)]

##        xdimnew, ydimnew = ydimnew, xdimnew
##        xoffind, yoffind = yoffind, xoffind

        # create new array
        aux = np.empty((self.numlay, xdimnew, ydimnew),
            dtype=float)
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

        self.set_chan(chan, aux.copy(), updatelayers=True)

