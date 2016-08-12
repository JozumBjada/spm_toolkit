from mpl_toolkits.axes_grid1 import make_axes_locatable    
from scipy.ndimage.interpolation import \
    shift, zoom, rotate, affine_transform
import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from numpy.linalg import norm
import copy

# MASK_VAL is a padding value for masks, MASK_LIM is a value such
# that all mask values below MASK_LIM are interpreted as invalid
MASK_VAL = 10
MASK_LIM = 1

def show_array(array, numlay, laylist=None, limit=8, rowlimit=4,
    aspect='equal', interpolation='none', origin='lower'):
    """plot data block for various layers
    
    laylist - which layers to plot, sequence containing indices,
        default is first 'limit' layers
    limit - maximum number of layers plotted
    rowlimit - maximum number of plots per row
    aspect - aspect of imshow
    interpolation - how to interpolate values in imshow
    origin - either 'lower' or 'upper'; determines where imshow
        should put origin of the image
    """

    # only first 'limit' (or less) entries plotted, layloc
    # stores layer indices
    if laylist:
        layloc = [i for i in set(laylist) if i < numlay]
        layloc = sorted(layloc)[:limit]
    if not laylist or not layloc:
        layloc = range(min(limit, numlay))
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
        num="Data block")
        
    # create plots
    l = 0
    for k in layloc:
        ind = tuple(divmod(l, rowlimit))

    	# create image
        img = ax[ind].imshow(array[k],
            aspect=aspect, interpolation=interpolation,
            origin=origin)
        ax[ind].set_title("layer no. " + str(k % numlay))
        
        # add colorbar
        divider = make_axes_locatable(ax[ind])
        cax = divider.append_axes("right", size=0.1,
            pad=(0.5 if aspect == 'auto' else 0.1))
        plt.colorbar(img, cax=cax)
        
        # increment counter
        l += 1

    # improve layout and launch the interactive window
    fig.tight_layout()	
    plt.show()

def show_match(arr1, arr2, rlay=0, threshold=0.5, alpha=0.5,
    fig=None, ax=None):
    """show both blocks of data to compare how the two are
    similar to each other
    
    rlay - which layer to show
    threshold - the threshold for values in arrays, under
        which the given value is "truncated" to zero, and
        above which the value is set to one; if None, then no
        threshold is used
    alpha - the amount of opacity of the overlap
    fig, ax - if at least one of them is None, then new figure 
        and axes are created; otherwise image is drawn into fig
        and ax
    """
    
    TOLERANCY = 1e-5
    
    # prepare arr1
    arr1 -= np.nanmin(arr1)
    if abs(np.nanmax(arr1)) > TOLERANCY:
        arr1 /= np.nanmax(arr1)      
    arr1[np.isnan(arr1)] = 0.0
    if threshold is not None:
        arr1[arr1 < threshold] = 0.0
        arr1[arr1 > threshold / 2] = 1.0

    # prepare arr2
    arr2 -= np.nanmin(arr2)
    if abs(np.nanmax(arr2)) > TOLERANCY:
        arr2 /= np.nanmax(arr2) 
    arr2[np.isnan(arr2)] = 0.0
    if threshold is not None:    
        arr2[arr2 < threshold] = 0.0
        arr2[arr2 > threshold / 2] = 1.0
        
    # make sure both arrays have the same dimensions
    shape1, shape2 = arr1.shape, arr2.shape
    xsh = min(shape1[0], shape2[0])
    ysh = min(shape1[1], shape2[1])
    arr1 = arr1[:xsh, :ysh]
    arr2 = arr2[:xsh, :ysh]
        
    # create figure
    if fig is None or ax is None:
        fig, ax = plt.subplots(1, num="Overlap")
        fig.show()

    # if 'threshold' is a valid threshold value, then show the 
    # overlap of the two "thresheld" arrays; otherwise, when 
    # 'threshold' is None, the two arrays are used as two color
    # channels in a RGB image
    if threshold is not None:
        ax.imshow(arr2, cmap=plt.cm.Blues, origin="lower")
        ax.imshow(arr1, cmap=plt.cm.Reds, alpha=alpha,
            origin="lower")
    else:
        arr1, arr2 = 255 * arr1, 255 * arr2

        rgbarray = np.zeros((xsh, ysh, 3), dtype='uint8')
        rgbarray[:, :, 0] = arr1.astype('uint8')
        rgbarray[:, :, 1] = arr2.astype('uint8')
        ax.imshow(rgbarray)

    fig.canvas.show()

def corr_overlap(reg1, reg2):
    """function calculating correlation coefficient of reg1 and
    reg2 registration blocks;
    even though this routine works for arrays with different
    dimensions, it is assumed that prior to call of this routine
    one calls make_identical_dims method of CompareRegBlocks class 
    
    returns:
    res - if masks of blocks have no overlap, res = np.nan is 
        returned; otherwise res is equal to correlation coeff.
    """

    # retrieve data
    arr1 = reg1.arrtemp.flat
    arr2 = reg2.arrtemp.flat

    # retrieve masks
    mask1 = reg1.arrmask.flat
    mask2 = reg2.arrmask.flat
    
    # np.corrcoef needs both input arrays to have the same length;
    # moreover, arrtemp.flat and arrmask.flat are NOT lists, so we 
    # have to convert them to lists; this conversion is done 
    # automatically when one takes a slice of these objects
    length = min(len(arr1), len(arr2))
    arr1 = arr1[:length]
    arr2 = arr2[:length]

    mask1 = mask1[:length]
    mask2 = mask2[:length]

    # the area where both data from reg1 and reg2 are valid
    mask = mask1 * mask2

    # transform mask into boolean array
    mask = (mask > MASK_LIM)
    
    # if no valid data in the intersection, return NaN, otherwise
    # calculate the correlation coefficient...
    if not mask.any():
        res = np.nan
    else:
        # get valid data
        arr1, arr2 = arr1[mask], arr2[mask]

        # get correlation coefficient
        res = np.corrcoef(arr1, arr2)
        res = res[0, 1]
    
    return res

def find_parameters(shape1, shape2, controlpoints):
    """find the best match for control points; transformations are
    done in x-y plane only! more general method should be
    written...
    
    shape1, shape2 - shapes of blocks of data; this parameter is
        needed in order to find the appropriate offset of points
        during rotation and shift
    controlpoints - sequence of the form
        controlpoints = [points1, points2], where
        points1 = [[zcoord1, xcoord1, ycoord1], ...,
            [zcoordN, xcoordN, ycoordN]], where
        N >= 2 and zcoord1 is z-coordinate of the first control
        point, zcoord2 is z-coord. of the second control point,
        etc. and analogously for points2
        
        controlpoints are UPDATED BY THIS ROUTINE
        
    returns:
    parvalues - list of best-match values for all relevant
        parameters
    """

    points1, points2 = copy.deepcopy(controlpoints)
    
    # first line segments of both sets of control data
    segpoi1 = points1[1] - points1[0]
    segpoi2 = points2[1] - points2[0]

    # check whether line segments are not of zero length
    TOL = 1e-8
    if norm(segpoi1) < TOL or norm(segpoi2) < TOL:
        print("find_parameters: Two identical points."
            " Halt.")
        return

    # --- AUXILIARY ROUTINES ---

    def rotate_xy_point(pnt, o1, o2, phi):
        """rotate point 'pnt', i.e. 3-element sequence, by
        angle phi in x-y plane
        
        pnt - point to be rotated
        o1, o2 - x, y-shape of the old array
        phi - angle by which to rotate in degrees
        """

        # from degrees to radians
        phi *= 2 * np.pi / 360
        phi %= (2 * np.pi)
        cosphi, sinphi = np.cos(phi), np.sin(phi)
        
        o1, o2 = o1 / 2, o2 / 2
        n1 = o1 * cosphi + o2 * sinphi
        n2 = o2 * cosphi + o1 * sinphi

#        n1, n2 = n2, n1 # ???

        p1, p2 = pnt[1], pnt[2]
        
        pnt[1] = cosphi * (p1 - o1) - sinphi * (p2 - o2) + n1
        pnt[2] = sinphi * (p1 - o1) + cosphi * (p2 - o2) + n2
        
    def rotate_xy_frame(o1, o2, phi):
        """
        o1, o2 - x, y-shape of the old array
        phi - angle by which to rotate in degrees
        """

        # from degrees to radians
        phi *= 2 * np.pi / 360
        phi %= (2 * np.pi)
        cosphi, sinphi = np.cos(phi), np.sin(phi)

        n1 = o1 * cosphi + o2 * sinphi
        n2 = o2 * cosphi + o1 * sinphi
        
        return n1, n2

    def shift_xy_point(point, xoff, yoff):
        """shift point by xoff in x-direction and by yoff in
        y-direction
        """
        point += np.array([0, xoff, yoff])

    def zoom_xy_point(point, zoomfac):
        """zoom point in x-y plane by zoomfac
        """
        # zeroth component, i.e. z-coord., is unchanged
        point[1] *= zoomfac
        point[2] *= zoomfac

    def find_angle(vec1, vec2):
        """find angle between vectors vec1 and vec2
        """
        
        cosphi = np.dot(vec1, vec2)
        cosphi /= (norm(vec1) * norm(vec2))
        cosphi = np.clip(cosphi, -1.0, 1.0)

        # scalar product does not care about sense of
        # rotation, so we have to cope with that manually
        phi = np.arccos(cosphi)
        if (vec1 - vec2)[2] < 0:
            phi *= -1

        # phi in degrees
        phi *= 360 / (2 * np.pi)

        return phi

    # --- ZOOM ---

    # find zoom factor
    zoomfac = norm(segpoi2) / norm(segpoi1)

    # update control points
    for point in points1:
        zoom_xy_point(point, zoomfac)
    
    # prepare auxiliary parameters for rotations and shifts
    xsh1, ysh1 = shape1[1:]
    xsh1, ysh1 = xsh1 * zoomfac, ysh1 * zoomfac
    
    xsh2, ysh2 = shape2[1:]
    
    # --- ROTATION ---
    
    # phi1 and phi2 are angles between 'vertical' and first
    # segments of 'points1' and 'points2', respectively;
    # assumes only rotation in x-y plane, hence the 0-th
    # component of vectors is unaffected

    # 'vertical' must be a unit vector!
    vertical = np.array([0, 1, 0])
    vertical /= norm(vertical)

    # find dphi, i.e. angle of rotation
    phi1 = find_angle(segpoi1, vertical)
    phi2 = find_angle(segpoi2, vertical)        
    dphi = phi2 - phi1

    xsh1, ysh1 = ysh1, xsh1 # WTF?!
    xsh2, ysh2 = ysh2, xsh2 # WTF?!

    # update control points
    for point in points1:
        rotate_xy_point(point, xsh1, ysh1, dphi)

    # --- SHIFT ---
    
    # find shift
    seqshift = points2[0] - points1[0]
    zoff, xoff, yoff = seqshift.copy()

    # update control points
    for point in points1:
        shift_xy_point(point, xoff, yoff)

#    xoff, yoff = yoff, xoff # ???

    # we want xoff and yoff to store only how the origins are
    # positioned with respect to each other; so we need to get rid
    # of offsets that are caused by enlarging the array due to
    # rotation
    xsh1new, ysh1new = rotate_xy_frame(xsh1, ysh1, dphi)    
    xoff = (xsh2 - xsh1new) / 2 - xoff
    yoff = (ysh2 - ysh1new) / 2 - yoff
    
#    xoff, yoff = yoff, -xoff #

    # --- FINALIZE ---

    # store resulting values of parameters
    parvalues = [zoomfac - 1, dphi, xoff, yoff]

    # save final form of control points;
    # we want to preserve pointer to original controlpoints and 
    # update it by new data, so we update each item separately
    controlpoints[0] = copy.deepcopy(points1)
    controlpoints[1] = copy.deepcopy(points2)

    return parvalues
