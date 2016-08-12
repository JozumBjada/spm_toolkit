from .registration_base import *
from .spm_data_inspect import inspect_channels
import copy

# BACHA: PORADI SLOZEK V arrorig A mask JE: z, x, y !!!

"""module containing definition of class RegBlock, which is a 
wrapper for data used for registration
"""

class RegBlock:
    """class containing data used for registrations
    
    internally, there are three np.ndarrays for given data;
    'arrorig' corresponds to the original data
    'arrtemp' corresponds to transformed data
    'arrmask' is a mask of arrtemp
    when the result of a given transformation is approved by a
    user, i.e. when the current form of arrtemp is approved, then 
    arrorig may be updated to contain data from arrtemp; at that
    point, new transformation may be applied to arrorig to obtain
    new arrtemp; if the new form of arrtemp is not good, it may be
    discarded and overriden by calling any other transformation
    on arrorig; when the user is satisfied with the outcome, (s)he
    may update arrorig by a current form of arrtemp; and so on and
    so forth... 
    
    """
    
    def __init__(self, inp=None, description="No description.",
        order=0):
        """initialize the registration block
        
        inp - if inp is string, then this string is expected to be
            a valid xsf-file name and data are loaded from this
            file;
            if inp is np.ndarray, then this array is used directly
            as the data
        description - sometimes it may be handy to mention where
            the data came from or what they represent; such
            metadata may be stored as a string of description
        order - order of interpolations done when transforming
            the data
        """
    
        self.order = order
    
        # initialize values of parameters used for reconstruction
        self.parsorig = [0]*4
        self.parstemp = [0]*4
        
        if isinstance(inp, str):
            self.load_xsf(inp)            
            self.filename = inp
            self.description = description
        
        elif isinstance(inp, np.ndarray):
        
            # original data
            self.arrorig = inp.copy()
            
            # temporary data used for comparison with other 
            # registration block
            self.arrtemp = inp.copy()
            
            # when one pads array with NaNs, then problems with
            # numpy routines occurs, since they cannot work with
            # NaNs; moreover, when one pads array with zeros, then
            # unwanted artifacts emerge near the border of valid
            # data; that is why the average value is used for 
            # padding
            paddingvalue = np.nanmean(self.arrtemp)
            self.arrtemp[np.isnan(inp)] = paddingvalue
            
            # mask of temporary data
            # instead of boolean values integer values are used
            # in order to prevent artifacts and unwanted behaviour
            # of results coming from transformation routines;
            # the padding value MASK_VAL should be large enough,
            # in calculation all data values, for which
            # corresponding mask value is above MASK_LIM are
            # taken to be valid
            self.arrmask = np.empty_like(inp, dtype=int)
            self.arrmask.fill(MASK_VAL)
            self.arrmask[np.isnan(inp)] = 0
            
            self.filename = ""
            self.description = description
            
        else:
            self.arrorig = None
            self.arrtemp = None
            self.arrmask = None
            self.filename = ""
            self.description = description

        # reference points used to mimic transformations of
        # arrays;
        # these points are used to mimic transformations only in
        # x-y plane; for more general transforms one should review
        # and generalize the code;
        # points form an L-like shape (z, x, y) situated
        # approximately in the centre of the arrays, with lines 
        # being VAL-pixels long
        VAL = 100

        if inp is None:
            o1, o2 = 0, 0            
        else:
            o1, o2 = self.arrorig.shape[1:]

        # due to imaging routines, swap the values
        o1, o2 = o2, o1

        pt1 = np.array([0, o1 / 2 + VAL, o2 / 2])
        pt2 = np.array([0, o1 / 2, o2 / 2])
        pt3 = np.array([0, o1 / 2, o2 / 2 + VAL])

        self.set_points([pt1, pt2, pt3], np.array([0, 0, 0]))
        
    def set_points(self, inppts, inporigin):
        """replace self.ptsorig and self.ptstemp by inppts and
        similarly for origin point
        
        inppts - sequence of three points; inppts is of the form:
            [np.array([z1, x1, y1]), ..., np.array([z3, x3, y3])]
        """

        # initial form of points attached to the initial value
        # of arrorig; unlike ptsorig, these points are usually
        # NOT changed in calculation; they serve as a reference
        self.ptsorig = copy.deepcopy(inppts)

        # points attached to arrtemp
        self.ptstemp = copy.deepcopy(inppts)

        # original form of point of origin
        self.oriorig = inporigin.copy()
        
        # temporary form of point of origin
        self.oritemp = inporigin.copy()
        
    # DATA AND METADATA RETRIEVAL
   
    @property
    def xy_pos(self):
        """return x- and y-position of the origin; these two
        values are used to correctly align one registration block
        with another in CompareRegBlocks routine
        'make_identical_dims'
        """
        return np.array(self.oritemp[-2:])
        
    @property
    def numlay(self):
        """number of layers
        """
        
        return self.arrtemp.shape[0]

    def print_pars(self):
        print("print_pars: parsorig = ", self.parsorig)
        print("            parstemp = ", self.parstemp)

    def print_shapes(self):
        print("print_shape: arrorig shape = ", self.arrorig.shape)
        print("             arrtemp shape = ", self.arrtemp.shape)
        print("             arrmask shape = ", self.arrmask.shape)
         
    def print_description(self):
        print("print_description: ", self.description)

    def print_points(self):
        print(("print_points:\n\tptsorig = {}\n\tptstemp = {}"
            ).format(self.ptsorig, self.ptstemp))

    def print_origin(self):
        print(("print_origin:\n\toriorig = {}\n"
            "\toritemp = {}").format(self.oriorig, self.oritemp))

    def value_range(self):
        """return the minimal and maximal values appearing in
        arrtemp
        """
        
        return np.nanmin(self.arrtemp), np.nanmax(self.arrtemp)

    @property
    def arr_with_nans(self):
        """return arr, which is equal to arrtemp, but zeros are
        replaced by NaNs
        """
        
        arr = self.arrtemp.copy()
        arr[self.arrmask < MASK_LIM] = np.nan
        
        return arr

    # GRAPHICAL OUTPUT
    
    def show(self, data='temp', **kwargs):
        """show data
        
        data - which data to show; data can assume one of values:
            temp - arrtemp is shown
            orig - arrorig is shown
            mask - arrmask is shown
        kwargs - all keyword arguments allowed by routine
            show_array
        """
    
        if data == 'temp':
            arr = self.arrtemp
        elif data == 'orig':
            arr = self.arrorig
        elif data == 'mask':
            arr = self.arrmask
        else:
            print("show: Invalid data specification.")
            return
    
        show_array(arr, arr.shape[0], **kwargs)

    def inspect_block(self, refpts=True):
        """inspect_channels routine modified for use with data
        from registration blocks; two blocks are shown,
        one correspond to arrtemp, the other to the mask        
        
        refpts - if True, then reference points are shown as well
        """
        
        # data to be shown
        data1 = self.arr_with_nans
        data2 = self.arrmask
        
        # names of blocks as depicted by inspect_channels
        NAME1, NAME2 = "data", "mask"
        
        # inspect_channels requires SPMdata-like structures as
        # inputs, so we create auxiliary class with atributes 
        # necessary for inspect_channels, which have no effective
        # values
        class AuxSPM:
            def __init__(self):
                self.channels = {NAME1:data1, NAME2:data2}
                self.xsteps = {NAME1:1, NAME2:1}
                self.ysteps = {NAME1:1, NAME2:1}
                self.units  = {NAME1:"", NAME2:""}
                self.numlay = data1.shape[0]
        
        # create SPMdata-like object
        aux1 = AuxSPM()
        
        # if refpts, include reference points
        if refpts:
            ctrlpts = [self.ptstemp.copy(), self.ptstemp.copy()]
        else:
            ctrlpts = []
        
        # show blocks
        inspect_channels([aux1, aux1],[NAME1, NAME2],
            showcut=True, cutlist=ctrlpts, transposectrlpts=False)

    # DATA MANIPULATION
        
    def strip_nans(self):
        """reshape arrays so that arrorig is stripped from NaNs
        """
        
        self.update_orig()
        
        # shrink along the z-axis
        m = ~(self.arrmask < MASK_LIM).all(axis=0)

        # find NaNs in rows and columns
        mrow, mclm = m.any(axis=1), m.any(axis=0)
        minrow, minclm = np.argmax(mrow), np.argmax(mclm)
        maxrow = len(mrow) - np.argmax(mrow[::-1])
        maxclm = len(mclm) - np.argmax(mclm[::-1])

        # update arrays
        self.arrorig = self.arrorig[:,
            minrow:maxrow, minclm:maxclm]        

        # update reference points
        for point in self.ptsorig:
            point -= [0, minclm, len(mrow) - maxrow]
    
        # update everything else
        self.reset_temp()
    
    def delete_layers(self, laylist, leave=False):
        """delete layers specified by 'laylist'
        
        leave - if False, then laylist is understood as a sequence
            of layers to be deleted;
            if True, then laylist is understood as a sequence of
            layers that are left in the structure, the rest is
            deleted
        """
        
        # delete or leave?
        if leave:
            laylist = [i for i in range(self.numlay)
                if i not in laylist]
        
        # list of layers to be preserved
        complist = [i for i in range(self.numlay)
            if i not in laylist]
                       
        self.arrorig = self.arrorig[complist]
        self.arrtemp = self.arrtemp[complist]
        self.arrmask = self.arrmask[complist]

    def clip(self, lim):
        """clip values of arrtemp according to limits in lim
        
        lim - sequence of the form (min, max)
        """
        
        self.arrtemp = np.clip(self.arrtemp, lim[0], lim[1])
        
    def tile(self, reps):
        """repeat data in arrays by number of times given
        by 'reps'
        """
        
        self.arrorig = np.tile(self.arrorig, reps)
        self.arrtemp = np.tile(self.arrtemp, reps)
        self.arrmask = np.tile(self.arrmask, reps)
        
        # TODO: NUTNO PŘIDAT TRANSFORMACI self.ptstemp, self.oritemp!, JE ZDE NUTNÁ?
            
    def pad_like(self, arr1):
        """pad arrtemp and arrmask so that x-y dimension is equal
        to that of 'arr1'; if one dimension is greater for self
        than for 'arr1', then this dimension is left untouched;
        padding is done symmetrically on both sides of a given 
        dimension
        
        return:
            xoff, yoff - offsets of padding
        """
        
        arr2 = self.arrtemp
        zsh, xsh, ysh = arr2.shape
        
        # get offsets
        # no floor function yet!!!
        xoff = (arr1.shape[1] - xsh) / 2
        yoff = (arr1.shape[2] - ysh) / 2

        xoff = 0 if xoff < 0 else xoff        
        yoff = 0 if yoff < 0 else yoff
        
        newshape = [zsh,
            np.ceil(xsh + 2 * xoff), np.ceil(ysh + 2 * yoff)]
        
        # "zero-checking"
        if xoff == 0 and yoff == 0:
            return 0, 0
        else:
            # floor function only now, after "zero-checking"
            # performed
            xoff = np.floor(xoff)
            yoff = np.floor(yoff)
          
        # update arrtemp
        paddingvalue = np.nanmean(self.arrtemp)
          
        arr2 = np.empty(newshape)
        arr2.fill(paddingvalue)
        arr2[:, xoff:xoff + xsh, yoff:yoff + ysh] = self.arrtemp
        self.arrtemp = arr2

        # update mask
        arr2 = np.empty(newshape)
        arr2.fill(0)
        arr2[:, xoff:xoff + xsh, yoff:yoff + ysh] = self.arrmask
        self.arrmask = arr2
        
#        print("xoff, yoff: ", xoff, yoff)
        
        # update points
        for pt in self.ptstemp:
            pt += np.array([0, xoff, yoff])
        
        return xoff, yoff

    def cut_data(self, indices):
        """cut data from the 3D arrays and leave only these;
        i.e. for all arrays it is done effectively:
        self.arrays = self.arrays[:, ixl:ixu, iyl:iyu], where
        indices = [ixl, ixu, iyl, iyu]
        
        indices - list of bounds for slicing of the form above
        """
        
        # take indices
        ixl, ixu, iyl, iyu = indices
        
        self.arrorig = self.arrorig[:, ixl:ixu, iyl:iyu]
        self.arrtemp = self.arrtemp[:, ixl:ixu, iyl:iyu]
        self.arrmask = self.arrmask[:, ixl:ixu, iyl:iyu]
            
    def sort(self, keylist):
        """sort individual layers of data arrays by 'keylist'
        
        keylist - sequence of integers corresponding to indices
            of newly sorted layers
        """
        
        keylist = list(keylist)
        
        self.arrorig = self.arrorig[keylist]
        self.arrtemp = self.arrtemp[keylist]
        self.arrmask = self.arrmask[keylist]
        
    # APPLY TRANSFORMATIONS, UPDATE AND RESET THESE APPLICATIONS

    def apply_transf(self, pars=None):
        """to arrtemp apply a series of transformations with 
        parameters taken from 'pars'; the order of these reads
        as follows:
            zoom_xy, rotate_xy, shift_x, shift_y
            
        pars - sequence of parameters of the form above; if None,
            then self.parstemp are used
        """
        
        if pars is None: pars = self.parstemp
        zoomfac, rotfac, xshfac, yshfac = pars
        
#        print("xshfac, yshfac: ", xshfac, yshfac)
        
        self.zoom_xy(zoomfac)
        self.rotate_xy(rotfac)
        self.shift_x(xshfac)
        self.shift_y(yshfac)
    
    def reset_temp(self):
        """from arrorig create brand new arrtemp and arrmask;
        moreover, from ptsorig create brand new ptstemp
        """
        
        # reset arrtemp
        self.arrtemp = self.arrorig.copy()
        
        # get rid of nans
        paddingvalue = np.nanmean(self.arrtemp)
        self.arrtemp[np.isnan(self.arrorig)] = paddingvalue

        # reset mask
        self.arrmask = np.empty_like(self.arrorig, dtype=int)
        self.arrmask.fill(MASK_VAL)
        self.arrmask[np.isnan(self.arrorig)] = 0
        
        # reset points
        self.ptstemp = copy.deepcopy(self.ptsorig)
        self.oritemp = copy.deepcopy(self.oriorig)

        # reset parameters
        self.parstemp = copy.deepcopy(self.parsorig)

    def update_orig(self, inp=None):
        """replaces 'arrorig' by 'inp' and updates 'arrtemp'
        and 'arrmask' accordingly; if inp=None, then 'arrtemp'
        is used as 'inp';
        moreover, ptsorig is replaced by ptstemp and similarly for
        origin point
        
        ASSUMES inp IS PADDED BY NaNs
        """
        
        if inp is None:
            inp = self.arrtemp
            inp[self.arrmask < MASK_LIM] = np.nan
        
        self.arrorig = inp.copy()
        self.ptsorig = copy.deepcopy(self.ptstemp)
        self.oriorig = copy.deepcopy(self.oritemp)
        self.parsorig = copy.deepcopy(self.parstemp)
        
        self.reset_temp()
    
    # INDIVIDUAL TRANSFORMATIONS
    
    def rotate_xy(self, angle, reshape=True):
        """rotate the block in x-y plane
        
        angle - rot. angle in degrees; or parameter of the form:
            angle = [o1, o2, n1, n2, phi], where phi is the
            rotation angle in degrees
        reshape - if True, then arrays are reshaped so as to
            cover the whole rotated data
        """

        # tolerancy
        TOL = 1e-6
    
        # angle = [o1, o2, n1, n2, phi], get phi
        if isinstance(angle, (list, tuple)):
            angle = angle[-1]

        # old shape
        o1, o2 = self.arrtemp.shape[1:]

        # if angle is zero, do nothing, otherwise rotate
        if np.abs(angle) < TOL:
            pass
        else:
            # arrtemp array is padded by paddingvalue to minimize
            # numeric artefacts
            paddingvalue = np.nanmean(self.arrtemp)
            self.arrtemp = rotate(self.arrtemp, angle,
                axes=(1, 2), reshape=reshape, order=self.order,
                cval=paddingvalue)
            self.arrmask = rotate(self.arrmask, angle,
                axes=(1, 2), reshape=reshape, order=self.order,
                cval=0)
            
        # new shape
        n1, n2 = self.arrtemp.shape[1:]

        # we have to swap values due to imaging issues, but don't
        # ask me why...WTF?!
        o1, o2 = o2, o1
        n1, n2 = n2, n1

        # rotate points in x-y plane
        # phi is equal to angle, but in radians        
        o1h, o2h = o1 / 2, o2 / 2
        phi = angle * (2 * np.pi) / 360
        cosphi, sinphi = np.cos(phi), np.sin(phi)
        for pt in self.ptstemp:
            p1, p2 = pt[1], pt[2]
            
            pt[1] = cosphi * (p1 - o1h) - sinphi * (p2 - o2h)
            pt[1] += n1 / 2
            
            pt[2] = sinphi * (p1 - o1h) + cosphi * (p2 - o2h)
            pt[2] += n2 / 2

        return [o1, o2, n1, n2, angle]

    def shift(self, seqshift=None):
        """shift arrays; real shift is not performed, but its
        value is stored and used by CompareRegBlocks routines
        'make_identical_dims' and 'reconstruct_arrays' to actually
        make arrays with appropriate shift
        """

        # why [seqshift[2], -seqshift[1]] instead of 
        # [seqshift[1], seqshift[2]]? WTF?
        self.oritemp += np.array(
            [seqshift[0], seqshift[2], -seqshift[1]])

    def shift_x(self, shift):
        self.shift([0, shift, 0])
        
    def shift_y(self, shift):
        self.shift([0, 0, shift])

    def zoom_xy(self, fac):
        """zoom arrays and points in x-y plane
        
        fac - fac is expected to be of the form:
            zoomfac = FACTOR**fac, i.e. if FACTOR = 2, then
            increasing the size of arrays eight times (i.e. 
            zoomfac = 8) means that fac = 3;
            values of this form are easier to work with
        """

        # tolerancy and zoom factor
        TOL = 1e-6 
        FACTOR = 2
        zoomfac = FACTOR**fac
        
        # if fac is zero, i.e. zoomfac is one, do nothing   
        if abs(fac) < TOL: return

#        print("zoom for: ", fac)

        # update arrays
        # ndimage routine zoom() raises a warning that is not 
        # relevant in this context, so we ignore it
        import warnings        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            
            paddingvalue = np.nanmean(self.arrtemp)
            
            # if there is only a single layer...
            if self.arrtemp.shape[0] == 1:
                zoomseq = [zoomfac, zoomfac]

                self.arrtemp = zoom(self.arrtemp[0], zoomseq,
                    cval=paddingvalue)[np.newaxis]
                self.arrmask = zoom(self.arrmask[0], zoomseq,
                    cval=0)[np.newaxis]

            else: # i.e. self.arrtemp.shape[0] > 1
                # for many-layered arrays...
                zoomseq = [1, zoomfac, zoomfac]
                
                self.arrtemp = zoom(self.arrtemp, zoomseq,
                    cval=paddingvalue)
                self.arrmask = zoom(self.arrmask, zoomseq,
                    cval=0)
        
        # update points
        for pt in self.ptstemp:
            pt *= zoomfac
        
    def skew_x(self, skewpar):

        TOL = 1e-6    
        if abs(skewpar) < TOL: return
    
#        transmat = np.array([[1, 0, 0], [skewpar, 1, 0], [0, 0, 1]])
        transmat = np.array([[1, 0, 0], [0, 1, skewpar], [0, 0, 1]])
        paddingvalue = np.nanmean(self.arrtemp)
        self.arrtemp = affine_transform(self.arrtemp, transmat,
            order=self.order, cval=paddingvalue)
        self.arrmask = affine_transform(self.arrmask, transmat,
            order=self.order, cval=0)
            
        # TODO: NUTNO PŘIDAT TRANSFORMACI self.ptstemp, self.oritemp!

    def skew_y(self, skewpar):

        TOL = 1e-6    
        if abs(skewpar) < TOL: return

        transmat = np.array([[1, 0, 0], [0, 1, 0], [skewpar, 0, 1]])
        paddingvalue = np.nanmean(self.arrtemp)
        self.arrtemp = affine_transform(self.arrtemp, transmat,
            order=self.order, cval=paddingvalue)
        self.arrmask = affine_transform(self.arrmask, transmat,
            order=self.order, cval=0)

        # TODO: NUTNO PŘIDAT TRANSFORMACI self.ptstemp, self.oritemp!

    # LOADING FROM EXTERNAL FILES (IMPORT AND EXPORT IS PARTIALLY
    # SUPPORTED BY CompareRegBlocks)

    def load_xsf(self, path):
        """load data from .xsf format file
        """
        try:        
            with open(path, "r") as f:
                # read the header
                line = f.readline().strip()
                while not line.startswith("BEGIN_DATAGRID_3D"):
                    line = f.readline().strip()
        
                # get array dimensions
                xdim, ydim, zdim = f.readline().strip().split()
                xdim, ydim, zdim = int(xdim), int(ydim), int(zdim)
                
                # skip the crystal grid specification
                for i in range(4):
                    line = f.readline()
                
                # read data
                NUMCHARLINE = 30 # approx. num. of char. per line
                lines = f.readlines(NUMCHARLINE*xdim*ydim*zdim)
                lines = list(map(float, lines[:-2]))
                
                # store data in arrays
                self.arrorig = np.array(lines)
                self.arrorig.shape = (zdim, ydim, xdim)
                self.arrorig = np.rollaxis(self.arrorig, 2, 1)
                
                self.arrtemp = self.arrorig.copy()
                paddingvalue = np.nanmean(self.arrorig)
                self.arrtemp[np.isnan(self.arrorig)] = paddingvalue
                
                self.arrmask = np.empty_like(
                    self.arrorig, dtype=int)
                self.arrmask.fill(MASK_VAL)
                self.arrmask[np.isnan(self.arrorig)] = 0
        except FileNotFoundError:
            print("load_xsf: No file to be found.")
            return

