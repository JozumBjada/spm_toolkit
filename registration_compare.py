import struct, os.path
from .spm_data_inspect import inspect_channels
from .registration_class import *

"""module containing definition of class CompareRegBlocks, which
takes two RegBlock objects and finds values of parameters of
several given transformations, application of which drives one
block to resemble the second block
"""

class CompareRegBlocks():
    """class responsible for finding the best settings of 
    parameters by which experimental and theoretical
    data are transformed to each other
    
    the instance of this class, say 'comp', is initialized by
    specifying two instances of RegBlock class; the second
    instance is kept (almost) unaffected while the first one is
    transformed appropriately by calling various methods of 'comp'
    
    example of a work with CompareRegBlocks may look as follows:
    >>> arr1 = ...np.ndarray of data...
    >>> arr2 = ...np.ndarray of data...
    >>> reg1, reg2 = RegBlock(arr1), RegBlock(arr2)
    >>> comp = CompareRegBlocks(reg1, reg2)
    >>> comp.find_match(...)
    >>> comp. ...any other method...
    >>> reg1.show()
    >>> reg2.show()
    at first, we create two data arrays and use them to create
    two registration blocks reg1 and reg2; then we create and 
    initialize CompareRegBlocks instance comp by these two blocks;
    after that we perform some adjustments and transformations on
    these two blocks via calling method find_match on comp; these
    transformations affect the form of reg1 and reg2, so the
    changes may be seen by calling method show on reg1 and reg2
    """

    def __init__(self, reg1=None, reg2=None,
        compfun=corr_overlap):
        """initialization
        
        reg1 - registration block to be transformed; if None, see
            remark below
        reg2 - registration block for comparison with reg1;
            if None, see remark below
        compfun - function used to calculate overlap,
            correlation between reg1 and reg2
            
        if reg1 is None or reg2 is None, then self.reg1 and 
        self.reg2 are NOT initialized and empty object is created,
        this can be filled by data imported from external file
        """
        
        if reg1 is None or reg2 is None:
            self.reg1 = RegBlock()
            self.reg2 = RegBlock()
        else:
            self.reg1 = reg1
            self.reg2 = reg2

            # make initial reference points identical for both
            # blocks of data
            self.reg2.set_points(
                self.reg1.ptsorig, self.reg1.oriorig)

        self.compfun = compfun

        self.fig = None
        self.ax = None
        self.overlaplist = []
        self.reconstructed = False

        # default values for dictpar for routine 'find_match'
        self.matchpars = {
            "zoom":[1, 0.2], "rotate":[1, 40],
            "xshift":[1, 4], "yshift":[1, 4],
            "xskew":[1, 0], "yskew":[1, 0]
            }

    # DATA MANIPULATION    

    def delete_layers(self, laylist, leave=False):
        """delete layers specified by 'laylist'
        
        leave - if False, then laylist is understood as a sequence
            of layers to be deleted;
            if True, then laylist is understood as a sequence of
            layers that are left in the structure, the rest is
            deleted
        """
        
        self.reg1.delete_layers(laylist, leave=leave)
        self.reg2.delete_layers(laylist, leave=leave)

        if self.reg1.numlay != self.reg2.numlay:
            print("delete_layers: Both blocks have unequal "
                "number of layers.")
        else:
            print(("delete_layers: There remain {} layer{}."
                ).format(self.reg1.numlay,
                "" if self.reg1.numlay == 1 else "s"))

    def strip_nans(self):
        """strip nans from both registration blocks SEPARATELY!,
        i.e. no mutual position of the blocks is taken into
        account
        """
        
        self.reg1.strip_nans()
        self.reg2.strip_nans()
        
    # "RECONSTRUCTION", "DE-RECONSTRUCTION"

    def reconstruct_arrays(self, pars=None, ctrlpts=None):
        """apply all necessary transformations onto reg1 and 
        suitably reallocate arrtmps in reg1 and reg2, so these
        two match together
        
        pars - parameter values to be applied to reg1; if None,
            then self.reg1.parstemp are used
        ctrlpts - control points
        """
        
        # if arrays are already reconstructed, do nothing
        if self.reconstructed: return
        
        # otherwise apply all necessary transformations to reg1
        # and reallocate arrays in both reg1 and reg2 so these 
        # have the same dimensions and correct positions of valid
        # data
        self.reg1.apply_transf(pars)
        self.make_identical_dims(ctrlpts=ctrlpts)
        
        # data are reconstructed, so set corresponding flag
        # accordingly
        self.reconstructed = True

#        self.reg1.print_points()        
#        self.reg1.print_origin()
#        self.reg1.print_pars()
        
    def dereconstruct_arrays(self):
        """the effective inverse to reconstruct_arrays; arrtmps
        are reset to have the same values as arrorigs and flag
        self.reconstructed is set to False if was True beforehand
        
        BEWARE: routine 'reconstruct_arrays' has optional argument 
            'ctrlpts' that gets updated during calling of this
            routine; nevertheless, no reset to original values is
            implemented by 'dereconstruct_arrays' to 'ctrlpts'
        """
  
        # if the arrays are not reconstructed, do nothing;
        # otherwise reset arrays to their original forms
        if self.reconstructed:
            self.reg1.reset_temp()
            self.reg2.reset_temp()
            self.reconstructed = False
            
    def make_identical_dims(self, ctrlpts=None):
        """pad reg1.arrtemp, reg2.arrtemp and their masks so they
        have identical dimensions and are with correct mutual 
        position
        
        ctrlpts - control points; if None, nothing is done
        """
        
        # shapes of arrays
        zsh1, xsh1, ysh1 = self.reg1.arrtemp.shape
        zsh2, xsh2, ysh2 = self.reg2.arrtemp.shape

        # how far is the origin of reg2 from origin of reg1
        dxor, dyor = self.reg2.xy_pos - self.reg1.xy_pos

        dxor, dyor = -dxor, -dyor # WTF?!

#        print("dxor: ", dxor)

        xsh1h, ysh1h = xsh1 / 2, ysh1 / 2
        xd2, yd2 = dxor + xsh2 / 2, dyor + ysh2 / 2

        # shape of new arrays
        xshape = max(xsh1h, xd2) + max(xsh1h, xsh2 - xd2)
        yshape = max(ysh1h, yd2) + max(ysh1h, ysh2 - yd2)
        
        # due to round-off errors we make sure that whole data
        # contained by newly allocated arrays by enlarging shape
        xshape += 1
        yshape += 1

        # if reg1 is more 'to the left' then reg2 in x-direction
        # set xoff1 equal to zero and calculate appropriate xoff2;
        # similarly for all other cases
        if np.argmax([xsh1h, xd2]) == 0:
            xoff1, xoff2 = 0, xsh1h - xd2
        else:
            xoff1, xoff2 = -(xsh1h - xd2), 0

        if np.argmax([ysh1h, yd2]) == 0:
            yoff1, yoff2 = 0, ysh1h - yd2
        else:
            yoff1, yoff2 = -(ysh1h - yd2), 0
        
        # xoff1 etc. are index offsets
        xoff1, yoff1 = np.round(xoff1), np.round(yoff1)
        xoff2, yoff2 = np.round(xoff2), np.round(yoff2)
            
        # allocate new arrays
        arr1  = np.empty((zsh1, xshape, yshape))
        arr2  = np.empty((zsh2, xshape, yshape))
        
        mask1 = np.empty((zsh1, xshape, yshape))
        mask2 = np.empty((zsh2, xshape, yshape))

        # fill arrays with appropriate padding values
        paddingvalue = self.reg1.arrtemp.mean()
        arr1.fill(paddingvalue)
        
        paddingvalue = self.reg2.arrtemp.mean()
        arr2.fill(paddingvalue)

        mask1.fill(0)
        mask2.fill(0)

        # fill arr1 with data in an appropriate position
        arr1[:, xoff1:xoff1 + xsh1, yoff1:yoff1 + ysh1
            ] = self.reg1.arrtemp.copy()
        self.reg1.arrtemp = arr1

        # fill mask1 with data in an appropriate position
        mask1[:, xoff1:xoff1 + xsh1, yoff1:yoff1 + ysh1
            ] = self.reg1.arrmask.copy()
        self.reg1.arrmask = mask1

        # fill arr2 with data in an appropriate position
        arr2[:, xoff2:xoff2 + xsh2, yoff2:yoff2 + ysh2
            ] = self.reg2.arrtemp.copy()
        self.reg2.arrtemp = arr2

        # fill mask2 with data in an appropriate position
        mask2[:, xoff2:xoff2 + xsh2, yoff2:yoff2 + ysh2
            ] = self.reg2.arrmask.copy()
        self.reg2.arrmask = mask2

#        print("xoff1, yoff1: ", xoff1, yoff1)
#        print("xoff2, yoff2: ", xoff2, yoff2)

        xoff1 = xshape - xoff1 - xsh1
#        yoff1 = yshape - yoff1 - ysh1
        xoff2 = xshape - xoff2 - xsh2
#        yoff2 = yshape - yoff2 - ysh2

        # update reference points
        for pt in self.reg1.ptstemp:
#            pt += [0, yoff1, xoff2] # WTF?!
            pt += [0, yoff1, xoff1] # WTF?!
#            pt += [0, yoff2, xoff2] # WTF?!
#            pt += [0, xoff1, yoff1] # WTF?!

        for pt in self.reg2.ptstemp:
#            pt += [0, yoff2, xoff1] # WTF?!
#            pt += [0, yoff1, xoff1] # WTF?!
#            pt += [0, yoff2, xoff2] # WTF?!
            pt += [0, yoff2, xoff2] # WTF?!
#            pt += [0, xoff2, yoff2] # WTF?!

#        print("xoff1: ", xoff1)

        # update control points
        if ctrlpts:
            for pt in ctrlpts[0]:
                pt[0] += 0
                pt[1] += yoff2 # WTF?
#                pt[1] += yoff1 # WTF?
                pt[2] += xoff2 # WTF?
                
            for pt in ctrlpts[1]:
                pt[0] += 0
                pt[1] += yoff2 # WTF?
#                pt[2] += xoff1 # WTF?
                pt[2] += xoff2 # WTF?
                
    # GRAPHICAL OUTPUT

    def inspect_blocks(self, refpts=True, data="temp",
        cutlist=[], preprocess=True):
        """inspect_channels routine modified for use with data
        from registration blocks
        
        refpts - if True, then reference points are shown as well
        data - which data to shown; three options possible:
            "temp" - arrtemp
            "orig" - arrorig
            "mask" - arrmask
        cutlist - auxiliary control points used e.g. by routine
            'find_parameters'; structure of cutlist is the 
            same as that of 'controlpoints' argument in routine
            'find_parameters'
        preprocess - if True, then routines 'apply_transf' and
            'make_identical_dims' are used to show reconstructed
            forms of arrays; if False, then arrays are shown how
            they look like with no transformations beign performed
            
        returns:
        adjpoints - adjusting points that can be used for finding
            the best match
        """
        
        cutlist = copy.deepcopy(cutlist)
        
        # preprocess data
        if preprocess: self.reconstruct_arrays(ctrlpts=cutlist)
        
        # data to be shown
        if data == "temp":
            data1 = self.reg1.arr_with_nans
            data2 = self.reg2.arr_with_nans
        elif data == "orig":
            data1 = self.reg1.arrorig
            data2 = self.reg2.arrorig
        elif data == "mask":
            data1 = self.reg1.arrmask
            data2 = self.reg2.arrmask
        else:
            print("inspect_blocks: Unknown type of data to show.")
            return

        # names of blocks as depicted by inspect_channels
        NAME1, NAME2 = "reg1", "reg2"
        
        # inspect_channels requires SPMdata-like structures as
        # inputs, so we create auxiliary class with atributes 
        # necessary for inspect_channels
        class AuxSPM:
            def __init__(self):
                self.channels = {NAME1:data1, NAME2:data2}
                self.xsteps = {NAME1:1,  NAME2:1}
                self.ysteps = {NAME1:1,  NAME2:1}
                self.units  = {NAME1:"", NAME2:""}
                self.numlay = data1.shape[0]
        
        # create SPMdata-like object
        aux = AuxSPM()
        
        # if refpts, include reference points
        if refpts:
            ctrlpts = [self.reg1.ptstemp, self.reg2.ptstemp]
        else:
            ctrlpts = []
        
        # if cutlist is not empty, include these control points
        if cutlist:
            if ctrlpts:
                ctrlpts[0].extend(list(cutlist[0]))
                ctrlpts[1].extend(list(cutlist[1]))
            else:
                ctrlpts = cutlist

        if ctrlpts:
            lenctrl = len(ctrlpts)
        else:    
            lenctrl = 0
        
        # show blocks
        inspect_channels([aux, aux], [NAME1, NAME2],
            showcut=True, cutlist=ctrlpts, showconstants=True,
            transposectrlpts=False, scaling=True)#,
            #adjustorient=True)
            
        # adjusting points; if additional points are chosen by
        # mouse in the environment of inspect_channels, the
        # variable ctrlpts gets updated by coords. of these points
        # so we take only newly added points
        adjpoints = ctrlpts[lenctrl:]
            
        # if data were reconstructed by setting
        # preprocess=True, then these should be "de-reconstructed"
        if preprocess: self.dereconstruct_arrays()
            
        return adjpoints

    def show_overlap(self, rlay=0, grid=True, show=True):
        """from reg1 and reg2 arrtemp's retrieve only such an
        area, for which both reg1.arrmask and reg2.arrmask yields
        True, i.e. take overlapping area of reg1.arrtemp and
        reg2.arrtemp and return it
        
        rlay - layer to process and show
        grid - if True, then grid is shown over the images
        show - if True, the result shown
        
        returns:
        arr1, arr2 - areas of reg1 and reg2 arrtmps, respectively,
            that correspond to overlapping masks
        """
        
        # reconstruct arrtmps so they have proper dimensions and
        # mutual positions of valid data
        self.reconstruct_arrays()
        
        # preprocess masks
        mask1 = self.reg1.arrmask > MASK_LIM 
        mask2 = self.reg2.arrmask > MASK_LIM
        
        # find area where both masks yields True
        m = mask1 * mask2
        
        # find lower and upper bounds for the area
        mrow, mclm = m[rlay].any(axis=1), m[rlay].any(axis=0)
#        mrow, mclm = m[rlay].all(axis=1), m[rlay].all(axis=0)
        minrow = np.argmax(mrow)
        minclm = np.argmax(mclm)
        maxrow = len(mrow) - np.argmax(mrow[::-1])
        maxclm = len(mclm) - np.argmax(mclm[::-1])

        # extract valid data from arrtmps
        arr1 = self.reg1.arr_with_nans[rlay,
            minrow:maxrow, minclm:maxclm]
        arr2 = self.reg2.arr_with_nans[rlay,
            minrow:maxrow, minclm:maxclm]

        # since we reconstructed data, we should "de-reconstruct"
        # them again
        self.dereconstruct_arrays()

        # create figure and axes
        fig, ax = plt.subplots(1, 2,
            num="Common area of both blocks")
        
        # show images
        if show:
            # create images
            ax[0].imshow(arr1, origin="lower")
            ax[1].imshow(arr2, origin="lower")
            
            # add grid over images
            if grid:
                ax[0].grid(which='major', axis='both', ls='-')
                ax[1].grid(which='major', axis='both', ls='-')
            
            # show
            plt.show()
        
        return arr1, arr2

    def show_match(self, rlay=0, threshold=0.5, alpha=0.5):
        """show both blocks of data to compare how the two are
        similar to each other
        
        rlay - which layer to show
        threshold - the threshold for values in arrays, under
            which the given value is "truncated" to zero, and
            above which the value is set to one; if None, then 
            no threshold is used
        alpha - the amount of opacity of the overlap
            
        BEWARE: if imshow window is closed or resized (by mouse?),
            then program run is aborted by en error, so it is wise
            to let the 'overlaps' window untouched during the
            calculation
        """

        arr1 = self.reg1.arr_with_nans[rlay].copy()                
        arr2 = self.reg2.arr_with_nans[rlay].copy()
        
        if self.fig is None:
            self.fig, self.ax = plt.subplots(1, num="Overlap")
            self.fig.show()    

        show_match(arr1, arr2, rlay=rlay, threshold=threshold,
            alpha=alpha, fig=self.fig, ax=self.ax)

    # FINDING THE BEST MATCH

    def find_match(self, pardict=None, maxsteps=1e1,
        tolerancy=1e-4, factor=2, facfreq=2,
        showoverlap=True, threshold=None, steptolerancy=1e-2):
        """iterate a sequence of transformations to find the best
        match of registration blocks; this routine assumes that
        data to match are close to each other already;
        in each step of iteration that is multiple of 'facfreq',
        the parameter limits and steps are divided by 'factor';
        initial values of parameters are stored in pardict;
        
        unlike weigh_corr, here we do not care about
        non-commutativity of transformations, so parlist contains
        only a single tuple with no need of post-processing
        
        pardict - dictionary of parameter values of the form
            {"zoom":(zoomneigh, zoomstep), ...}, where
            zoomneigh is a number of steps in the neighbourhood
            of current values to investigate and zoomstep is an
            initial step;           
            dictionary of parameter values may also be of the form
            {"zoom":zoomstep, ...}, in which case by default
            zoomneigh is put equal to 1;
            if some or all necessary parameters are missing or 
            pardict is None, then default values from
            self.matchpars are used
        maxsteps - maximum number of iterations
        tolerancy - threshold to determine whether the corr.
            coefficient in a given step of iterations attains its
            maximum
        steptolerancy - threshold to determine whether to stop
            iterations
        showoverlap - if True, then overlap of the two RegBlocks
            is plotted after each iteration step
        factor - factor by which all parameter ranges are divided
            in order to refine ranges for which correlation is 
            calculated; values are refined each 'facfreq' steps
        facfreq - integer stating how often values of parameters
             are refined; e.g. if facfreq=3, then values are
             divided by 'factor' every third iteration
        threshold - see 'threshold' parameter of the routine
            'show_match'
        """
        
        TOL_ZERO = 1e-5
        NUM_STEPS = 1

        # we do not put default value of 'pardict' to {},
        # because in that case multiple callings of 'find_match'
        # update 'pardict', so when 'find_match' is called for the
        # second time, the default form of 'pardict' is no longer
        # empty! we prevent this by setting default value of 
        # pardict to None
        if pardict is None: pardict = {}
        
        # retrieve initial values of parameters
        for key, val in pardict.items():
            if isinstance(val, int) or isinstance(val, float):
                pardict[key] = (NUM_STEPS, val)
            elif len(list(val)) == 2:
                pass
            else:
                print("find_match: Invalid form of "
                     "parameter dictionary. Halt.")
                return
            
        # process input values for all parameters and create 
        # sequences of values; if some or all values are not
        # specified, use the default ones
        lims = {}
        
        print("pardict: ", pardict)
        
        defaultvalues = self.matchpars
        for key, (dneigh, val) in defaultvalues.items():

            if key not in pardict.keys():   
                pardict[key] = (dneigh, val)
                print("find_match: Values (neigh, val) for "
                    "'{:6}' are not specified, so default values "
                    "{:10} are used instead.".format(
                    key, str((dneigh, val))))
                
            neigh, step = pardict[key]
            lims[key] = np.linspace(-neigh * step, neigh * step,
                num=(1 + 2 * neigh))
                
        # get sequences of parameter values to be searched through
        zoomlim,   rotatelim = lims["zoom"],   lims["rotate"]
        xshiftlim, yshiftlim = lims["xshift"], lims["yshift"]
        xskewlim,  yskewlim  = lims["xskew"],  lims["yskew"]

        # NOTE: up to now only initial values were preprocessed;
        # nevertheless, which transformations should be applied 
        # to the data is determined by limtuple below, which in 
        # turn is determined by the order of transformations 
        # present in the routine apply_transf
        limtuple = (zoomlim, rotatelim, xshiftlim, yshiftlim)

#        print("limtuple: ", limtuple)

        # two auxiliary lists for storing temporary data
        parlisttmp  = np.zeros(len(limtuple), dtype=float)
        corrlisttmp = np.zeros(len(limtuple), dtype=float)

        # reset temporary arrays in registration blocks
        self.reg1.reset_temp()
        self.reg2.reset_temp()
        self.reconstructed = False
        
        # reconstruct arrays in order to calculate initial corr.
        # coeff.
        self.reconstruct_arrays()

        # calculate correlation coefficient
        corrold = self.compfun(self.reg1, self.reg2)       

#        # value ranges for clipping
#        clip1 = self.reg1.value_range()
#        clip2 = self.reg2.value_range()

        # show overlap
        if showoverlap: self.show_match(threshold=threshold)

        # dereconstruct arrays so they can be used anew in a loop
        self.dereconstruct_arrays()

        # initial values of parameters
        pars = copy.deepcopy(self.reg1.parstemp)
        pars = np.array(pars, dtype=float)

        # iterations are performed until maxsteps is reached
        # or the difference between values of parameters of
        # two successive steps is below 'steptolerancy'
        for i in range(int(maxsteps)):
            parlisttmp  = np.zeros(len(limtuple), dtype=float)
            corrlisttmp = np.zeros(len(limtuple), dtype=float)

            print(("find_match: Step no. {} with corr. "
                "coeff. = {:.5}.").format(i, corrold))

            # for each transformation calculate corresponding
            # value of corr. coefficient
            for j, lim in enumerate(limtuple):
                corrlist = np.zeros(len(lim), dtype=float)

                # for a given transformation get corr. coeff. 
                # for every choice of parameter value in range 
                for k, val in enumerate(lim):
                    # reset 'parstemp' to original values
                    parstemp = pars.copy()
                    
                    # modify j-th parameter by value 'val'
                    parstemp[j] += val
                    
                    # reconstruct array when all transformations
                    # with parameters 'parstemp' are applied to it
                    self.reconstruct_arrays(parstemp)
                    
                    # calculate corr. coeff. 
                    corrlist[k] = self.compfun(
                        self.reg1, self.reg2)
                        
                    # reset arrays to their original forms
                    self.dereconstruct_arrays()

                # corr. coeff. of the best match
                corrlisttmp[j] = np.nanmax(corrlist)
                
                # incremented value of j-th parameter for which
                # the best match was obtained
                parlisttmp[j] = lim[np.nanargmax(corrlist)]

            # how corr. coeff. changed from previous iteration;
            # take only positive changes, negative differences are
            # discarded
            corrlisttmp = np.clip(corrlisttmp - corrold, 0.0, 1.0)

            # normalize 'corrlisttmp'
            if np.nansum(corrlisttmp) < TOL_ZERO:
                pass
            else:
                corrlisttmp /= np.nansum(corrlisttmp)
                
            # weight parameter increments 'parlisttmp' by
            # 'corrlisttmp' and add these increments to original
            # values of parameters
            pars += parlisttmp * corrlisttmp

            # reconstruct the form of arrays when all transforms
            # are applied with the best-match values weighted by
            # corresponding corr. coeff.
            self.reconstruct_arrays(pars)

            # show overlap
            if showoverlap: self.show_match(threshold=threshold)

            # calculate correlation coefficient
            corr = self.compfun(self.reg1, self.reg2)
            
            # reset arrays to their original forms
            self.dereconstruct_arrays()
            
            # if corr. coef. is no more increasing, make the 
            # steps finer by dividing them by 'factor'
            corrdiff = np.abs(corr / corrold - 1)
            corrold = corr
            
            if corrdiff < tolerancy:
                zoomlim   /= factor; rotatelim /= factor
                xskewlim  /= factor; yskewlim  /= factor
                xshiftlim /= factor; yshiftlim /= factor
                                  
            # if the steps are too small, stop iterations
            if any([abs(lim[0]) < steptolerancy
                for lim in limtuple]):
                print(("find_match: Iterations terminated "
                    "in step no. {} with corr. coeff. "
                    "= {:.5}.").format(i, corr))
                break
        else:
            print(("find_match: Maximum number {} of steps "
                "reached. Current value of correlation "
                "coefficient reads: corr = {:.5}.").format(
                maxsteps, corr))
                
        # store final values of parameters
        self.reg1.parstemp = pars.copy()
        self.reg1.parsorig = pars.copy()
        
        print("self.reg1.parstemp: ", self.reg1.parstemp)
        
        # reset figure and axes
        self.reset_figure()

    def reset_figure(self):
        """effectively destroy auxiliary window; this routine is
        especially helpful when killing routine 'find_match', 
        which shows real-time overlaps; in such a case the figure
        needs to be destroyed manually
        """
#        plt.close(self.fig)        
        self.fig, self.ax = None, None

    def find_parameters(self, controlpoints):
        """find the best match for control points and apply
        appropriate transformations on registration blocks;
        transformations are done in x-y plane only! more general
        method should be written...
        
        controlpoints - sequence of the form
            controlpoints = [reg1list, reg2list], where
            reg1list = [[zcoord1, xcoord1, ycoord1], ...,
                [zcoordN, xcoordN, ycoordN]], where
            N >= 2 and zcoord1 is z-coordinate of the first
            control point for reg1, zcoord2 is z-coord. of the
            second control point for reg1, etc.
            and analogously for reg2
        """
        
        shape1 = self.reg1.arrtemp.shape
        shape2 = self.reg2.arrtemp.shape
        pars = find_parameters(shape1, shape2, controlpoints)
        
        # this part specifies how reg1 has to be transformed in
        # order to get aligned with reg2
        self.reg1.parsorig = pars.copy()
        self.reg1.parstemp = pars.copy()

    # EXPORT AND IMPORT

    def export_data(self, name="regist_data",
        headingsuffix=".hdr", datasuffix=".dat"):
        """export data to the external files; two files are 
        created:
            text file 'name'.'headingsuffix' contains all
            metadata
            binary file 'name'.'datasuffix' contains all data,
            i.e. arrorigs, arrtmps, points and masks for both
            blocks       
        """
        
        headingname = name + headingsuffix
        filename = name + datasuffix

        # construct attrlist
        attrlist = [
        "# METADATA FOR REGISTRATION DATA\n",
        "# coresponding raw data file", filename,
        "# description of reg1 block", self.reg1.description,
        "# description of reg2 block", self.reg2.description,
        "# comparing function", self.compfun.__name__,
        "\n"
        ]

        # arraylist has to have the same ordering
        # as arraylist in import_data routine!
        arraylist = [
            ("reg1", "ptsorig"), ("reg1", "ptstemp"),
            ("reg2", "ptsorig"), ("reg2", "ptstemp")]

        # take values of reference points
        for regnm, datnm in arraylist:
            string = "# reference points " + regnm + "." + datnm
            vallist = ""
            
            reg = getattr(self, regnm)
            triple = getattr(reg, datnm)
            for array in triple:
                point = list(map(str, array))
                vallist += ("\t".join(point) + "\n")

            attrlist.extend([string, vallist])

        # orlist has to have the same ordering as orlist in
        # import_data routine!
        orlist = [
            ("reg1", "oriorig"), ("reg1", "oritemp"),
            ("reg2", "oriorig"), ("reg2", "oritemp")]

        # take values of origin points
        for regnm, datnm in orlist:
            string = "# origin points " + regnm + "." + datnm

            reg = getattr(self, regnm)
            data = getattr(reg, datnm)
            val = list(map(str, data))
            val = "\t".join(val) + "\n"
            
            attrlist.extend([string, val])

        # datalist has to have the same
        # ordering as datalist and in import_data routine!
        datalist = [
            ("reg1", "arrorig"), ("reg1", "arrtemp"),
            ("reg1", "arrmask"), ("reg2", "arrorig"),
            ("reg2", "arrtemp"), ("reg2", "arrmask")]

        # take shapes of arrays
        for regnm, datnm in datalist:
            string = "# shape of array " + regnm + "." + datnm

            reg = getattr(self, regnm)
            data = getattr(reg, datnm)
            shape = "\t".join(list(map(str, data.shape)))
            
            attrlist.extend([string, shape])

        attrlist.append(" ")

        # parslist has to have the same
        # ordering as parslist in import_data routine!
        parslist = [
            ("reg1", "parsorig"), ("reg1", "parstemp"),
            ("reg2", "parsorig"), ("reg2", "parstemp")]

        # parameters for reconstruction
        for (regnm, datnm) in parslist:
            string = "# parameters for reconstruction " 
            string += regnm + "." + datnm
            
            reg = getattr(self, regnm)
            pars = getattr(reg, datnm)
            pars = "\t".join(list(map(str, pars)))
            attrlist.extend([string, pars])

        # save attrlist to heading text file
        with open(headingname, "w") as f:
            f.write('\n'.join(attrlist))
        
        # save raw data for channels in binary data file
        with open(filename, "wb") as f:

            ENDIAN = '>'
            DATA_TYPE = 'f'

            for regnm, datnm in datalist:
                reg = getattr(self, regnm)
                val = getattr(reg, datnm)

                spec = ENDIAN + str(val.size) + DATA_TYPE
                data = struct.pack(spec, *val.flat)
                f.write(data)
            
        print("export_data: Raw data and metadata exported to"
            " files:\n\t{}\n\t{}".format(filename, headingname))

    def import_data(self, headingname, filename=None):
        """import data from external files
        
        headingname - name of the text file containing metadata
        filename - name of the binary file containing raw data;
            if None, then this name is taken from the name
            stored in headingname; if None, it assumes that binary
            file is in the same folder as text file
        """

        # --- HEADING ---
        
        # import from the text file
        try:
            with open(headingname, "r") as f:
                heading = f.read().splitlines()
        except FileNotFoundError:
            print("import_data: No text file to be found. Halt.")
            return

        # name of the file with binary data
        if filename is None:
            # get path of the directory
            namedir = os.path.dirname(headingname)
            
            # get name of the binary file
            filename = heading[3]
            
            # prepend the path to the name
            filename = os.path.join(namedir, filename)

        # import descriptions
        self.reg1.description = heading[5]
        self.reg2.description = heading[7]
        
        # import the name of comparing function
        compfun = heading[9]
        
        if compfun in globals():
            self.compfun = globals()[compfun]
            print(("import_data: Comparing function is chosen as "
                "'{}'.").format(compfun))                
        else:
            print(("import_data: Unknown comparing function "
                "'{}'. Halt.").format(compfun))
            return

        initname = "# reference points reg1.ptsorig"
            
        arraylist = [
            ("reg1", "ptsorig"), ("reg1", "ptstemp"),
            ("reg2", "ptsorig"), ("reg2", "ptstemp")]
            
        # initial index
        initidx = heading.index(initname) + 1
        
        # number of lines containing data about one set of points
        # i.e. Lines Per Set
        LPS = 5
        
        # extract reference points...
        for i, (regname, arrname) in enumerate(arraylist):
            pt1 = heading[LPS * i + initidx + 0].split("\t")
            pt2 = heading[LPS * i + initidx + 1].split("\t")
            pt3 = heading[LPS * i + initidx + 2].split("\t")

            pt1 = np.array(list(map(float, pt1)))
            pt2 = np.array(list(map(float, pt2)))
            pt3 = np.array(list(map(float, pt3)))
            
            triple = [pt1, pt2, pt3]
            
            # get registration block
            reg = getattr(self, regname)
            
            # set imported values of reference points to reg
            setattr(reg, arrname, triple)

        orlist = [
            ("reg1", "oriorig"), ("reg1", "oritemp"),
            ("reg2", "oriorig"), ("reg2", "oritemp")]

        initidx += LPS * (i + 1)
        LPS = 3

        # extract origin points...
        for i, (regname, arrname) in enumerate(orlist):
            pt = heading[LPS * i + initidx].split("\t")
            pt = np.array(list(map(float, pt)))
            
            # get registration block
            reg = getattr(self, regname)
            
            # set imported values of reference points to reg
            setattr(reg, arrname, pt)

        datanameslist = [
            ("reg1", "arrorig"), ("reg1", "arrtemp"),
            ("reg1", "arrmask"), ("reg2", "arrorig"),
            ("reg2", "arrtemp"), ("reg2", "arrmask")]

        # new initial index
        initidx += LPS * (i + 1)

        # number of lines containing data about one dataset
        # i.e. Lines Per Set
        LPS = 2
        shapelist = [None]*len(datanameslist)

        # extract shapes of arrays...
        for i, (regname, arrname) in enumerate(datanameslist):
            shape = heading[LPS * i + initidx].split("\t")
            shapelist[i] = list(map(int, shape))

        # extract parameters for reconstruction...
        parslist = [
            ("reg1", "parsorig"), ("reg1", "parstemp"),
            ("reg2", "parsorig"), ("reg2", "parstemp")]

        # new initial index
        initidx += LPS * (i + 1) + 1

        for i, (regname, arrname) in enumerate(parslist):        
            pt = heading[LPS * i + initidx].split("\t")
            pt = list(map(float, pt))
            
            # get registration block
            reg = getattr(self, regname)
            
            # set imported values of parameters to reg
            setattr(reg, arrname, pt)
    
        # --- BINARY DATA ---

        ENDIAN = '>'
        DATA_TYPE = 'f'
        
        datasize = [np.prod(shape) for shape in shapelist]
        spec = ENDIAN + str(sum(datasize)) + DATA_TYPE

        # import
        try:
            with open(filename, "rb") as f:
                data = f.read()
        except FileNotFoundError:
            print("import_data: No binary file to be found. "
                "The import is incomplete. Halt.")
            return

        data = struct.unpack(spec, data)
        data = np.array(data)

        # update individual arrays
        low = 0
        for i, (regname, arrname) in enumerate(datanameslist):
            high = low + datasize[i]
            arr = data[low:high]
            arr = arr.reshape(shapelist[i])
            
            reg = getattr(self, regname)
            setattr(reg, arrname, arr.copy())            
            low = high
        
        print(("import_data: Metadata and raw data successfully"
            " imported from files:\n\t{}\n\t{}."
            ).format(headingname, filename))

