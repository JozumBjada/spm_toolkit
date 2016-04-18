from spm_data_cut1d import *
import copy
from mpl_toolkits.axes_grid1 import make_axes_locatable    
from scipy.ndimage import rotate as ndrotate


def _modify_pos(nvec, npos, arrshape):
    """return control points for cut through entire plane;
    method used is quite crude, since it does not determine
    intersection of cutting plane and data cuboid, but
    just adds the maximum possible range, i.e. space diagonal,
    in both directions; resulting cut has to be then strip_nans-ed
    
    nvec - normal vector of the form (x, y, z)
    npos - position vector of the form (x, y, z)
    arrshape - shape of the data cuboid
    """
    
    # normal vector of the cutting plane
    nvec /= np.linalg.norm(nvec)
    n1, n2, n3 = nvec

    # first vector generating the cutting plane
    if abs(n3) < 1e-9:
        vec1 = np.array([0, 0, 1])
    else:
        vec1 = np.array([n3, 0, -n1])
        vec1 /= np.linalg.norm(vec1)

    # second vector generating the cutting plane
    vec2 = np.cross(nvec, vec1)
    
    # space diagonal of a data cuboid
    diag = np.linalg.norm(arrshape)

    # control points
    point1 = npos - diag*vec1 - diag*vec2
    point2 = npos - diag*vec1 + diag*vec2
    point3 = npos + diag*vec1 + diag*vec2

#        print("vec1, vec2: ", vec1, vec2)
#        print("point1, point2, point3: ", point1, point2, point3)    

    return point1, point2, point3

def average_layers(spmdata, laylist=None):
    """average layers of SPMdata structure and return Cut2D structure containing results
    
    spmdata - input SPMdata structure
    laylist - which layers to choose for averaging; if None, then all layers are used
    """
    
    res = Cut2D()
    laylist = range(spmdata.numlay) if laylist is None else laylist
    for chan, val in spmdata.channels.items():
        res.channels[chan] = val[laylist].sum(axis=0) / len(laylist)
        res.minmax[chan] = (np.nanmin(val), np.nanmax(val))

    return res
    

class Cut2D:

    def __init__(self):
        self.channels = {}
        self.posl = []
        self.minmax = {}

    def rotate(self, angle):
        """rotate cuts counter-clockwise by 'angle' specified in degrees
        """
        
        for chan in self.channels:
            self.channels[chan] = ndrotate(self.channels[chan], angle, cval=np.nan)

    def _process_cut_coords(self, posx, posy, posz, delta, deltaz, arrshape=None):
        """based on control points stored in posz, posx, posy create sampling points
        and store them in 2D arrays posxl, posyl, poszl
        
        posx - list containing x-coordinates of three control points
        posy - list containing y-coordinates of three control points
        posz - list containing z-coordinates of three control points
        delta - approx. step in x- and y-directions
        deltaz - approx. step in z-direction
        arrshape - shape (3-tuple) of the data cuboid; if None, then 
            cutting is done in the selection parallelogram, otherwise
            the cut is performed through the whole cutting plane
        """
           
        minlen = min(len(posx), len(posy), len(posz))
        
        if minlen < 3:
            print("_process_cut_coords: Not enough control points.")
            return [], [], []
        
        if minlen != max(len(posx), len(posy), len(posz)):
            print("_process_cut_coords: Unequal lengths of coordinate lists. "
                "Only first 3 elements used.".format(minlen))

        # reorder the three control points so that selection is independent of
        # the original order picked by e.g. inspect_channels (assumes the middle
        # point of selection is also the middle point in lists):
        point1, point2, point3 = np.array([posx[:3], posy[:3], posz[:3]]).T        
        point1, point2, point3, nvec = vector_handedness(point1, point2, point3)            
        if arrshape is not None:
            point1, point2, point3 = _modify_pos(nvec, point2, arrshape)
            point1, point2, point3 = vector_handedness(point1, point2, point3)[:3]
        posx, posy, posz = np.array([point1, point2, point3]).T

#        print("point1, point2, point3: ", point1, point2, point3)    
#        print("ha, nvec: ", np.cross(point1 - point2, point3 - point2), nvec)
       
        # unlike in _process_cut_coords for Cut1D, here these lists contain
        # two lists of coordinates
        posxl, posyl, poszl = [], [], []   
        
        # sampling of lines connecting control points
        for i in range(2):
        
            # how many sampling points in each dimension?
            auxnumz = round(abs(posz[i+1] - posz[i])/deltaz)
            auxnumx = round(abs(posx[i+1] - posx[i])/delta)
            auxnumy = round(abs(posy[i+1] - posy[i])/delta)
            
            # we need the same number of sampling points in each dimension
            num = int(max(auxnumx, auxnumy, auxnumz)) # int(...) here is important
            
            if num == 0:
                print("_process_cut_coords: Same control points.")
                continue
            
            # calculate consistent sampling steps
            dx = (posx[i+1] - posx[i]) / num
            dy = (posy[i+1] - posy[i]) / num
            dz = (posz[i+1] - posz[i]) / num
             
            # calculate sampling points
            auxlz = [posz[i] + dz*j for j in range(num)]
            auxlx = [posx[i] + dx*j for j in range(num)]
            auxly = [posy[i] + dy*j for j in range(num)]

            # update lists of points
            poszl.append(auxlz)
            posxl.append(auxlx)
            posyl.append(auxly)

        fstline = np.array([poszl[0], posxl[0], posyl[0]]) # first line
        scnline = np.array([poszl[1], posxl[1], posyl[1]]) # second line
        deltaline = scnline - scnline[:,0, np.newaxis] # delta line
    
        array = fstline.T[:, :, np.newaxis] + deltaline
        return array[:,1,:], array[:,2,:], array[:,0,:], point1, point2, point3

    def create_cut_rotate(self, spmdata, posz, posx, posy):
        """create 2D cut through spmdata, cut covers the whole intersection
        of the cutting plane and spmdata cuboid
        
        spmdata - SPMdata instance through which the cut is done
        posx - x-coordinate 2D array of sampling points
        posy - y-coordinate 2D array of sampling points
        posz - z-coordinate 2D array of sampling points
        """

        minlen = min(len(posx), len(posy), len(posz))
        
        if minlen < 3:
            print("create_cut_rotate: Not enough control points.")
            return
        
        if minlen != max(len(posx), len(posy), len(posz)):
            print("create_cut_rotate: Unequal lengths of coordinate lists. "
                "Only first 3 elements used.".format(minlen))

        # reorder the three control points so that selection is independent of
        # the original order picked by e.g. inspect_channels (assumes the middle
        # point of selection is also the middle point in lists):
#        points = np.array([posz[:3], posx[:3], posy[:3]]).T # each row corresponds to one point
        points = np.array([posx[:3], posy[:3], posz[:3]]).T # each row corresponds to one point
        vec21 = points[0] - points[1]
        vec23 = points[2] - points[1]
        
        # normal to the cutting plane and its base
        nvec = np.cross(vec21, vec23)
        nvec /= np.linalg.norm(nvec)
        npos = points[1]
#        print("nvec: ", nvec)        
        print("nvec, npos: ", nvec, npos)
        
        zvec = np.array([0, 0, -1])
#        zvec = np.array([-1, 0, 0])

        # if vectors connecting the three control points form with z-vector
        # a left-hand basis, then reorder them 
        if np.dot(nvec, zvec) < 0:
            print("create_cut_rotate: Control points formed a left-hand basis,"
                " so they were reordered.")
            nvec *= -1

        self.xpos = None
        self.ypos = None
        self.zpos = None

        x, y, z = nvec

        rho = np.arccos(z) * 360 / (2 * np.pi)          
        if abs(rho % np.pi) < 1e-9:
            phi = 0
        else:
            phi = np.pi / 2 if abs(x) < 1e-9 else np.arctan(y / x)                     
            phi *= 360 / (2 * np.pi)          
        
        print("x, y, z: ", x, y, z)
        print("phi, rho: ", phi, rho)
        
        self.minmax = {chan: (np.nanmin(val), np.nanmax(val)) for chan, val in spmdata.channels.items()}

#        npos = 0 # !!!
                
        for chan, arr in spmdata.channels.items():

            par = 0
            for z in [0, arr.shape[0]]:
                for x in [0, arr.shape[1]]:
                    for y in [0, arr.shape[2]]:
#                        parloc = np.dot(nvec, npos - np.array([x, y, z]))
                        parloc = np.dot(nvec, npos + np.array([x, y, z]))
                        if parloc > par: par = parloc
            
#            par = spmdata.numlay - par
            print("par: ", par)
            
            # TODO
            
            plt.imshow(arr[par])
            plt.show()
            
            arr[np.isnan(arr)] = 0
            arr = ndrotate(arr, phi, axes=(2,1))#, cval=np.nan) #(1,0)
#            print("arr.shape: ", arr.shape)

            plt.imshow(arr[par])
            plt.show()

            arr = ndrotate(arr, rho, axes=(0,1))#, cval=np.nan)
            self.channels[chan] = arr[par]

            print("arr.shape: ", arr.shape)
            plt.imshow(arr[par])
            plt.show()



    def create_cutzxy_raw(self, spmdata, posz, posx, posy, interpolate=0):
        """create 2D cut through spmdata
        
        spmdata - SPMdata instance through which the cut is done
        posx - x-coordinate 2D array of sampling points
        posy - y-coordinate 2D array of sampling points
        posz - z-coordinate 2D array of sampling points
        interpolate - the order of interpolation ranging from 0 to 5
        """
        
#        minlen = min(len(posx), len(posy), len(posz))
#        if minlen == 0:
#            print("create_cutzxy_raw: No positions specified.")
#            return
#        
        shapes = set([posx.shape, posy.shape, posz.shape])
        if len(shapes) != 1:
            print("create_cutzxy_raw: Unequal shapes of coordinate arrays.")
            return
                
        shapes = shapes.pop()
        if len(shapes) != 2:
            print("create_cutzxy_raw: Coordinate arrays are not 2D arrays.")
            return
            
        self.minmax = {chan: (np.nanmin(val), np.nanmax(val)) for chan, val in spmdata.channels.items()}
        
        self.posx = posx
        self.posy = posy
        self.posz = posz        

        for chan, val in spmdata.channels.items():
            # BACHA, BACHA, BACHA!!! i zde se projevuje prehozeni xove a yove souradnice
            self.channels[chan] = map_coordinates(val,
                np.array([self.posz, self.posy, self.posx]), 
                order=interpolate, mode='constant', cval=np.nan)
            # nearest-neighbour approximation
#            self.channels[chan] = val[self.posz, self.posy, self.posx]        
        
    def align_channels(self, point1, point2, point3, interpolation=1):
        """pad array for all channels by NANs, NANs are padded
        in such a way that resulting array has the same geometrical
        shape as the original cut of spm data
        
        point1, point2, point3 - control points,
            for 2D cut there are three control points, the first one
            marks the origin of the cutting parallelogram
        interpolation - order of interpolating polynomial, it can
            range from 0 to 5, but for order >= 2 I encountered problems
            for parallelograms determined by segments whose angle was
            greater than 90 degrees
        """

        # compute cosphi, i.e. cosine of the angle between line segments
        # determined by control points
        veca = point1 - point2
        vecb = point3 - point2
        cosphi = np.dot(veca, vecb) / (np.linalg.norm(veca) * np.linalg.norm(vecb))
        
        # align each channel
        for chan, val in self.channels.items():
            xnum, ynum = val.shape            

            dynum = xnum * abs(cosphi)
            auxf = np.empty((xnum, ynum + math.floor(dynum)))
            auxf.fill(np.nan)

            if cosphi >= 0:            
                fun = lambda i: dynum*(xnum - i - 1)/xnum
            else:
                fun = lambda i: dynum*(i)/xnum

            # x and y coordinates
            px = np.fromfunction(lambda i, j: i,
                (xnum, ynum), dtype=float)
            py = np.fromfunction(lambda i, j: j - (fun(i) % 1),
                (xnum, ynum), dtype=float)
                
            # new values
            aux = map_coordinates(val, [px, py], order=interpolation, cval=np.nan)
            # OPET DILEMA x vs. y                

            # reallocate cut array
            for xi in range(xnum):
                yi = math.floor(fun(xi))
                auxf[xi, yi:yi + ynum] = aux[xi]

            # save reallocated array
            self.channels[chan] = auxf
            
#    def align_channels(self, posx, posy, posz, interpolation=1):
#        """pad array for all channels by NANs, NANs are padded
#        in such a way that resulting array has the same geometrical
#        shape as the original cut of spm data
#        
#        posx, posy, posz - x-, y- and z-coordinates of control points,
#            for 2D cut there are three control points, the first one
#            marks the origin of the cutting parallelogram
#        interpolation - order of interpolating polynomial, it can
#            range from 0 to 5, but for order >= 2 I encountered problems
#            for parallelograms determined by segments whose angle was
#            greater than 90 degrees
#        """

#        # compute cosphi, i.e. cosine of the angle between line segments
#        # determined by control points
#        points = np.array([posz, posx, posy]).T # each row corresponds to one point
#        veca = points[0] - points[1]
#        vecb = points[2] - points[1]
#        cosphi = np.dot(veca, vecb) / (np.linalg.norm(veca) * np.linalg.norm(vecb))
#        
#        # align each channel
#        for chan, val in self.channels.items():
#            xnum, ynum = val.shape            

#            dynum = xnum * abs(cosphi)
#            auxf = np.empty((xnum, ynum + math.floor(dynum)))
#            auxf.fill(np.nan)

#            if cosphi >= 0:            
#                fun = lambda i: dynum*(xnum - i - 1)/xnum
#            else:
#                fun = lambda i: dynum*(i)/xnum

#            # VERSION 1 - START            
##            aux = val
#            # VERSION 1 - END
#            
#            # VERSION 2 - START
#            px = np.fromfunction(lambda i, j: i,
#                (xnum, ynum), dtype=float)
#            py = np.fromfunction(lambda i, j: j - (fun(i) % 1),
#                (xnum, ynum), dtype=float)
#                
#            # OPET DILEMA x vs. y
#                
#            aux = map_coordinates(val, [px, py], order=interpolation)
##             VERSION 2 - END

#            for xi in range(xnum):
#                yi = math.floor(fun(xi))
#                auxf[xi, yi:yi + ynum] = aux[xi]

#            self.channels[chan] = auxf

        
    def create_cutzxy(self, spmdata, posz, posx, posy, delta=1, deltaz=1, interpolate=0, align=False, entire=True):
        """create 2D cut through spmdata
        
        spmdata - SPMdata instance through which the cut is done
        posx - x-coordinate list of control points
        posy - y-coordinate list of control points
        posz - z-coordinate list of control points
        interpolate - the order of interpolation ranging from 0 to 5
        delta - sampling step in x- and y-direction, e.g. 'delta' = 3
            means that approx. every third data point will be 
            taken, value of 'delta' is only an advice for the routine to 
            calculate actual sampling step, which might differ in the end
        deltaz - sampling step in z-direction, for more see 'delta'
        align - if True, arrays corresponding to channels are aligned,
            so that they are not skewed
        entire - if True, then the cut is performed through the entire bulk
            of data cuboid, otherwise the cut is a parallelogram defined by
            posx, posy, posz only 
        """
        
        if not posx or not posy or not posz:
            print("create_cutzxy: Not enough control points.")
            return
        
        # creation of sampling 2D arrays

        # 2D cut dimensions are determined with help of the shape of the first array
        # appearing in spmdata.channels.keys()
        shape = list(spmdata.channels.values())[0].shape if entire else None
        self.posx, self.posy, self.posz, *points = self._process_cut_coords(
                                                        posx, posy, posz, delta, deltaz, shape)

        # creation of the cut itself
        self.create_cutzxy_raw(spmdata, self.posz, self.posx, self.posy, interpolate=interpolate)  
        
        # aligning
        if align: self.align_channels(*points)

    def show_channels(self, *cl, rescale=True, interpolation='none'):
        """plot 2D cut of channels
        
        cl - list of channels to show; if empty, then all channels are chosen
        rescale - if False, imshow uses implicit colouring according to the data,
            if True, than plot is rescaled according to original spm data, so that
            colours of the cut are identical to that in the original data
        interpolation - how to interpolate values in imshow
        """
        
        # pick only valid channels
        cl = cl if cl else self.channels.keys()
        cl = [chan for chan in cl if chan in self.channels.keys()]
        if len(cl) == 0:
            print("show_channels: No valid channels.")
            return
            
        # create figure
        fig, subfiglist = plt.subplots(1, len(cl),
                            num="1D Cut - channel: " + ", ".join(cl))
        if len(cl) == 1: subfiglist = [subfiglist]
        
        # plot each valid channel
        for i, chan in enumerate(cl):            
            # use the scale of the original SMPdata?
            vmin, vmax = self.minmax[chan] if rescale else (None, None)
            
            # create image with title
            img = subfiglist[i].imshow(self.channels[chan], aspect='equal',
                vmin=vmin, vmax=vmax, interpolation=interpolation)                
            subfiglist[i].set_title("Channel: {}[units: {}]".format(chan, UNITS[chan]))
            
            # colorbar
            divider = make_axes_locatable(subfiglist[i])
            ax_cb = divider.append_axes("right", size="20%", pad=0.2)
            plt.colorbar(img, cax=ax_cb)
            
#            subfiglist[i].set_xlabel("No. of steps")
#            subfiglist[i].set_ylabel("{}".format(UNITS[chan]))
        
        # show
        fig.tight_layout()
        plt.show()

    def strip_nans(self, ref=None):        
        """strip NaNs from arrays corresponding to all channels, so that
        all valid data are preserved
        
        ref - reference channel which is used for determination
            of a clip frame
        """
            
        # reference channel
        ref = list(self.channels)[0] if ref is None else ref
        
        mask = ~np.isnan(self.channels[ref])
        mask0 = mask.any(axis=0)
        mask1 = mask.any(axis=1)
                
#        print("mask, mask0, mask1: ", mask, mask0, mask1)
#        print("ref: ", ref)
#        
        for chan in self.channels:
            self.channels[chan] = self.channels[chan][mask1][:, mask0]
#            print("shape 2: ", self.channels[chan].shape)

        return ref
        
    def export_data(self, filename):
        """export data to external file
        """
        
        with open(filename, "wb") as f:
            # TODO
            pass
            
        # nparray ma sve procedury pro export dat
            
        print("export_data: Data exported to file: {}".format(filename))
        
    def import_data(self, filename):
        """import data from external file
        """
        
        with open(filename, "rb") as f:
            buff = f.read().split(b"\x05\x89") # VYMYSLET NEJAKY VLASTNI ZNAK, NAPR. \x05\x89
            
        buff[0] # header
        
        # TODO

