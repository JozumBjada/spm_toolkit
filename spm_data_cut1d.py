#from spm_data_struct import *
from spm_data_base import *
from scipy.ndimage import map_coordinates

class Cut1D:

    def __init__(self):
        self.posx, self.posy, self.posz = None, None, None
        self.dx, self.dy, self.dz = None, None, None
        self.channels = {}
        
    def __len__(self):
        return len(self.posl)
        
    def print_channels(self):
        print("1D cut channels: ", list(self.channels.keys()))
        
    @property
    def poslist(self):
        return list(zip(self.posx, self.posy, self.posz))
        
    def _process_cut_coords(self, posz, posx, posy, delta, deltaz):
        """based on control points stored in posz, posx, posy create sampling points
        and store them in posxl, posyl, poszl
        
        delta - approx. step in x- and y-directions
        deltaz - approx. step in z-direction
        """
           
        minlen = min(len(posx), len(posy), len(posz))
        if minlen != max(len(posx), len(posy), len(posz)):
            print("_process_cut_coords: Unequal lengths of coordinate lists. "
                "Only first {} elements used.".format(minlen))

        posxl, posyl, poszl = [], [], []   
        for i in range(minlen - 1):
        
            # how many sampling point in each dimension?
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
            poszl.extend(auxlz)
            posxl.extend(auxlx)
            posyl.extend(auxly)
    
        return posxl, posyl, poszl
        
    def create_cutzxy(self, spmdata, posz, posx, posy, delta=1, deltaz=1, interpolate=0):
        """create 1D cut through spmdata
        
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
        """
        
        # sampling points creation
        self.posx, self.posy, self.posz = self._process_cut_coords(posz, posx, posy, delta, deltaz)

        # creation of the cut itself
        self.create_cutzxy_raw(spmdata, self.posz, self.posx, self.posy, interpolate=interpolate)    
        
    def create_zcut(self, spmdata, positions, deltaz=1, interpolate=0):
        """create 1D cut through spmdata is z-direction
        
        spmdata - SPMdata instance through which the cut is done
        positions - unlike in create_cut, in this function 'positions' is
            a sequence of 2 numbers [x, y] specifying position of the z-cut
        delta - sampling step in x- and y-direction, e.g. 'delta' = 3
            means that approx. every third data point will be 
            taken, value of 'delta' is only an advice for the routine to 
            calculate actual sampling step, which might differ in the end
        deltaz - sampling step in z-direction, for more see 'delta'
        interpolate - the order of interpolation ranging from 0 to 5,
            for 'interpolate' > 1 there might be problems with interpolation
            and NaNs are returned instead of numbers
        """
        
        if len(positions) != 2:
            print("create_zcut: Invalid positions.")
            return
        
        self.create_cut(spmdata, [positions], append=True, deltaz=deltaz, interpolate=interpolate)
        
    def create_cut(self, spmdata, positions, delta=1, deltaz=1, append=True, interpolate=0):
        """create 1D cut through spmdata
        
        spmdata - SPMdata instance through which the cut is done
        positions - sequence of 2-tuples or 3-tuples of point
            coordinates, these coordinates express indices of respective
            points in the spmdata array, but need not be integers, in which
            case interpolation is used, see 'interpolate', BEWARE: these are
            NOT sampling point, but only control points between which the
            points should be sampled, see 'delta' and 'deltaz'
            2-tuples case: (x, y) - x- and y-coordinate, only first
                len(spmdata) elements of 'positions' is used, the first
                2-tuple is assumed to lie in the first layer, the second
                2-tuple is assumed to lie in the second layer, etc.,
                in this case, sampling of z-coordinates is not performed
            3-tuples case: (z, x, y) - the first coordinate is the
                z-coordinate, no truncation is performed
        delta - sampling step in x- and y-direction, e.g. 'delta' = 3
            means that approx. every third data point will be 
            taken, value of 'delta' is only an advice for the routine to 
            calculate actual sampling step, which might differ in the end
        deltaz - sampling step in z-direction, for more see 'delta'
        append - option valid only for 2-tuple positions, if True, then
            the last position point is repeated so as the lenght of the cut
            is equal to the number of layers present in spmdata array
        interpolate - the order of interpolation ranging from 0 to 5,
            0 case: effectively no interpolation, take the nearest point
            1 case: linear interpolation
            etc.
        """
        
        # is 'positions' a mixture of n-tuples with different n?
        dims = set(map(len, positions))
        if len(dims) != 1:
            print("create_cut: Invalid positions.")
            return
            
        # if not, then treat 2-tuple and 3-tuple cases
        if len(positions[0]) == 2:
            # 2-tuple case
            numl = len(spmdata)
            posl = positions[:numl]
            if append and len(posl) < numl:
                posl += (numl - len(posl))*[posl[-1]]                
            posl = [[i, x, y] for i, (x, y) in enumerate(posl)]
            
        elif len(positions[0]) == 3:
            # 3-tuple case            
            posl = sorted(positions)
            
        else:
            print("create_cut: Invalid positions.")
            return

        # calculate sampling points
        poszl = [item[0] for item in posl]
        posxl = [item[1] for item in posl]
        posyl = [item[2] for item in posl]
        
        self.posx, self.posy, self.posz = self._process_cut_coords(poszl, posxl, posyl, delta, deltaz)

        # creation of the cut itself
        self.create_cutzxy_raw(spmdata, self.posz, self.posx, self.posy, interpolate=interpolate)
        
    def create_cutzxy_raw(self, spmdata, posz, posx, posy, interpolate=0):
        """create 1D cut through spmdata
        
        spmdata - SPMdata instance through which the cut is done
        posx - x-coordinate list of sampling points
        posy - y-coordinate list of sampling points
        posz - z-coordinate list of sampling points
        interpolate - the order of interpolation ranging from 0 to 5
        """
        
        minlen = min(len(posx), len(posy), len(posz))
        if minlen != max(len(posx), len(posy), len(posz)):
            print("create_cutzxy_raw: Unequal lengths of coordinate lists. "
                "Only first {} elements used.".format(minlen))
                
#        if minlen == 0:
#            print("create_cutzxy_raw: No positions specified.")
#            return

        self.minmax = {chan: (np.nanmin(val), np.nanmax(val)) for chan, val in spmdata.channels.items()}

        self.posx = posx[:minlen]
        self.posy = posy[:minlen]
        self.posz = posz[:minlen]

        for chan, val in spmdata.channels.items():
            # BACHA, BACHA, BACHA!!! i zde se projevuje prehozeni xove a yove souradnice
            self.channels[chan] = ndimage.map_coordinates(val, [self.posz, self.posy, self.posx], 
                order=interpolate, mode='constant', cval=np.nan)
#            self.channels[chan] = ndimage.map_coordinates(val, [self.posz, self.posx, self.posy], 
#                order=interpolate, mode='constant', cval=np.nan)
            # nearest-neighbour approximation
#            self.channels[chan] = val[self.posz, self.posy, self.posx]        

    def show_channels(self, *cl, tdata=None, linespec="b-", rescale=False):
        """plot 1D cuts of channels
        
        cl - list of channels to plot
        tdata - attribute of self which should be used as x-axis
        linespec - format specification of lines
        """
        
        # process tdata
        if tdata is not None:
            if hasattr(self, tdata):
                tdatal = getattr(self, tdata)
            else:
                print("show_channels: Invalid tdata.")
                return
        else:
            tdatal = None

        # pick only valid channels
        chanlist = [chan for chan in cl if chan in self.channels]
            
        if len(chanlist) == 0:
            print("show_channels: No valid channels.")
            return
            
        # create figure
        fig, subfiglist = plt.subplots(1, len(chanlist), squeeze=False, num="1D Cut - channel: " + ", ".join(chanlist))
        subfiglist = subfiglist[0]
        
        # plot each valid channel
        for i, chan in enumerate(chanlist):
            if tdatal is None:
                # no tdata was specified
                subfiglist[i].plot(self.channels[chan], linespec)
                subfiglist[i].set_xlabel("No. of steps")
            else:
                # tdata specified and valid
                ydatal = [y for x, y in sorted(zip(tdatal, self.channels[chan]))]
                subfiglist[i].plot(sorted(tdatal), ydatal, linespec)
                subfiglist[i].set_xlabel(tdata)
            
            if rescale: subfiglist[i].set_ylim(self.minmax[chan])    
            subfiglist[i].set_title("Channel: {}".format(chan))
            subfiglist[i].set_ylabel("{}".format(UNITS[chan]))
        
        # show
        fig.tight_layout()
        plt.show()

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
