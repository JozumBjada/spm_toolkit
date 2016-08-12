from .spm_data_base import *

class InvalidFormatError(Exception): pass

class SPMdataLayer:
    """single layer of SPM data, i.e. single scan/session
    """

    def __init__(self, filename=None):

        # number of data points in x- and y-direction
        self.xnum, self.ynum = 0, 0
        
        # height of the probe 
        self.height = None
        
        # resolution of the data,
        # i.e. (data range) / (number of data points)
        self.xstep, self.ystep = 0, 0

        # self.channels = {
        # "nameOfTheChannel, e.g. ELEV":pointerToChannel, ...}
        self.channels = {}

        # offsets in real physical dimensions
        self.xoffset, self.yoffset = None, None
        
        # corresponding offsets in terms of number of data points
        self.xoffind, self.yoffind = {}, {}
        
        # ranges of probed surface field
        self.xran, self.yran = None, None

        self.filename = ""
        self.scanangle = 0
        self.metadata = {}
        self.units = {}

        if filename is not None:
            self.load_data(filename)

    def copy(self, ind, cl):
        """return a copy of a layer
        
        ind - index of the new layer in new SPMdata structure
        cl - an external dictionary of channels to be copied
        """
        
        copylay = SPMdataLayer()
        itemlist = ["xnum", "ynum", "xran", "yran", "height",
            "xoffset", "yoffset", "filename", "metadata",
            "xstep", "ystep"]
        # more universal would be to use itemlist=locals() or
        # something like that, but the explicit form used above
        # makes it clear what is actually copied
        
        for item in itemlist:
            setattr(copylay, item, getattr(self, item))
        
        copylay.xoffind = {chan: val
            for chan, val in self.xoffind.items() if chan in cl}
        copylay.yoffind = {chan: val
            for chan, val in self.yoffind.items() if chan in cl}
        
        xnum, ynum = copylay.xnum, copylay.ynum
        for chan in cl:
            if chan not in self.channels:
                continue
            xoff = copylay.xoffind[chan]
            yoff = copylay.yoffind[chan]
            copylay.set_chan(chan, cl[chan][ind,
                xoff:xoff + xnum, yoff:yoff + ynum])
        
        return copylay

    def del_chan(self, chan):
        """delete channel from the structure
        """
        
        if chan not in self.channels:
            return

        del self.channels[chan]
        del self.xoffind[chan]
        del self.yoffind[chan]            
        
        if hasattr(self, chan): delattr(self, chan)

    def set_chan(self, chan, val):
        """set channel to the structure
        """
        
        if chan.isidentifier():
            setattr(self, chan, val)
        self.channels[chan] = val

    def load_data(self, filename):
        """load data from external file

        filename...name of the file, valid filename is assumed
        """

        self.filename = filename

        if filename.endswith(".sxm"):
            self._load_data_nanonis(filename)
#        elif filename.endswith("another format..."):
#           self._load...
#        etc.
        else:
            print("load_data: Unknown data format.")
            
        # data resolution
        self.xstep = self.xran / self.xnum
        self.ystep = self.yran / self.ynum

    @property
    def time(self):
        """return string containing the time and data of the
        experiment; this attribute may be used e.g. to sort layers
        according to acquisition time
        """
        
        return self.metadata["time"] + "," + self.metadata["date"]

#    @property
#    def scanangle(self):
#        """return scanning angle
#        """
#        
#        return self.metadata["angle"]

    def _load_data_nanonis(self, filename):
        """load data from nanonis format file
        """

        with open(filename, "rb") as file:

            # READ HEADING
            while True:
                linepar = file.readline().strip()
                lineval = file.readline().strip()

                if linepar == b":SCANIT_TYPE:":
                    DATA_FILE_TEMP, ENDIAN_TEMP = lineval.split()
#                        print(DATA_FILE_TEMP, ENDIAN_TEMP)
                    if DATA_FILE_TEMP != b"FLOAT" or \
                        ENDIAN_TEMP != b"MSBFIRST":
                        print("_load_data_nanonis: Invalid data"
                            " format.")
                        return

                elif linepar == b":Z-CONTROLLER:":
                    dataz = []
                    while lineval:
                        filepos = file.tell()
                        lineval = file.readline()
                        if lineval.startswith(b":"):
                            file.seek(filepos)
                            break
                        else:
                            dataz.append(lineval.split())
#                        print(dataz)
                elif linepar == b":REC_DATE:":
                    self.metadata["date"] = lineval.decode("utf-8")
                elif linepar == b":REC_TIME:":
                    self.metadata["time"] = lineval.decode("utf-8")
                elif linepar == b":Z-Controller>Z (m):":
                    self.height = float(lineval.decode("utf-8"))
                elif linepar == b":COMMENT:":
                    while lineval:
                        filepos = file.tell()
                        lineval = file.readline()
                        if lineval.startswith(b":"):
                            file.seek(filepos)
                            break
                elif linepar == b":SCAN_PIXELS:":
                    lineval = lineval.split()
                    self.xnum = int(lineval[0])
                    self.ynum = int(lineval[1])
                elif linepar == b":SCAN_DIR:":
                    incr = 1   # !!! KDYZTAK incr = -1
                    if lineval == b"down":
                        incr *= -1    
                elif linepar == b":SCAN_ANGLE:":
                    self.scanangle = float(
                        lineval.decode("utf-8"))
                # elif linepar == b":Scan>pixels/line:":
                #     self.xnum = int(lineval)
                # elif linepar == b":Scan>lines:":
                #     self.ynum = int(lineval)
                elif linepar == b":SCAN_OFFSET:":
                    self.xoffset, self.yoffset = lineval.split()
                    self.xoffset = float(self.xoffset)
                    self.yoffset = float(self.yoffset)
                elif linepar == b":SCAN_RANGE:":
                    self.xran, self.yran = lineval.split()
                    self.xran = float(self.xran)
                    self.yran = float(self.yran)
                elif linepar == b":DATA_INFO:":
                    dataval = []
                    while True:
                        filepos = file.tell()
                        lineval = file.readline()
                        if lineval.startswith(b":"):
                            file.seek(filepos)
                            dataval.pop()
                            break
                        else:
                            dataval.append(lineval.split())
                elif linepar == b":SCANIT_END:":
                    break
                else:
                    pass

            # READ BINARY DATA
            
            # format specification
            ENDIAN = '>'
            DATA_TYPE = 'f'
            bckwardscanning = [(x[3] == b"both") for x in dataval]
            numchandir = len(dataval) + sum(bckwardscanning)
            spec = ENDIAN + str(self.xnum*self.ynum*numchandir)
            spec += DATA_TYPE

            # move to the beginning of valid data
            space = file.read(2)
            while space != b'\x1a\x04':
                space = space[1:] + file.read(1)
                
            # read and process valid data
            datalist = file.read()
            datalist = struct.unpack(spec, datalist)
            data = np.array(datalist)
            data.resize((numchandir, self.xnum, self.ynum))

            # read forward and backward scans
            i = 0
            for item in dataval:
                # name of forward scan
                channame = item[1].decode("utf-8")

                # name of backward scan
                channamebck = get_backward_scan(channame) 

                # set units
                self.units[channame] = item[2].decode("utf-8")
                self.units[channamebck] = item[2].decode("utf-8")

                # load data                                
                if item[3] == b"both":
                    self.channels.setdefault(channame,
                        data[i][::incr].copy())
                    self.channels.setdefault(channamebck,
                        data[i + 1][::incr,::-1].copy())
                    i += 2
                else:
                    self.channels.setdefault(channame,
                        data[i][::incr].copy())
                    self.channels.setdefault(channamebck, None)
                    i += 1
                                        
            # initilize offinds
            self.xoffind = {chan:0 for chan in self.channels}
            self.yoffind = {chan:0 for chan in self.channels}

