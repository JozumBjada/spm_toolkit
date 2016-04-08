from spm_data_load import *

class InvalidFormatError(Exception): pass

class SPMdataLayer:
    """single layer of SPM data, i.e. single scan/session
    """

    def __init__(self, filename=None):

        self.xnum, self.ynum = 0, 0
        self.height = None

        # self.channels = {"nameOfTheChannel, e.g. ELEV":pointerToChan, ...}
        self.channels = {}

#        self.elev, self.curr, self.freq, self.exci, self.omeg = None, None, None, None, None
#        self.ampl, self.phas = None, None

        self.xoffset, self.yoffset = None, None
        self.xoffind, self.yoffind = {}, {}
        self.xran, self.yran = None, None

        self.filename = ""
        self.metadata = {}

        if filename is not None:
            self.load_data(filename)

    def copy(self, ind, cl):
        """return a copy of a layer
        
        ind - index of the new layer in new SPMdata structure
        cl - an external dictionary of channels to be copied
        """
        
        copylay = SPMdataLayer()
        itemlist = ["xnum", "ynum", "xran", "yran", "height",
            "xoffset", "yoffset", "filename", "metadata"]
        # more universal would be to use itemlist=locals() or
        # something like that, but the explicit form used above
        # makes it clear what is actually copied
        
        for item in itemlist:
            setattr(copylay, item, getattr(self, item))
        
        copylay.xoffind = {chan: val for chan, val in self.xoffind.items() if chan in cl}
        copylay.yoffind = {chan: val for chan, val in self.yoffind.items() if chan in cl}
        
        xnum, ynum = copylay.xnum, copylay.ynum
        for chan in cl:
            if chan not in self.channels:
#                print("copy: Channel '{}' not present in the layer. Skipped.".format(chan))
                continue
            xoff, yoff = copylay.xoffind[chan], copylay.yoffind[chan]
            copylay.setchan(chan, cl[chan][ind, xoff:xoff + xnum, yoff:yoff + ynum])
        
        return copylay

    def delchan(self, chan):
        """delete channel from the structure
        """
        
        if chan not in self.channels:
            return
            
        del self.channels[chan]
#        del self.shapes[chan]
        del self.xoffind[chan]
        del self.yoffind[chan]
        
        delattr(self, chan)

    def setchan(self, chan, val):
        """set channel to the structure
        """
        
#        if chan not in self.chan_ident:
#            print("setchan: Invalid channel '{}' for given layer.".format(chan))
#            return
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

    def _load_data_nanonis(self, filename):
        """load data from nanonis format file
        """

        with open(filename, "rb") as file:

            # read the heading

            try:
                while True:
                    linepar = file.readline().strip()
                    lineval = file.readline().strip()

                    if linepar == b":SCANIT_TYPE:":
                        DATA_FILE_TEMP, ENDIAN_TEMP = lineval.split()
#                        print(DATA_FILE_TEMP, ENDIAN_TEMP)
                        if DATA_FILE_TEMP != b"FLOAT" or ENDIAN_TEMP != b"MSBFIRST":
                            raise InvalidFormatError
                    elif linepar == b":Z-CONTROLLER:":
#                        datasetZNames = lineval.strip()
                        # print("datasetZNames {}".format(datasetZNames))
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
                        self.xnum, self.ynum = int(lineval[0]), int(lineval[1])
                    elif linepar == b":SCAN_DIR:":
                        incr = 1   # !!! KDYZTAK incr = -1
                        if lineval == b"down":
                            incr *= -1    
                    # elif linepar == b":Scan>pixels/line:":
                    #     self.xnum = int(lineval)
                    # elif linepar == b":Scan>lines:":
                    #     self.ynum = int(lineval)
                    elif linepar == b":SCAN_OFFSET:":
                        self.xoffset, self.yoffset = lineval.split()
                        self.xoffset, self.yoffset = float(self.xoffset), float(self.yoffset)
                    elif linepar == b":SCAN_RANGE:":
                        self.xran, self.yran = lineval.split()
                        self.xran, self.yran = float(self.xran), float(self.yran)
                    elif linepar == b":DATA_INFO:":
#                        datasetNames = lineval.strip()
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
#                        print(dataval)
                    elif linepar == b":SCANIT_END:":
#                        print(lineval)
                        break
                    else:
                        pass
            except InvalidFormatError:
                    print("_load_data_nanonis: Invalid data format.")
                    return

            # read the binary data

            # NUM_BYTES_PER_NANONIS_FLOAT = 4
            ENDIAN = '>'
            DATA_TYPE = 'f'

            space = file.read(2)
            while space != b'\x1a\x04':
                space = space[1:] + file.read(1)
                
            datalist = file.read()

            backwardscanning = [(x[3] == b"both") for x in dataval]
            numchandir = len(dataval) + sum(backwardscanning)

            datalist = struct.unpack(ENDIAN + DATA_TYPE*self.xnum*self.ynum*numchandir, datalist)

            data = np.array(datalist)
            data.resize((numchandir, self.xnum, self.ynum))

            chanstrdict = {
                "Z":ELEV, "Current":CURR, "Frequency_Shift":FREQ, "Excitation":EXCI,
                "LIY_1_omega":OMEG, "Amplitude":AMPL, "Phase":PHAS
                }

            i = 0
            for item in dataval:
                itm = item[1].decode("utf-8")
                if itm in chanstrdict:
                    if item[3] == b"both":
                        # print(chanstrdict[itm])
                        setattr(self, chanstrdict[itm], data[i][::incr].copy())
                        setattr(self, get_backward_scan(chanstrdict[itm]), data[i + 1][::incr,::-1].copy())
                        i += 2
                    else:
                        setattr(self, chanstrdict[itm], data[i][::incr].copy())
                        i += 1
                    self.channels.setdefault(chanstrdict[itm], getattr(self, chanstrdict[itm]))
                    self.channels.setdefault(get_backward_scan(chanstrdict[itm]), getattr(self, get_backward_scan(chanstrdict[itm])))
                else:
                    print("_load_data_nanonis: {} not known.".format(itm))
                    
#            self.shapes = {chan:(self.xnum, self.ynum) for chan in self.channels}
            self.xoffind = {chan:0 for chan in self.channels}
            self.yoffind = {chan:0 for chan in self.channels}



