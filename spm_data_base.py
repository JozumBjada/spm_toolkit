import struct, pickle
import numpy as np
import matplotlib.pyplot as plt
import tkinter, os, math
from scipy.ndimage import map_coordinates

class InvalidFormatError(Exception): pass

# file format suffixes
SUFFIXES = ("sxm",)

# suffix for channel names that have backward direction
BCKSUFF = 'b'
BCKSUFF_VAR = 'B'

## multiples
#NM = 10e-9  # nano

# DYNAMICAL DEFINITION OF CHANNEL NAME CONSTANTS AND 
# APPROPRIATE UNITS:

# (forward) scan channels and units
CHAN_DICT = {
    "ELEV":("elev", "m"  ), # Z
    "CURR":("curr", "A"  ), # Current
    "FREQ":("freq", "Hz" ), # Frequency Shift
	"EXCI":("exci", "V"  ), # Excitaion
	"OMEG":("omeg", "V?" ), # LIY_1_omega
	"AMPL":("ampl", "m"  ), # Amplitude
	"PHAS":("phas", "deg"), # Phase
    "UNKNOWN":("unknown", "?") # Unknown channel
	}

UNITS = {}

# create forward and backward scan channels and units
for key, (val, unit) in CHAN_DICT.items():
    globals()[key] = val
    globals()[key + BCKSUFF_VAR] = val + BCKSUFF
    UNITS[val] = unit
    UNITS[val + BCKSUFF] = unit

# STATIC DEFINITION OF CHANNEL NAME CONSTANTS AND APPROPRIATE
# UNITS:

## forward scan channels
#ELEV = "elev"   # Z
#CURR = "curr"   # Current
#FREQ = "freq"   # Frequency Shift
#EXCI = "exci"   # Excitation
#OMEG = "omeg"   # LIY_1_omega
#AMPL = "ampl"   # Amplitude
#PHAS = "phas"   # Phase
#ZREL = "zrel"   # Z rel

UNKNOWN = "unknown"

## backward scan channels
#ELEVB = "elevb"   # Z backward
#CURRB = "currb"   # Current backward
#FREQB = "freqb"   # Frequency Shift backward
#EXCIB = "excib"   # Excitation backward
#OMEGB = "omegb"   # LIY_1_omega backward
#AMPLB = "amplb"   # Amplitude backward
#PHASB = "phasb"   # Phase backward
#ZRELB = "zrelb"   # Z rel backward

# units
UNITS = {ELEV:"m", CURR:"A", FREQ:"Hz", EXCI:"V", OMEG:"V?",
    AMPL:"m", PHAS:"deg", ELEVB:"m", CURRB:"A", FREQB:"Hz", 
    EXCIB:"V", OMEGB:"V?", AMPLB:"m", PHASB:"deg",
    UNKNOWN:"?"
}

# FORWARD AND BACKWARD CHANNEL INQUIRY ROUTINES

def is_backward_scan(chan):
    """decides whether string 'chan' denotes backward scan
    """
    
    return chan.endswith(BCKSUFF)
    
def get_backward_scan(chan):
    """returns name of the backward-scan version of channel 'chan'
    """
    
    return chan + BCKSUFF

def get_forward_scan(chan):
    """returns name of the forward-scan version of channel 'chan',
    assumes that chan is a backward-scan
    """
    
    return chan[:-1]

# RESAMPLING

def gen_resample_data_2d(inparr, xno, yno, xn, yn, order=0):
    """perform upsampling in 2D
    """
   
    x = np.linspace(0, xno - 1, xn)
    y = np.linspace(0, yno - 1, yn)
    coords = np.meshgrid(x, y, indexing='ij')
    arr = map_coordinates(inparr,
        coords, order=order, cval=np.nan)
        
    return arr

def gen_resample_data_3d(inparr, xno, yno, zno, xn, yn, zn,
    order=0):
    """perform upsampling in 3D
    """
   
    x = np.linspace(0, xno - 1, xn)
    y = np.linspace(0, yno - 1, yn)
    z = np.linspace(0, zno - 1, zn)
    coords = np.meshgrid(z, x, y, indexing='ij')
    arr = map_coordinates(inparr,
        coords, order=order, cval=np.nan)
        
    return arr

# PICKLING

def pickle_data(data, filename):
    """export data to external file by mechanism provided by
    pickle
    
    data - data to export
    filename - where to export the data
    """
    
    with open(filename, 'wb') as f:
        pickle.dump(obj, f)

def unpickle_data(datalist, filename):
    """import data from external file by mechanism provided by
    pickle
    
    datalist - list whose first element will be assigned data
        from file
    filename - where to import data from
    """

    with open(filename, 'rb') as f:
        spmdatalist[0] = pickle.load(f)

