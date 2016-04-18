import struct
import numpy as np
import matplotlib.pyplot as plt
import tkinter, os, math

class InvalidFormatError(Exception): pass

#CHAN_DICT = {"ELEV":"elev", "CURR":"curr", "FREQ":"freq",
#	"EXCI":"exci", "OMEG":"omeg", "AMPL":"ampl", "PHAS":"phas"}
#
#for key, val in CHAN_DICT.items():
#    globals()[key] = val
#    globals()[key + 'B'] = val + 'b'

## multiples
#NM = 10e-9  # nano

UNKNOWN = "unknown"

# forward scan channels
ELEV = "elev"   # Z
CURR = "curr"   # Current
FREQ = "freq"   # Frequency Shift
EXCI = "exci"   # Excitation
OMEG = "omeg"   # LIY_1_omega
AMPL = "ampl"   # Amplitude
PHAS = "phas"   # Phase
ZREL = "zrel"   # Z rel
CHAN_IDENT = (ELEV, CURR, FREQ, EXCI, OMEG, AMPL, PHAS, ZREL)

# backward scan channels
ELEVB = "elevb"   # Z backward
CURRB = "currb"   # Current backward
FREQB = "freqb"   # Frequency Shift backward
EXCIB = "excib"   # Excitation backward
OMEGB = "omegb"   # LIY_1_omega backward
AMPLB = "amplb"   # Amplitude backward
PHASB = "phasb"   # Phase backward
ZRELB = "zrelb"   # Z rel backward
CHAN_IDENT += (ELEVB, CURRB, FREQB, EXCIB, OMEGB, AMPLB, PHASB, ZRELB)

def is_backward_scan(chan):
    """decides whether string 'chan' denotes backward scan
    """
    
    return chan.endswith('b')
    
def get_backward_scan(chan):
    """returns name of the backward-scan version of channel 'chan'
    """
    
    return chan + 'b'

def get_forward_scan(chan):
    """returns name of the forward-scan version of channel 'chan',
    assumes that chan is a backward-scan
    """
    
    return chan[:-1]

# file format suffixes
SUFFIXES = ("sxm",)

# units
UNITS = {ELEV:"m", CURR:"A", FREQ:"Hz", EXCI:"V", OMEG:"V?", AMPL:"m", PHAS:"deg", ZREL:"m",
ELEVB:"m", CURRB:"A", FREQB:"Hz", EXCIB:"V", OMEGB:"V?", AMPLB:"m", PHASB:"deg", ZRELB:"m",
UNKNOWN:"?"
}


def find_offind(indi, array, tolerancy):
    """finds appropriate offind by skipping initial NaNs
    
    indi - initial offind
    array - array with respect to which to find offind
    tolerancy - how big difference between initial and
        final offind is allowed
    """

    flag = False
    for line in array:
        for i, item in enumerate(line):
            if not np.isnan(item):
                res = i
                flag = True
                break
        if flag: break
    else:
        print("find_offind: Whole array full of NaNs.")
        res = -1

    if abs(res - indi) > tolerancy:
        print("find_offind: Too big difference between offinds: final - initial = {}.".format(res - indi))
#        res = indi

    return res

def get_offinds_corr(arr1, arr2):            
    """return argmax indices of cross-correlation array, which is
    calculated by fft
    """
    corrarr = np.fft.ifftn(np.fft.fftn(arr1)*np.conj(np.fft.fftn(arr2)))        
    corrarr = np.roll(corrarr, int(np.floor(corrarr.shape[0] / 2)), axis=0)
    corrarr = np.roll(corrarr, int(np.floor(corrarr.shape[1] / 2)), axis=1)
    relx, rely = np.unravel_index(np.argmax(corrarr), corrarr.shape)

    return relx, rely

def count_nans(array):
    """counts how many NaNs are present in 'array'
    """
    
    k = 0
    for item in array:
        if np.isnan(item): k += 1
    print("count_nans: len, nans, len - nans: ", len(array), k, len(array) - k) 

def vector_handedness(point1, point2, point3):
    """ if vectors (point1 - point2) and (point3 - point2)
    form a left-handed basis with vector (0, 0, -1), return
    the right-handed counterpart of (point1, point2, point3),
    otherwise return (point1, point2, point3)
    """
    
    vec12 = point1 - point2
    vec32 = point3 - point2
    crossvec = np.cross(vec12, vec32)

    rvec = np.array([0, 0, 1])
    scalprod = np.dot(crossvec, rvec)      
#    print("scalprod: ", scalprod)
              
    if abs(scalprod) < 1e-9:
        rvec = np.array([0, 1, 0])
        scalprod = np.dot(crossvec, rvec)
        if scalprod < 0:        
            print("vector_handedness: Control points formed a left-hand basis,"
                " so they were reordered.")
            point3, point1 = point1, point3
            crossvec *= -1    
    elif scalprod < 0:
        print("vector_handedness: Control points formed a left-hand basis,"
            " so they were reordered.")
        point3, point1 = point1, point3
        crossvec *= -1    
    
    return point1, point2, point3, crossvec


