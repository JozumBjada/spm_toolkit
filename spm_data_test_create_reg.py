import sys
sys.path[0] = '..'

from code.spmdata.spm_data_aux import *
from code.registration.registration_class import *
from code.registration.registration_compare import *

FIRST_SET_OF_DATA  = "/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration/data_uncharged"
SECOND_SET_OF_DATA = "/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration/data_charged/shrange"
EXPORT_FILE = "compdata"

if __name__ == "__main__":

    # THIS SCRIPT LOADS DATA FOR THE TWO BLOCKS TO UNDERGO REGISTRATION AND EXPORTS PREPARED CompareRegBlocks STRUCTURE TO EXTERNAL FILE 

    # retrieve the first set of data by creating SPMdata structure from files stored in directory FIRST_SET_OF_DATA
    spm = create_spm(paths=FIRST_SET_OF_DATA)
    
    # we would like to use frequency shift data, so let us see which channels are present in spm
    spm.print_channels()
    
    # we see that 'freq' is present, how does the data look like?
    spm.inspect_channels(FREQ)
    
    # we do the same thing with the second set of data, but here we have to sort it first according to file names
    spmc = create_spm(paths=SECOND_SET_OF_DATA)
    spmc.sort_spm(key="filename")
    spmc.inspect_channels(FREQ)
    
    # both sets of data have too many layers, let us focus on first 'numlay' layers
    numlay = 2
    spm.delete_layers(range(numlay), leave=True)  
    spmc.delete_layers(range(numlay), leave=True)  

    # the data are not aligned in context of layers, so let us align them according to channel 'curr'
    fac = 2
    spm.align_correlate(rchan="Current", fac=fac)    
    spmc.align_correlate(rchan="Current", fac=fac)    

    # OK, now we have our data prepared, so put them into registration blocks    
    order = 0
    reg1 = RegBlock(spm.Frequency_Shift, order=order)
    reg2 = RegBlock(spmc.Frequency_Shift, order=order)

    # we want to do registration, so we put both registration blocks into the class CompareRegBlocks, which is responsible for registration process
    comp = CompareRegBlocks(reg1, reg2)
    
    # we get rid of redundant NaNs
    comp.strip_nans()
    
    # since we want to reuse our data, we export our 'comp' object containing reg1 and reg2
    comp.export_data(EXPORT_FILE)
    
    # it is good to help correlation routine to find the first rough guess of mutual position of both blocks, we can do so by feeding control point to the routine 'find_parameters'
    # here we set control points via graphical interface and print these point on the screen;
    # inspect_blocks launches a window with two blocks shown in an interactive graphical mode, control points are then set as follows:
    #   1) by left-mouse button choose two points on the left-hand-side block of data
    #   2) press CTRL + P and then CTRL + D
    #   3) by left-mouse button choose two points on the right-hand-side block of data; 
    #       these two points should represent the same feature on the second block as is represented by points in 1) in the first block
    #   4) press CTRL + P and then CTRL + D
    #   5) close the window 
    controlpoints = comp.inspect_blocks()
    print("controlpoints: {}".format(controlpoints))   

