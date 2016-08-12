import sys
sys.path[0] = '..'

from code.spmdata.spm_data_aux import *
from code.registration.registration_class import *
from code.registration.registration_compare import *

IMPORT_FILE = "compdata.hdr"

if __name__ == "__main__":

    # THIS SCRIPT IMPORT PREPARED CompareRegBlocks STRUCTURE AND CONTROL POINTS AND USE THESE TO MAKE A PRELIMINARY ROUGH REGISTRATION

    # controlpoints may be obtained with help of routine 'inspect_blocks'

#    controlpoints = [
#        [np.array([   0.        ,   79.1622945 ,  188.77805475]), np.array([   0.        ,  164.62990317,   76.08742999])],
#        [np.array([   0.        ,   70.45991106,  116.60540743]), np.array([   0.        ,  148.33039895,  198.90754911])]
#        ]
        
    controlpoints = [
        [np.array([   0.        ,  164.52339757,   72.36804463]), np.array([   0.        ,  144.72063566,  218.4134137 ])],
        [np.array([   0.        ,  141.78375426,  185.61508929]), np.array([   0.        ,   26.68020068,  151.57909226])]
        ]
        
    # import CompareRegBlocks structure
    comp = CompareRegBlocks()
    comp.import_data(IMPORT_FILE)
    comp.reg2.set_points(comp.reg1.ptsorig, comp.reg1.oriorig)

    # make a rough registration
    comp.find_parameters(controlpoints)
    
    # let us see how this rough registration works
    comp.inspect_blocks(cutlist=controlpoints)

    # now we let the automatic procedure to find the best match for registrations
    # first of all we say the procedure in which range of parameters to search for the best match
    zoomstep = 0.01
    rotatestep = 1
    skewstep = 0.5
    xshiftstep = 1
    yshiftstep = 1

#    pardict = {
#        "zoom":  (0, zoomstep),
#        "rotate":(2, rotatestep),
#        "xskew": (skewstep),
#        "yskew": (skewstep),
#        "xshift":(3, xshiftstep),
#        "yshift":(3, yshiftstep)
#        }
#    
    pardict = {}
    
    # find the best match
    comp.find_match(pardict, maxsteps=20, showoverlap=True, facfreq=4)

    # which values of parameters were needed to transform one block onto the other?
    comp.reg1.print_pars()

    # let us see how the best-matching result looks like
    comp.inspect_blocks()

    # let us see the overlap of both blocks of data
#    comp.show_overlap()

