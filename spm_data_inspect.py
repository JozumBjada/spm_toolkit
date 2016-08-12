from .spm_data_base import *
from mpl_toolkits.axes_grid1 import make_axes_locatable    
import matplotlib.ticker as tck

"""this module implements procedure "inspect_channels",
all the other routines in this file are auxiliary routines
for the main procedure "inspect_channels"
"""

# initialize graphics
def _initialize_plots(arrs, unitdir, vextr, layer, interpolation,
    xind, yind, linestyle, linecolor, pointcolor, size, aspect,
    cuts, numl, rightlim, bottolim):
    """initialize main plot, cut plots, horizontal and vertical
    lines and text
    
    arrs - list of channels, i.e.
        #arrs={(i, id(spmdata), chan):array, ...}
        arrs = [(id(spmdata), chan, array), ...]
    unitdir - dictionary of corresponding units, i.e.
        unitdir={(i, id(spmdata), chan):units[chan]}
    vextr - dictionary of (min, max) values
    layer - initial layer to be drawn
    interpolation - interpolation order for imshow
    xind, yind - initial position of the cursor
    linestyle - style of lines in the images
    linecolor - color of lines in the images
    pointcolor - color of the point marking the z-cut
    size - size of cut images
    aspect - aspect of cut images
    cuts - data for cuts
    numl - number of layers
    rightlim, bottolim - real dimensions of the arrays;
        if None, then no scaling is used
    """

    # create figure and layout
    chanlist = [chan for (idnum, chan, val) in arrs]
    fig, subfiglist = plt.subplots(1, len(arrs),
        num="Channel: " + ", ".join(chanlist),
        squeeze=True)
    if len(arrs) == 1: subfiglist = [subfiglist]

    # create auxiliary lists
    axe1list, axe2list, axe3list = [], [], []
    poi1list, poi2list, poi3list = [], [], []        
    img1list, img2list, img3list = [], [], []        
    ver1list, ver2list, ver3list = [], [], []
    hor1list, hor2list, hor3list = [], [], []
    zcu1list, zcu2list, zcu3list = [], [], []
    tex1list, tex2list = [], []

    iddict, k = {}, 0
    numxticks, numyticks = 6, 10

    # create subfigures, all created in a single row
    for i, (idnum, chan, val) in enumerate(arrs):
        key = (idnum, chan)
        vmin, vmax = vextr[key]

        x, y = xind, yind

        # identification number of given SPMdata
        # structure
        if idnum not in iddict:
            iddict[idnum] = k
            k += 1
        idnum = iddict[idnum]
            
        # main image
        axes1 = subfiglist[i]
        divider = make_axes_locatable(axes1)
        axes1.set_adjustable('box-forced')
        extent = (0, val[layer].shape[1], val[layer].shape[0], 0)
#        extent = (0, val[layer].shape[1], 0, val[layer].shape[0])
        imag1 = axes1.imshow(val[layer], extent=extent,
            vmin=vmin, vmax=vmax, interpolation=interpolation,
            origin='lower')            
        vlin1 = axes1.axvline(x=x, ls=linestyle, c=linecolor)
        hlin1 = axes1.axhline(y=y, ls=linestyle, c=linecolor)
        poin1, = axes1.plot([], [], ls=linestyle, lw=2,
            c=linecolor)
        zcut1, = axes1.plot([], [], ls='', marker='o', lw=2, 
            c=pointcolor)
        plt.setp(axes1.get_xticklabels(), visible=False)

        # physical dimension labels
        if rightlim is not None and bottolim is not None:
            xticks = np.linspace(
                0, val[layer].shape[0] - 1, numxticks)
            yticks = np.linspace(
                0, val[layer].shape[1] - 1, numyticks)
            xticklabels = np.linspace(0, rightlim[key], numxticks)
            yticklabels = np.linspace(0, bottolim[key], numyticks)
            yticklabels = ["{:.4}".format(i) for i in yticklabels]
            
            axes1.set_xticks(xticks)
            axes1.set_yticks(yticks)
            axes1.set_xticklabels(xticklabels)
            axes1.set_yticklabels(yticklabels)
       
        # side cut
        axes2 = divider.append_axes("right", size=size, pad=0.2,
            sharey=axes1)
        axes2.set_adjustable('box-forced')
        hlin2 = axes2.axhline(y=y, ls=linestyle, c=linecolor)
        vlin2 = axes2.axvline(x=layer, ls=linestyle, c=linecolor)
        poin2, = axes2.plot([], [], ls=linestyle, lw=2,
            c=linecolor)
        zcut2 = axes2.axhline(y=-1, ls=linestyle, c=pointcolor)
        extent = (0, numl, val[layer].shape[0], 0)
        imag2 = axes2.imshow(cuts[key][0], vmin=vmin, vmax=vmax,
            interpolation=interpolation, aspect=aspect,
            extent=extent)
            
        # for numl == 1 there are errors, since Python has 
        # problems with xticks = [0, 0, 0], that is why set(...)
        # is used below
        ticks = list(set([0, numl // 2, numl - 1]))
        axes2.xaxis.set_ticks(ticks)
        plt.setp(axes2.get_yticklabels(), visible=False)

        # colorbar
        # if imag1 is full of the same value (i.e. the
        # range of values has zero length) the plt.colorbar issues
        # 'RuntimeWarning: invalid value encountered in less', so
        # we suppress it by 'warning' package tools
        axes5 = divider.append_axes("right", size="15%", pad=0.2)
        axes5.set_adjustable('box-forced')

        import warnings        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            plt.colorbar(imag1, cax=axes5)
        
        # bottom cut
        axes3 = divider.append_axes("bottom", size=size, pad=0.2,
            sharex=axes1)
        axes3.set_adjustable('box-forced')
        vlin3 = axes3.axvline(x=x, ls=linestyle, c=linecolor)
        hlin3 = axes3.axhline(y=layer, ls=linestyle, c=linecolor)
        poin3, = axes3.plot([], [], ls=linestyle, lw=2,
            c=linecolor)
        zcut3 = axes3.axvline(x=-1, ls=linestyle, c=pointcolor)
        extent = (0, val[layer].shape[1], numl, 0)
        imag3 = axes3.imshow(cuts[key][1], vmin=vmin, vmax=vmax,
            interpolation=interpolation, aspect=aspect,
            extent=extent)

        # for numl == 1 there are errors, since Python has 
        # problems with xticks = [0, 0, 0], that is why set(...)
        # is used below
        ticks = list(set([0, numl // 2, numl - 1]))
        axes3.yaxis.set_ticks(ticks)
        
        # metadata
        axes4 = divider.append_axes("top", size=0.9)
        axes4.set_axis_off()

        text0 = axes4.text(0, .9,
            "SPM: {} \nChannel: {}".format(idnum, chan))
        text1 = axes4.text(0, .7,
            "layer no. {}".format(layer % numl))
        text2 = axes4.text(0, .5,
            _str_value_at_cursor(
            val[layer], unitdir[key], xind, yind))

        # update lists        
        ver1list.append(vlin1); ver2list.append(vlin2);
        ver3list.append(vlin3)
        hor1list.append(hlin1); hor2list.append(hlin2);
        hor3list.append(hlin3)

        axe1list.append(axes1); axe2list.append(axes2);
        axe3list.append(axes3)
        img1list.append(imag1); img2list.append(imag2);
        img3list.append(imag3)
        
        poi1list.append(poin1); poi2list.append(poin2);
        poi3list.append(poin3)
        zcu1list.append(zcut1); zcu2list.append(zcut2);
        zcu3list.append(zcut3)
        
        tex1list.append(text1); tex2list.append(text2)

    # gather lists
    axes = [axe1list, axe2list, axe3list]
    lins = [ver1list, ver2list, ver3list,
            hor1list, hor2list, hor3list]
    imgs = [img1list, img2list, img3list]
    texs = [tex1list, tex2list]
    pois = [poi1list, poi2list, poi3list]
    zcus = [zcu1list, zcu2list, zcu3list]
    
    return fig, axes, lins, imgs, texs, pois, zcus

# string pertaining to the value at cursor position
def _str_value_at_cursor(chanarr, chanunit, xind, yind):

    if xind is None or yind is None:
        val = 0
    else:
        # TODO: MOZNA PREHOZENE xind A yind???
#        val = chanarr[xind, yind]

        # for origin='lower' in imshow we have to reverse it [???]:
#        val = chanarr[xind, -yind - 1]
        val = chanarr[-xind - 1, yind]
        
    return "value: {:.5g} {}".format(val, chanunit)

# function for creation z-direction cuts
def _create_cuts(chanarr, ind1, ind2):
    """create cut1 and cut2 from array chanarr corresponding to
    channel 'chan' in SPMData structure
    """
    
#    ind1, ind2 = ind2, ind1 # ???
    cut1 = chanarr[:, -ind1 - 1]
    cut2 = chanarr[:, :, ind2]
    
    return cut2.transpose()[::-1], cut1 # !!!

# event handler for choosing a point through which z-cut
# is to be done
class _CutZedPoint:
    def __init__(self, fig, axes1, linelists, cutzlist):
        self.linelist, self.linelist2, self.linelist3 = linelists
        self.axes1, self.cutzlist, self.fig = axes1, cutzlist, fig
        self.xpos, self.ypos = [], []
                        
        self.cid1 = fig.canvas.mpl_connect(
            'button_press_event', self)
        self.cid2 = fig.canvas.mpl_connect(
            'key_press_event', self)

    def __call__(self, event):
        if event.name == 'button_press_event':
            # choose control point by left mouse button
            if event.inaxes in self.axes1 and \
                event.button == 3 and \
                event.xdata is not None and \
                event.ydata is not None:
                self.xpos.append(event.xdata) # event.xdata ???
                self.ypos.append(event.ydata) # event.ydata ???
                    
                # even though linelist is represented by a plot,
                # meaningful are only individual points in a plot,
                # they specify x- and y-coordinate of the z-cut    
                for line in self.linelist:
                    line.set_xdata(self.xpos)
                    line.set_ydata(self.ypos)        
                    
                for line in self.linelist2:
                    line.set_ydata([event.ydata, event.ydata])        

                for line in self.linelist3:
                    line.set_xdata([event.xdata, event.xdata])        
                    
        elif event.name == 'key_press_event':
            # delete selection of points by ctrl+d
            if event.key == 'ctrl+d' or event.key == 'ctrl+D':
                self.xpos.clear()
                self.ypos.clear()

                for line in self.linelist:
                    line.set_xdata(self.xpos)
                    line.set_ydata(self.ypos)                        

                for line in self.linelist2:
                    line.set_ydata([-1, -1])        

                for line in self.linelist3:
                    line.set_xdata([-1, -1])

            elif event.key == 'ctrl+p' or event.key == 'ctrl+P':
            # save selection of points by ctrl+p
                if self.cutzlist is None:
                    print("inspect_channels: No list for"
                        " storing z-cut points specified.")
                else:
                    print("inspect_channels: List of z-cut"
                        " points updated.")
                    self.cutzlist.extend(
                        list(zip(self.xpos, self.ypos)).copy())
            
        # redraw canvas
        self.fig.canvas.draw()
            
# show cuts given from the input, provided that showcut=True
def _show_init_cuts(axes, cutlist, linestyle, linecolor,
    transposectrlpts):            
        
    if not cutlist:
        print("inspect_channels: Invalid cuts to show.")
        return
    elif len(axes[0]) != len(cutlist):
        print("inspect_channels: Unequal number of cuts and "
            "channels to show. No cuts therefore shown.")
        return

    # show each cut        
    for i, cut in enumerate(cutlist):
        
        if transposectrlpts:
            zpos, xpos, ypos = cut
        else:
            zpos = [z for z, x, y in cut]
            xpos = [x for z, x, y in cut]
            ypos = [y for z, x, y in cut]
    
        axes1, axes2, axes3 = axes[0][i], axes[1][i], axes[2][i]
        
        # show lines of the cuts
        axes1.plot(xpos, ypos, ls=linestyle, lw=2,
            c=linecolor, scalex=False, scaley=False)
        axes2.plot(zpos, ypos, ls=linestyle, lw=2,
            c=linecolor, scalex=False, scaley=False)
        axes3.plot(xpos, zpos, ls=linestyle, lw=2,
            c=linecolor, scalex=False, scaley=False)

        # show annotation    
        for z, x, y in zip(zpos, xpos, ypos):
            axes1.annotate(int(z), (x, y))    

# event handler for picking points in the image
class _CutArbPoint:
    def __init__(self, fig, axes1, linelists, layer, numl,
        cutlist, transposectrlpts):
        self.linelist, self.linelist2, self.linelist3 = linelists
        self.xpos, self.ypos, self.zpos = [], [], []
        self.axes1, self.cutlist = axes1, cutlist
        self.layer, self.fig = layer, fig
        self.transposectrlpts = transposectrlpts
        
        # list of lists of annotations for all imag1 images
        self.ans1 = [[]]*numl
        
        self.cid1 = fig.canvas.mpl_connect(
            'button_press_event', self)
        self.cid2 = fig.canvas.mpl_connect(
            'key_press_event', self)

    def __call__(self, event):
        if event.name == 'button_press_event':
            # choose control point by left mouse button
            if event.inaxes in self.axes1 and \
                event.button == 1 and \
                event.xdata is not None and \
                event.ydata is not None:
                self.xpos.append(event.xdata) # event.xdata ???
                self.ypos.append(event.ydata) # event.ydata ???
                self.zpos.append(self.layer[0])
                    
                for line in self.linelist:
                    line.set_xdata(self.xpos)
                    line.set_ydata(self.ypos)        
                    
                for line in self.linelist2:
                    line.set_xdata(self.zpos)
                    line.set_ydata(self.ypos)        

                for line in self.linelist3:
                    line.set_xdata(self.xpos)
                    line.set_ydata(self.zpos)        
                
                for an1, ax in zip(self.ans1, self.axes1):
                    an = ax.annotate(self.layer[0],
                        (event.xdata, event.ydata))    
                    an1.append(an)
                    
        elif event.name == 'key_press_event':
            # delete selection of points by ctrl+d
            if event.key == 'ctrl+d' or event.key == 'ctrl+D':
                self.xpos.clear()
                self.ypos.clear()
                self.zpos.clear()

                for line in self.linelist:
                    line.set_xdata(self.xpos)
                    line.set_ydata(self.ypos)                        

                for line in self.linelist2:
                    line.set_xdata(self.zpos)
                    line.set_ydata(self.ypos)        

                for line in self.linelist3:
                    line.set_xdata(self.xpos)
                    line.set_ydata(self.zpos)        
                
                for lay in self.ans1:
                    for an1 in lay:
                        an1.remove()
                    lay.clear()
            elif event.key == 'ctrl+p' or event.key == 'ctrl+P':
            # save selection of points by ctrl+p
                if self.cutlist is None:
                    print("inspect_channels: No list for storing"
                        " cut control points specified.")
                else:
                    print("inspect_channels: List of cut control"
                        " points updated.")

                    if self.transposectrlpts:
                        newcut = [
                            self.zpos.copy(),
                            self.xpos.copy(),
                            self.ypos.copy()]
                    else:
                        newcut = [np.array([z, x, y])
                            for z, x, y in
                            zip(self.zpos, self.xpos, self.ypos)]
                    self.cutlist.append(newcut)

        # redraw canvas
        self.fig.canvas.draw()
        
# event handler for motion between layers
class _MoveCutsVer:
    def __init__(self, fig, axes, arrs, xind, yind, layer,
        unitdir, t2, t1, imglist, v2, h3, numl, rangelayer):
        self.imglist, self.arrs = imglist, arrs

        self.xind, self.yind = xind, yind
        self.t1, self.t2 = t1, t2                
        self.v2, self.h3 = v2, h3
        self.layer = layer
        self.axes, self.numl = axes, numl
        self.fig, self.rangelayer = fig, rangelayer
        self.unitdir = unitdir
        
#        self.ccc = 1
        
        self.cid = fig.canvas.mpl_connect('scroll_event', self)

    def __call__(self, event):
        # adjust layer according to mouse scroll
        
#        print("NO TY VOLE!!! ", self.ccc)
#        self.ccc += 1
        
        # update layer
        self.layer[0] += event.step
        
        # if self.layer[0] moves out of valid bounds, do nothing
        if self.layer[0] >= self.numl:
            self.layer[0] = self.numl - 1
            return
        if self.layer[0] < 0:
            self.layer[0] = 0
            return

        # update layer number
        for text in self.t1:
            text.set_text(
                "layer no. {}".format(self.layer[0] % self.numl))
            
        # update value at cursor
        for text, (idnum, chan, val) in zip(self.t2, self.arrs):
            key = (idnum, chan)
            xind = -1 if self.xind[0] >= val.shape[1] \
                else self.xind[0]
            yind = -1 if self.yind[0] >= val.shape[2] \
                else self.yind[0]
            text.set_text(_str_value_at_cursor(
                val[self.layer[0]],
                self.unitdir[key], xind, yind))

        # adjust scaling of values for given layer
        if self.rangelayer == 'each':     
            for img, (idnum, chan, val) \
                in zip(self.imglist, self.arrs):              
                imgdata = val[self.layer[0]]
                img.set_data(imgdata)
                img.set_clim(vmin=np.nanmin(imgdata))
                img.set_clim(vmax=np.nanmax(imgdata))
        else:          
            for img, (idnum, chan, val) \
                in zip(self.imglist, self.arrs):   
                img.set_data(val[self.layer[0]])

        # update lines in cut images marking the layer position
        for verl, horl in zip(self.v2, self.h3):
            verl.set_xdata([self.layer[0], self.layer[0]])
            horl.set_ydata([self.layer[0], self.layer[0]])

        # TODO: místem největšího zpomalení vykreslování je následující
        # příkaz; nejdéle tedy trvá věc překreslit na obrazovku; samotné
        # vytváření řezů či aktualizace dat pro čáry není v porovnání tak
        # náročné:
        self.fig.canvas.draw()

# event handler for motion within a single layer
class _MoveCutsHor:
    def __init__(self, fig, axes, arrs, xind, yind, layer,
        unitdir, t2, img2, img3, lines):
        self.t2, self.arrs = t2, arrs
        self.xind, self.yind = xind, yind
        self.vl, self.hl, self.hl2, self.vl3 = lines
        self.img2, self.img3 = img2, img3
        self.axes1, self.axes2, self.axes3 = axes
        self.layer = layer
        self.unitdir, self.fig = unitdir, fig

        self.cid = fig.canvas.mpl_connect(
            'motion_notify_event', self)

    def __call__(self, event):
        # mouse cursor position
        if event.inaxes in self.axes1:
            self.yind[0] = math.floor(event.xdata)  # BEWARE! yind vs. xdata
            self.xind[0] = math.floor(event.ydata)  # BEWARE! xind vs. ydata
        elif event.inaxes in self.axes2:
            self.xind[0] = math.floor(event.ydata)  # BEWARE! xind vs. ydata
        elif event.inaxes in self.axes3:
            self.yind[0] = math.floor(event.xdata)  # BEWARE! yind vs. xdata
        
        # update value at cursor and cuts
        for text, img2, img3, (idnum, chan, val) \
            in zip(self.t2, self.img2, self.img3, self.arrs):
            key = (idnum, chan)
            xind = -1 if self.xind[0] >= val.shape[1] \
                else self.xind[0]
            yind = -1 if self.yind[0] >= val.shape[2] \
                else self.yind[0]
        
            # update value at cursor
            text.set_text(
                _str_value_at_cursor(
                val[self.layer[0]], self.unitdir[key],
                xind, yind))
        
            # update cuts
            cuts = _create_cuts(val, xind, yind)
            img2.set_data(cuts[0])
            img3.set_data(cuts[1])

        # update lines
        for line in self.vl:
            line.set_xdata([self.yind[0], self.yind[0]])
            
        for line in self.hl:
            line.set_ydata([self.xind[0], self.xind[0]])

        for line in self.vl3:
            line.set_xdata([self.yind[0], self.yind[0]])
            
        for line in self.hl2:
            line.set_ydata([self.xind[0], self.xind[0]])

        # redraw canvas
        self.fig.canvas.draw()

# main procedure
def inspect_channels(spmdatalist, chanlist, cutlist=None,
    cutzlist=None, layer=0, xind=None, yind=None,
    rangelayer='all', linestyle='-', linecolor='w',
    pointcolor='g', aspect='auto', size="30%",
    interpolation='none', scaling=True, showcut=False,
    showconstants=False, transposectrlpts=True,
    adjustorient=False):
    """show one or more channels for various SPMData structures
    together with its cuts 

    Instead of SPMData structure it suffices to use any class
    with attributes:
        channels - dictionary {chan:np.ndarray, ...}
        xsteps - dictionary {chan: , ...}
        ysteps - dictionary {chan: , ...}
        units - dictionary {chan: , ...]
        numlay - integer, number of layers in np.ndarrays

    spmdatalist - list of SPMData structures used for plotting
    chanlist - list of channels to be drawn; for each index 'i'
        channel chanlist[i] is shown for SPMData structure
        spmdatalist[i]
    cutlist - list into which cut control points are stored;
        if None, then nothing is stored, otherwise after each
        selection is made by mouse and saved by ctrl+p new
        triple [zpos, xpos, ypos] is appended to cutlist so
        that these data are accessible even after the
        termination of inspect_channels routine; cutlist is 
        finally of the form (if transposectrlpts=True):
        cutlist = [[zpos1, zpos2, ...], [xpos1, xpos2, ...], \
            [ypos1, ypos2, ...]]
    cutzlist - list into which z-cut points are stored;
        if None, then nothing is stored, otherwise after each
        selection is made by mouse right-button new double
        [xpos, ypos] is appended to cutzlist so that these
        data are accesssible even after the termination of
        inspect_channels routine
    layer - which layer of all the channels in chanlist to
        draw at the beginning, layer is identical for all 
        channels in chanlist
    xind, yind - initial position of the intersection point
    linestyle - format specification of the intersecting lines
    linecolor - color specification of the intersecting lines
    pointcolor - color specification of the z-cut points and
        corresponding lines
    rangelayer - which layer to choose for image min and max
        determination;
        if 'all', then global minima and maxima are used;
        if 'each', then maxima and minima are determined
        independently for each layer;
        if rangelayer is not a valid number of a layer,
        then 'layer' is used instead
    aspect - aspect of the cuts in the figure
    size - relative size of cuts in the figure
    interpolation - how to interpolate values in imshow
    scaling - if True, then image labels are chosen
        to correspond to real physical dimensions of the data;
        if False, then array indices are used instead
    showcut - if True and cutlist is not empty, then
        their values are shown in the images
    showconstants - if False, then each channel, whose array has
        every element equal to the same number, is excluded     
    transposectrlpts - if True, then control points (both input
        points and output points) are treated in a transposed
        form [zlist, xlist, ylist], where zlist stores z-coords.
        of all control points etc.; if False, then control points
        are treated in a form [point1, ..., pointN], where each
        points is of the form [zcoord, xcoord, ycoord]
    adjustorient - if True, then the main images are transposed
        and rotated so that it is easy to retrieve data from them
        DODELAT TOHLE, ABYCH PAK MOHL SPRAVNE ZOBRAZOVAT DATA A
        PRITOM KONZISTENTNE PRACOVAT S REFERENCNIMI BODY
    """

    # if 'spmdatalist' is empty, halt
    if not spmdatalist:
        print("inspect_channels: No channels specified.")
        return

    # if 'chanlist' contains invalid channels, halt
    for spmdata, chan in zip(spmdatalist, chanlist):
        if chan not in spmdata.channels.keys():
            print(("inspect_channels: Invalid or absent "
                "channel '{}'.").format(chan))
            return

    # if 'layer' out of valid range, halt
    numl = min([spmdata.numlay for spmdata in spmdatalist])
    if -numl > layer or layer >= numl:
        print(("inspect_channels: Invalid layer number "
            "{}.").format(layer))
        return

    # create image data and dic. of corresponding units
    arrs = [(id(spmdata), chan, spmdata.channels[chan])
        for spmdata, chan in zip(spmdatalist, chanlist)]
    unitdir = {(id(spmdata), chan):spmdata.units[chan]
        for spmdata, chan in zip(spmdatalist, chanlist)}
        
    # if some array is constant, e.q. each value is identical (up
    # to a round-off errors of the order of TOLERANCY), exclude
    # this array from arrs list, since color function does not
    # like these arrays and RuntimeWarnings are thrown
    TOLERANCY = 1e-15
    arrstodel = [False]*len(arrs)
    if not showconstants:
        for i, item in enumerate(arrs):        
            minaux = np.nanmin(item[2])
            arraux = np.abs(item[2] - minaux)
            arraux[np.isnan(arraux)] = 0
            if (arraux < TOLERANCY).all():
                print(("inspect_channels: Array for channel '{}' "
                    "full of values '{}' and is excluded from "
                    "showing.").format(item[1], minaux))
                arrstodel[i] = True
        arrs = [item for i, item in enumerate(arrs)
            if not arrstodel[i]]
        
        if showcut and cutlist:
            cutlist = [item for i, item in enumerate(cutlist)
                if not arrstodel[i]]
                
        if not arrs:
            print("inspect_channels: No channel to show.")
            return
        
    # real physical limits as labels for the image
    if scaling:
        rightlim = {(id(spmdata), chan):(
            spmdata.xsteps[chan] * spmdata.channels[chan].shape[1]
            )
            for spmdata, chan in zip(spmdatalist, chanlist)}
        bottolim = {(id(spmdata), chan):(
            spmdata.ysteps[chan] * spmdata.channels[chan].shape[2]
            )
            for spmdata, chan in zip(spmdatalist, chanlist)}
    else:
        rightlim = None
        bottolim = None
    
    # image limits
    xlim = min([val.shape[1] for idnum, chan, val in arrs])
    ylim = min([val.shape[2] for idnum, chan, val in arrs])

    # initial position of the point of intersection
    if xind is None or -xlim > xind or xind >= xlim:
        xind = xlim // 2
    if yind is None or -xlim > yind or yind >= ylim:
        yind = ylim // 2

    # create initial cuts
    cuts = {(idnum, chan): _create_cuts(val, xind, yind)
        for idnum, chan, val in arrs}

    # determine min and max
    if rangelayer == 'all':
        vextr = {(idnum, chan): (np.nanmin(val), np.nanmax(val))
            for idnum, chan, val in arrs}
    elif rangelayer == 'each':
        vextr = {(idnum, chan): (None, None)
            for idnum, chan, val in arrs}
    else:
        if -numl > rangelayer or rangelayer >= numl:
            print("inspect_channels: Invalid rangelayer number.")
            return
        vextr = {(idnum, chan):(
            np.nanmin(val[rangelayer]),
            np.nanmax(val[rangelayer]))
            for idnum, chan, val in arrs}

    # create figure and layout    
    fig, axes, lines, imgs, texs, pois, zcus = _initialize_plots(
        arrs, unitdir, vextr, layer, interpolation, xind, yind,
        linestyle, linecolor, pointcolor, size, aspect, cuts,
        numl, rightlim, bottolim)
    
    layer = [layer] # layer, now the LIST, is shared by handlers
    xind, yind = [xind], [yind] # ...and similarly for xind, yind
    
    # initialize event handlers
    _CutZedPoint(fig, axes[0], zcus, cutzlist)   
    _CutArbPoint(fig, axes[0], pois, layer, numl, cutlist,
        transposectrlpts)
    _MoveCutsHor(fig, axes, arrs, xind, yind, layer, unitdir,
        texs[1], imgs[1], imgs[2], [lines[0], lines[3], lines[4],
        lines[2]])
    _MoveCutsVer(fig, axes, arrs, xind, yind, layer, unitdir,
        texs[1], texs[0], imgs[0], lines[1], lines[5], numl,
        rangelayer)

    # show initial control points
    if showcut and cutlist:
        _show_init_cuts(axes, cutlist, linestyle, linecolor,
            transposectrlpts)

    # improve layout and launch the interactive window
    fig.tight_layout()
    plt.show()

