from spm_data_cut2d import *
from spm_data_struct import *
from spm_data_inspect import inspect_channels

if __name__ != "__main__":
    
    plt.ion()

else:
    print("Test mode")
    b = SPMdata()
    bb = SPMdata()
#    b.load_data("/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration_inequal", key=lambda x: x.filename)
#    b.load_data("/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration/data_charged/lorange", key=lambda x: x.filename)

    b.load_data("/home/jaroslav/Plocha/FZU/MartinSvecProgramy/spm_data_3/1/test")#, key=lambda x: x.height)

#    b.load_data("/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration/data_uncharged", key=lambda x: x.filename)
#    b.load_data("/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration/data_charged_2/", key=lambda x: x.filename)
#    bb.load_data("/home/jaroslav/Plocha/FZU/MartinSvecProgramy/MartinSvecData/4jarda/registration/data_charged/shrange", key=lambda x: x.filename)
    
    
    
#    bb = a.copy(CURR, FREQ, CURRB)
#    bb.align_corr_chan_fac(CURRB, FREQ, ichan=CURRB, downsample=True, fac=3, incr=True)#, show=False, rellay=range(5, 10), mislay=range(11, 21))
#    

    a=b.copy(FREQ, FREQB, CURR, CURRB)#, ll=range(10))
    a.show_layer(FREQ,ilay=3)
    
#    del b

#    a.add_layers([1,8,9,10])

    print(a.xoffinds)
    a.align_corr_chan_fac(CURRB, FREQ, ichan=CURRB, downsample=True, fac=3)#, show=False, rellay=range(5, 10), mislay=range(11, 21))


    a.inspect_channels(FREQ)
    
    
#    ADHOC = "ahoj"

##    arr=np.ones_like(a.freq)
##    arr=np.zeros((len(a), a.layers[0].xnum, a.layers[0].ynum + 20))
##    for i, lay in enumerate(a.layers):
##        xoff, yoff = lay.xoffind[FREQ], lay.yoffind[FREQ]
##        xnum, ynum = lay.xnum, lay.ynum
##        arr[i, xoff:xoff + xnum, yoff:yoff + ynum] = 1

##    a.add_channel(arr, ADHOC)
#    
##    aver = average_layers(a, range(5))
##    print("aver.channels.keys(): ", aver.channels.keys())    
##    aver.show_channels(FREQ, CURR, PHAS)
##    a.channels[ADHOC] -= a.freq

##    a.add_channel(-a.freq +10e9* a.curr, "ahh")
##    print("a.xoffinds: ", a.xoffinds)


#    aa = create_spm(numlay=a.numlay)
#    aa.print_channels()
#    aa.print_nums()
#    
##    print("aa.layers: ", aa.layers)
#    
##    for i, lay in enumerate(aa):
###        lay.height = -i
##        setattr(aa, "layers[{}].height".format(i), i)
#        
#    for i in range(len(aa)):
##        print("i = ", i)
#        aa.layers[i].height = 50-i
##        print("HA: ", aa.layers[i].height)

##    for i in range(len(aa)):
##        print("HA: ", aa.layers[i].height)
##        
##    for i, lay in enumerate(aa):
##        print("lay.height: ", lay.height)
#    
##    view = Cut2D()            
##    a.view_channel(FREQ, view, key=lambda x: -x.height)
##    view.show_channels(FREQ)
##    a.inspect_channels(FREQ)                
#                
#    print("a.xoffinds: ", a.xoffinds[FREQ])

#    a.show_layer(FREQ, ilay=3, rangelayer='all')
#                
#    aa.add_channel(ADHOC, -a.freq, xofflist=a.xoffinds[FREQ], yofflist=a.yoffinds[FREQ], xnum=a.xnums, ynum=a.ynums)
#    
#    print("aa.xoffinds: ", aa.xoffinds[ADHOC])
#    aa.sort_spm()
#    
##    aa.print_nums()
##    
##    aa.rewrite_channel(ADHOC, np.floor(10*a.freq))
#    aa.rewrite_channel(ADHOC, np.floor(10*aa.ahoj))
#    
#    a = aa
#    a.check_valid_frames()
#    
##    print("a.xoffinds: ", a.xoffinds)
##    a.print_heights()
##    a.print_ranges()

##    print("heigh", a.layers[0].height)

##    a = stack_spm(a, bb)
##    a.align_corr_chan_fac(ichan=CURR, downsample=True, ilay=2, fac=2, show=False, order=1)#, actlay=range(0,8))#, rellay=range(0, 8), mislay=range(8, 16))
#    
#    a.show_layer(ADHOC, ilay=3, rangelayer='all')
##    a.show_layer(FREQ, ilay=-3)
#    a.inspect_channels(ADHOC, rangelayer='all', xind=10, yind=10)
#    a.delete_layers(list(range(0,10)))
##    print("len(a): ", len(a))
##    a.inspect_channels(ADHOC)
##    a.show_layer(ADHOC, ilay=3, rangelayer='all')
#    
##    plt.imshow(a.layers[0].freq)
##    plt.show()

##    a.show_layer(FREQ, ADHOC, 5)
#    
#    laylist = list(range(10, 26))
#    cl, clz = [], []
#        
##    c = stack_spm(a, bb)
##    c.show_channel(CURR, laylist, limit=16)
##    c.show_channel(EXCIB, laylist, limit=16)
##    c.show_channel(ADHOC, laylist, limit=16)
##    c.inspect_channels(FREQ, ADHOC, rangelayer='all', cutlist=cl)
##    inspect_channels([bb,c],[FREQ, ADHOC])

##    c.insp(FREQ, ADHOC, rangelayer='all', cutlist=cl, cutzlist=clz)

#    print("cl, clz: ", cl, clz)
#    
##    help(SPMdata)
##    print(dir())

##    a.show_layer(FREQ, ilay=0)

##    a.print_channels()
##    a.align_forback(ichan=FREQ, ilay=0, pad=True)    
##    a.align_offsets()
##    a.align_corr_chan_fac(FREQ, FREQB, ichan=CURR, downsample=True, fac=3)#, show=False, rellay=range(5, 10), mislay=range(11, 21))



##    a.align_corr_chan_fac(ichan=CURR, downsample=True, ilay=2, fac=2, show=True, order=1, nlay=8, rnlay=7, actlay=range(8,16))#, rellay=range(0, 8), mislay=range(8, 16))


##    a.align_corr_chan_fac(ichan=CURR, downsample=True, fac=3, show=True)#, rellay=range(5, 10), mislay=range(11, 21))
##    a.align_trend(ichan=CURR, fac=3, show=True, order=3, rellay=range(1,9), mislay=range(11,38))
##    a.hew_nans()
#    
##    a.align_corr_chan_fac(FREQ, FREQB, ichan=CURR, downsample=False, fac=4, show=False, rellay=range(5, 10), mislay=range(11, 21))
##    a.align_corr_chan_fac(ichan=CURR, downsample=False, fac=4, show=False, rellay=range(5, 10), mislay=range(11, 21))
##    a.show_layer(FREQ, ilay=0)
#    
##    for i in range(len(a)):
##        count_nans(a.channels[FREQ][i][50])
#    
##    a.align_correlate(EXCI, FREQ, ilay=0)
##    a.inspect_channels(FREQ, FREQB)

##    a.align_forback(ilay=0, same=True)

#    cutlist, cutzlist = [], []
##    a.align_forback(FREQ, PHAS, EXCI, CURR, ichan=FREQ, ilay=0, pad=True)


##    a.inspect_channels(FREQ, ADHOC, cutlist=cutlist, cutzlist=cutzlist, rangelayer=2)
##    a.inspect_channels(CURR, FREQ, ADHOC, cutlist=cutlist, cutzlist=cutzlist, rangelayer='all')

#    print("cutlist: ", cutlist)
#    
##    a.strip_nans(same=True)
##    a.hew_nans(same=False)
#    
##    a.align_forback(FREQ, PHAS, EXCI, CURR, ichan=FREQ, ilay=0, pad=True)
##    print("a.xdims: ", [chan.shape for chan in a.channels.values()])

##    print("a.arrshape: ", a.arrshape)

##    a.align_corr_chan_fac(FREQB, ichan=CURR, fac=1)#, show=True, mislay=range(5, 11))#, corrfunc="_align_fft_upsample")#, incr=True, xcoefs=(-0.5, 0), ycoefs=(0.5, 0))
##    a.align_corr_chan_fac(FREQ, CURR, ichan=CURR, fac=4, downsample=True, show=True, rellay=range(5,10), mislay=range(11, 21))#, corrfunc="_align_fft_upsample")#, incr=True, xcoefs=(-0.5, 0), ycoefs=(0.5, 0))

##    gg=copy.deepcopy(a)
##    gg.show_channels(FREQ, CURR)
#        
#    cutlist = [[[0, 9, 15], [26.510958704009717, 74.315435667310027, 116.46257807874049], [4.3459164566023674, 110.98667275934925, 14.529118649968126]], [[0, 9, 15], [26.510958704009717, 74.315435667310027, 116.46257807874049], [4.3459164566023674, 110.98667275934925, 14.529118649968126]]]

##    a.inspect_channels(FREQ, FREQB, cutlist=cutlist, cutzlist=cutzlist)#, interp='bilinear')
#    a.inspect_channels(PHAS, PHASB, cutlist=cutlist, cutzlist=cutzlist)#, interp='bilinear')
#    a.inspect_channels(EXCI, EXCIB, cutlist=cutlist, cutzlist=cutzlist)#, interp='bilinear')
#    a.inspect_channels(CURR, CURRB, cutlist=cutlist, cutzlist=cutzlist)#, interp='bilinear')
    
#    e=Cut1D()
###    e.create_cutzxy(a, cutlist[0][0], cutlist[0][1], cutlist[0][2])
#    e.create_cutzxy(a, *cutlist[0])
##    e.show_channels(FREQ, CURR, rescale=True)    
#    e.create_zcut(a, cutzlist[0], deltaz=1, interpolate=0)
###    print("e.channels[FREQ]: ", e.channels[FREQ])
####    print("min, max: ", min(e.channels[FREQ]), max(e.channels[FREQ]))
#    e.show_channels(FREQ, CURR)

#    print("2D CUTS")
#    ff=Cut2D()
#    ff.create_cutzxy(a, *cutlist[0], interpolate=1, align=True, delta=0.5, entire=True)
##    ff.rotate(10)
##    ff.create_cut_rotate(a, *cutlist[0])
####    ff.create_cut(a, cutlist)

#    ff.show_channels(FREQ, CURR, rescale=True)
#    ff.strip_nans(ref=FREQ)
#    ff.show_channels(FREQ, CURR, rescale=True)
#    
##    print("ff.channels[FREQ]: ", ff.channels[FREQ])
#    
##    pickle_data(a, "/home/jaroslav/Plocha/FZU/MartinSvecProgramy/spm_data/ahoj.data")
##    del a    
##    b = [None]
##    unpickle_data(b, "/home/jaroslav/Plocha/FZU/MartinSvecProgramy/spm_data/ahoj.data")
#    b = b[0]
#    b.show_channel(PHAS)
    
#    a.inspect_channels(FREQ)
#    print([(b,c) for (a, b, c) in a.channels_fb.values()])
#    a.show_channel(PHAS,l=[1,2,3])
