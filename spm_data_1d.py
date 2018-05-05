from .spm_data_cuts import *

print("ahoj")

class LineData:
    """1D data
    """
    
    def __init__(self):
        self.channels = {}
        self.num = 0
        self.data = None
        self.metadata = {}
        
    def load_data(self, filename, header=True):
        """load data from external file
        
        header - if True, then only a header is read
        """
    
        with open(filename, 'r') as f:
            nrows = 1
            line = f.readline().strip()
            while line != '[DATA]': 
                line = line.split()
                if not line:
                    pass
                elif line[0] == 'Date':
                    self.metadata['date'] = line[1]
                    self.metadata['time'] = line[2]
                elif line[0] == 'Z' and \
                    line[1] == 'Spectroscopy>Channels':
                    line = ' '.join(line[2:]).split(';')
                    self.metadata['channels'] = line
#                elif ...
#                 SEM KDYZTAK PRIDAT DALSI POLOZKY

#                try:
#                    self.metadata[key] = float(val)
#                except ValueError:
#                    print("load_data: Cannot convert value into float.")
                    
                nrows += 1
                line = f.readline().strip()

            chans = f.readline().strip().split('\t')
            nrows += 1
        
        if not header:
            self.data = np.loadtxt(filename, skiprows=nrows).T                 
            quandict = {'Z rel (m)':ZREL, 'Current (A)':CURR, 'Z (m)':ELEV,
                'Phase (deg)':PHAS, 'Amplitude (m)':AMPL, 'Frequency Shift (Hz)':FREQ,
                'Excitation (V)':EXCI, 'Current [bwd] (A)':CURRB, 'Z [bwd] (m)':ELEVB,
                'Phase [bwd] (deg)':PHASB, 'Amplitude [bwd] (m)':AMPLB,
                'Frequency Shift [bwd] (Hz)':FREQB, 'Excitation [bwd] (V)':EXCIB
                }
                     
            for i, chan in enumerate(chans):
                if chan not in quandict:
                    print("load_data: Invalid channel.")
                else:
                    self.channels[quandict[chan]] = self.data[i]
            
    def show_channels(self, *cl, linespec="b-"):
        """plot 1D data of channels
        
        cl - list of channels to plot
        linespec - format specification of lines
        """

        # pick only valid channels
        chanlist = [chan for chan in cl if chan in self.channels]
            
        if len(chanlist) == 0:
            print("show_channels: No valid channels.")
            return
            
        # create figure
        fig, subfiglist = plt.subplots(1, len(chanlist), squeeze=False, num="1D Data - channels: " + ", ".join(chanlist))
        subfiglist = subfiglist[0]
        
        # plot each valid channel
        for i, chan in enumerate(chanlist):
            datal = range(len(self.channels[chan]))
            subfiglist[i].plot(datal, self.channels[chan], linespec)
            subfiglist[i].set_xlim([0, len(datal)])
            subfiglist[i].set_title("Channel: {} [units: {}]".format(chan, UNITS[chan]))
            subfiglist[i].set_xlabel("No. of steps")
            subfiglist[i].set_ylabel("{}".format(UNITS[chan]))
        
        # show
        fig.tight_layout()
        plt.show()







if __name__ == "__main__":
    print("LINE TEST MODE")
    
    filename='/home/jaroslav/Plocha/FZU/MartinSvecProgramy/spm_data_3/Spectroscopy/Z_Spectroscopy_044.dat'
    
    a = LineData()
    a.load_data(filename)

#    print(a.data[:5])
#    print(a.channels[ZREL][:5])
    
#    a.show_channels(AMPL, ZREL, PHAS)
    a.show_channels(AMPL)
