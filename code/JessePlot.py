# Python3 divison Int/Int = Float
from __future__ import division
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import os


#==============================================================================
# Plotting functions
#==============================================================================



def rc(style='ppt_1', close=True):
    '''A function to set the default parameters for plots, just how I like 'em
    
    Keywords
    ------
    style : String, Options include 'ppt_1', 'ppt_2'
    
    close : Bool, Close all plots
    '''        
    #This all could also be done with plt.rc(), ex:
    #plt.rc('font', **{'family':'sans-serif', 'sans-serif':'Helvetica'})
    
    
    if close:
        plt.close('all')
    else:
        plt.figure()

    #Common starting point for all styles
    plt.rcdefaults()
    
    p = plt.rcParams    
    
    p['lines.linewidth'] = 2
    p['axes.linewidth'] = 2
    p['xtick.minor.width'] = 1
    p['ytick.minor.width'] = 1
    p['xtick.major.width'] = 2
    p['ytick.major.width'] = 2
    p['xtick.minor.size'] = 5
    p['ytick.minor.size'] = 5
    p['xtick.major.size'] = 10
    p['ytick.major.size'] = 10
    p['figure.dpi'] = 80 
    p['savefig.dpi'] = 300 


    if style == 'ppt_1':
        p['font.size'] = 20
        p['axes.labelsize'] = 24
        p['figure.figsize'] = (8, 6) #w,h

    elif style == 'ppt_2':
        p['font.size'] = 16
        p['axes.labelsize'] = 20
        p['figure.figsize'] = (5, 6) #w,h  

    elif style == 'ppt_4':
        p['font.size'] = 16
        p['axes.labelsize'] = 20
        p['figure.figsize'] = (5, 3) #w,h  

    elif style == 'ppt_6':
        p['font.size'] = 16
        p['axes.labelsize'] = 20
        p['figure.figsize'] = (3, 3) #w,h 
        
    elif style == 'word_3':
        p['lines.linewidth'] = 1.5
        p['axes.linewidth'] = 1.5
        p['xtick.minor.width'] = 1
        p['ytick.minor.width'] = 1
        p['xtick.major.width'] = 1.5
        p['ytick.major.width'] = 1.5
        p['xtick.minor.size'] = 2
        p['ytick.minor.size'] = 2
        p['xtick.major.size'] = 5
        p['ytick.major.size'] = 5
        p['figure.dpi'] = 300 
        p['savefig.dpi'] = 900 
        p['font.family'] = 'sans-serif'
        p['font.sans-serif'] = ['Arial']
        p['font.weight'] = 'bold'
        p['axes.labelweight'] = 'bold'
        p['font.size'] = 10
        p['legend.fontsize']= 10
        p['axes.labelsize'] = 12
        p['figure.figsize'] = (29/12, 2) #w,h 
        
    elif style == 'word_2/3':
        p['lines.linewidth'] = 1.5
        p['axes.linewidth'] = 1.5
        p['xtick.minor.width'] = 1
        p['ytick.minor.width'] = 1
        p['xtick.major.width'] = 1.5
        p['ytick.major.width'] = 1.5
        p['xtick.minor.size'] = 2
        p['ytick.minor.size'] = 2
        p['xtick.major.size'] = 5
        p['ytick.major.size'] = 5
        p['figure.dpi'] = 300 
        p['savefig.dpi'] = 900 
        p['font.family'] = 'sans-serif'
        p['font.sans-serif'] = ['Arial']
        p['font.weight'] = 'bold'
        p['axes.labelweight'] = 'bold'
        p['font.size'] = 10
        p['legend.fontsize']= 10
        p['axes.labelsize'] = 12
        p['figure.figsize'] = (29/6, 2) #w,h 


    else:
        print 'Wrong style: [ppt_1, ppt_2, ppt_4, ppt_6, word_3]'


def saveFig(filename, transparent=True, format='png', path='./../figs/'):
    '''Convenience function to make a transparent plot
    Creates a folder at 'path', and puts a figure in it.
    '''
    if path:
        if not os.path.exists(path):
            os.mkdir(path)
        filename = path + filename
    
    plt.savefig(filename, transparent=transparent, format=format)




def legend(title=''):
    plt.legend(shadow = False, loc ='upper left', bbox_to_anchor = (1.0, 1.0), 
    title=title, fontsize = 14)

# Plotting Helper Funcs
def expandAxes(ratio=1.1):
    ax = plt.gca()
    xlimits = ax.get_xlim()
    ylimits = ax.get_ylim()
    dx = ((xlimits[1]-xlimits[0]) * (ratio-1) ) / 2.
    dy = ((ylimits[1]-ylimits[0]) * (ratio-1) ) / 2.
    ax.set_xlim(xlimits[0]-dx, xlimits[1]+dx)
    ax.set_ylim(ylimits[0]-dy, ylimits[1]+dy)






def shiftedColorMap(cmap, data=False, start=0, midpoint=0.5, stop=1.0, 
                    set=True, register=False, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    if type(data) == np.ndarray:
        vmax = data.max()
        vmin = data.min()
        midpoint = 1 - vmax/(vmax + abs(vmin))


    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)

    if set:
        plt.set_cmap(newcmap)

    if register:
        plt.register_cmap(cmap=newcmap)

    return newcmap


#==============================================================================
# Data Loading / Manipulation
#==============================================================================
# Anything that uses globals() can only be run from the module of interest


def loadNPZ(filename, suffix = ''):
    '''
    Load arrays from a .npz pickled file to the global namespace
    '''
    print "Importing...\n"  
    with np.load(filename) as data:
        for key in data.keys():
            name = key + suffix
            print name
            globals()[name] = data[key]


def loadPulseData(filename, suffix = ''):
    """Given a pulse text file, puts the data into the global namespace
    """
    data = np.genfromtxt(filename+'.txt', skip_header=3, names=True,
                   dtype='i8,f8,S5,f8,f8,f8,f8,f8,f8')
    print "Importing...\n"
    for key in data.dtype.fields.keys():
        name = key + suffix
        print name
        globals()[name] = data[key]


def loadDCData(filename):
    """Given a DC sweep text file, puts the data into the global namespace
    """
    data = np.genfromtxt(filename+'.txt', skip_header=3, names=True)
    globals()[filename] = data


def getTxtFileList(folder_path='.'):
    """Returns a list of '.txt' filenames in the folder 'folder_path'
    """
    files = [os.path.splitext(f)[0] for f in os.listdir(folder_path) if os.path.isfile(f) and os.path.splitext(f)[1] == '.txt']
    return files




#==============================================================================
# Data  Manipulation
#==============================================================================

        
def scaleArray(array, limits=(0, 1)):
    return ((array - array.min()) / (array.max()-array.min())) * limits[1] + limits[0]    



#==============================================================================
# Data Plots 
#==============================================================================

    
def plotDC(filename):
    '''Needs loadDCData filename'''
    data = globals()[filename]    
    plt.plot(data['V'], data['I'])
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (A)')
    plt.gca().ticklabel_format(axis='y', style='sci', scilimits=(-2,2))

def semilogyDC(filename):
    '''Needs loadDCData filename'''
    data = globals()[filename]    
    plt.semilogy(data['V'], abs(data['I']))
    plt.xlabel('Voltage (V)')
    plt.ylabel('Current (A)')

def plotPulses( V, R, pulse_starts=np.zeros(1), pulse_stops=np.zeros(1), 
               RVplot=False):
    """Plots a series of gradual set/reset pulses
    """
    if not pulse_stops.any():
        pulse_stops = np.array([R.size])
        
    lengths = pulse_stops - pulse_starts
    print lengths
    

    if RVplot:
        for l, pstart, pstop in zip(lengths, pulse_starts, pulse_stops):
            plt.semilogy(V[pstart:pstop], R[pstart:pstop], 'o-', linewidth=1)
            plt.xlabel('Voltage WL (V)')
            plt.ylabel('Resistance (Ohms)')

    else:
        for l, pstart, pstop in zip(lengths, pulse_starts, pulse_stops):
            plt.semilogy(np.arange(l)+1, R[pstart:pstop])
            plt.xlabel('Pulse Number (n)')
            plt.ylabel('Resistance (Ohms)')


def createRArray(R, pulse_starts=np.zeros(1), pulse_stops=np.zeros(1), array_length=None):
    '''
    Creates an array with a row for each pulse sequence ((pulse_starts.size, array_length))
    Sequences defined by pulse_starts, and pulse_stops (arrays of indices).
    '''

    if not pulse_stops.any():
        pulse_stops = np.array([R.size])
        
    lengths = pulse_stops - pulse_starts
    print lengths

    if not array_length:        
        array_length = max(lengths)

    R_array = np.zeros( (lengths.size, array_length))

    for i, l, pstart, pstop in zip(np.arange(lengths.size), lengths, pulse_starts, pulse_stops):
       R_array[i, 0:l] = R[pstart:pstop]

    R_array[R_array == 0] = np.nan        
    return R_array
        





#==============================================================================
# Specific Plots
#==============================================================================
    
def plotDensity( Px, Pyx, limits=[-5,15], no_ticks=False, fill=True,
                color_list=['r', 'y', 'g', 'c'], plot_Pyx=False, cmap=False):
    '''
    Plots distributions given Px[array], and Pyx[CapacityTools.FunctionList]
    '''
    Py = Pyx.dot(Px)
    PxPyx = Pyx * Px    
    
    fig, ax = plt.subplots()
    
    x = np.linspace(limits[0], limits[1], 1000)
    plt.plot(x, Py(x), 'b')
    
    if plot_Pyx:
        Pplot = Pyx
    else:
        Pplot = PxPyx
        
    Pplot = [p for i, p in enumerate(Pplot) if Px[i] != 0]        
        
    for i, p in enumerate(Pplot):
        if cmap:
            color = cmap(i/len(Pplot))
        else:
            color = color_list[i%len(color_list)]
            
        plt.plot(x, p(x), '--', color=color)

        if fill:
            plt.fill_between(x, p(x), color=color, alpha=0.1)

    if no_ticks:
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        
    plt.show()
    
    
    
def plotDividers(dividers):
    '''
    Add dividers to a plot
    '''
    ax = plt.gca()
    ylim0 = ax.get_ylim()
    
    for d in dividers:
        plt.plot([d,d], [0,1e5], 'k', linewidth=1.5)

    ax.set_ylim(ylim0)                
    

def plotDI(dI, fill=True, color_list=['r', 'y', 'g', 'c']):
    '''
    Add dI overlay to a plot
    '''
    ax = plt.gca()
    xlim0 = ax.get_xlim()    
    x = np.linspace(xlim0[0], xlim0[1], 1000)

    print dI
    for i, p in enumerate(dI):
        plt.plot(x, p(x), '-', color=color_list[i%len(color_list)], linewidth=1.5)
        if fill:
            plt.fill_between(x, p(x), hatch='/', color=color_list[i%len(color_list)], alpha=0.4)
        
    plt.show()



    
def plotDensityDiscrete( Px_discrete, Pyx_discrete ):
    '''
    Plot Discrete Pyx Matrices
    '''
    plt.matshow(np.log10(Pyx_discrete), cmap=plt.cm.Blues)
    plt.xlabel('Write')
    plt.ylabel('Read')
    plt.colorbar(label='$log_{10}$ P(Read|Write)', fraction=0.15, shrink=0.8)
    plt.show()

    