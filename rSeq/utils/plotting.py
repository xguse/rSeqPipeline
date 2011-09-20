import os
from scipy.stats import gaussian_kde
from numpy.random import normal
from numpy import arange
import numpy as np
from gusPyCode.defs import packings
import matplotlib as mpl



def autoPickColors(number):
    """Generate a set of colors with near optimal contrasts."""
    packingsPath = os.path.dirname(os.path.realpath(packings.__file__))
    
    coords = map(lambda l: l.strip('\n'), open("%s/pack.3.%s.txt" % (packingsPath,number), 'rU'))
    coords = np.array(coords, dtype=np.float)
    coords = coords.reshape(number,3)
    coords = coords * 0.5
    coords = coords + coords.max()
    
    colors = []
    for c in coords:
        colors.append(mpl.colors.rgb2hex(np.absolute(c)))
        
    return tuple(colors)
    
def setTickSizes(axObj,fontSize):
    for label in axObj.xaxis.get_ticklabels():
        # label is a Text instance
        #label.set_color('red')
        #label.set_rotation(45)
        label.set_fontsize(fontSize)
        
    for label in axObj.yaxis.get_ticklabels():
        # label is a Text instance
        #label.set_color('red')
        #label.set_rotation(45)
        label.set_fontsize(fontSize)

def violin_plot(ax,data,pos, bp=False):
    '''(http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html)
    create violin plots on an axis
    '''
    dist = max(pos)-min(pos)
    w = min(0.15*max(dist,1.0),0.5)
    for d,p in zip(data,pos):
        if d == 0: # WAD: handles cases of no data
            continue
        k = gaussian_kde(d) #calculates the kernel density
        m = k.dataset.min() #lower bound of violin
        M = k.dataset.max() #upper bound of violin
        x = arange(m,M,(M-m)/100.) # support for violin
        v = k.evaluate(x) #violin profile (density curve)
        v = v/v.max()*w #scaling the violin to the available space
        ax.fill_betweenx(x,p,v+p,facecolor='y',alpha=0.3)
        ax.fill_betweenx(x,p,-v+p,facecolor='y',alpha=0.3)
    if bp:
        ax.boxplot(data,notch=1,positions=pos,vert=1)