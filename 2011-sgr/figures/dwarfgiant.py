import atpy
import pylab as plt
from numpy import *
import os
from scipy import interpolate
import matplotlib
plt.rcParams['text.usetex'] = True

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: Number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = linspace(0,1.,N)
    # N+1 indices
    indices = linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)
  

def cmap_map(function,cmap):
    """ Applies function (which should operate on vectors of shape 3:
    [r, g, b], on colormap cmap. This routine will break any discontinuous     points in a colormap.
    """
    cdict = cmap._segmentdata
    step_dict = {}
    # Firt get the list of points where the segments start or end
    for key in ('red','green','blue'):         step_dict[key] = map(lambda x: x[0], cdict[key])
    step_list = sum(step_dict.values(), [])
    step_list = array(list(set(step_list)))
    # Then compute the LUT, and apply the function to the LUT
    reduced_cmap = lambda step : array(cmap(step)[0:3])
    old_LUT = array(map( reduced_cmap, step_list))
    new_LUT = array(map( function, old_LUT))
    # Now try to make a minimal segment definition of the new LUT
    cdict = {}
    for i,key in enumerate(('red','green','blue')):
        this_cdict = {}
        for j,step in enumerate(step_list):
            if step in step_dict[key]:
                this_cdict[step] = new_LUT[j,i]
            elif new_LUT[j,i]!=old_LUT[j,i]:
                this_cdict[step] = new_LUT[j,i]
        colorvector=  map(lambda x: x + (x[1], ), this_cdict.items())
        colorvector.sort()
        cdict[key] = colorvector

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)
    
    
table = atpy.Table('combined.xml', verbose=False)

table = table.where(table.Candidates)

umg = table.data['(g-i)']
mgh = table.data['EW_{MgH}']

plt.close('all')
x = arange(0, 4, 0.1)
y = 3 * x + 0.5 # dividing range


gmiBound = [0.65, 1.6]
ewmgBound = [0, 19.9]


fig = plt.figure()

ax1 = fig.add_subplot(223)

# synthetic spectra
logg, feh, ewmg, gmi = loadtxt('synthetic.data', usecols=(2, 3, 4, 6), unpack=True)
cmap = cmap_discretize(plt.cm.RdBu_r, 6)

g_ewmg = []
g_gmi = []
g_feh = []
d_ewmg = []
d_gmi = []
d_feh = []
for gmi_, ewmg_, feh_, logg_ in zip(gmi, ewmg, feh, logg):
    if logg_ < 1.:
        g_ewmg.append(ewmg_)
        g_gmi.append(gmi_)
        g_feh.append(feh_)
        
    else:
        d_ewmg.append(ewmg_)
        d_gmi.append(gmi_)
        d_feh.append(feh_)
    

scat = ax1.scatter(gmi, ewmg, c=feh, s=15, marker='o', cmap=cmap, vmin=-2.5, vmax=0.5)
scat2 = ax1.scatter(g_gmi, g_ewmg, c=g_feh, marker='s', s=40, cmap=cmap, vmin=-2.5, vmax=0.5)
#scat = ax1.scatter(array(d_gmi), array(d_ewmg), c=array(d_feh), marker='o', s=15, cmap=cmap, vmin=-2.5, vmax=0.5,)
#scat2 = ax1.scatter(g_gmi, g_ewmg, c=g_feh, marker='s', s=30, cmap=cmap, vmin=-2.5, vmax=0.5, edgecolor='none')

ax1.set_xlim(gmiBound)
ax1.set_ylim(ewmgBound)
#ax1.set_xlabel(r'g \--- i')
ax1.set_ylabel(r'Mg I EW (\AA{})')

#ax1.set_yticklabels(['0', '5', '10', '15', ' '])

#clbar.ax.invert_yaxis()

ax2 = fig.add_subplot(224)

ax2.scatter(umg, mgh, s=6, edgecolor='none', facecolor='k')
ax2.set_xlim(gmiBound)
ax2.set_ylim(ewmgBound)
ax1.plot(x,y, 'k-', color='#666666', linewidth=2)
ax2.plot(x,y,'k-', color='#666666', linewidth=2)
ax2.set_xlabel(r'$(g - i)_o$')
ax2.set_ylabel(' ')
ax1.set_xlabel(r'$(g - i)_o$')
ax2.yaxis.set_visible(False)

subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, wspace=0.10, hspace=0.20)

ax3 = fig.add_subplot(221)
clbar = plt.colorbar(scat, cax=ax3, ticks=[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, +0.5], orientation='horizontal')#, shrink=0.8)
clbar.set_ticklabels(['$-2.5$', '$-2.0$', '$-1.5$', '$-1.0$', '$-0.5$', '\,$0.0$', '$+0.5$'], update_ticks=True)
#clbar.set_label('[Fe/H]')
ax3.title.set_text('[Fe/H]')
ax3.title.set_fontsize(12)
ax3.set_position([0.25, 0.55, 0.5, 0.04], which='both')


#clbar2 = plt.colorbar(scat, ax=ax2, ticks=[-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, +0.5], orientation='horizontal')
#clbar2.ax.set_visible(False)


os.system('rm -f dwarfgiant.eps')
plt.savefig('dwarfgiant.eps')
