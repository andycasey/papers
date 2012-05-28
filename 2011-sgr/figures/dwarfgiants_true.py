
import atpy
import matplotlib
import matplotlib.pyplot as plt
from numpy import *
from scipy import interpolate

fontsize=15.
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rcParams['font.size'] = fontsize

plt.rcParams['text.usetex'] = True
#plt.rcParams['text.dvipnghack'] = True
plt.close('all')


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
  
  

table = atpy.Table('combined.xml', verbose=False)

subset = '(table.Candidates == True) & (-350. < table.Vgsr) & (table.Vgsr < 350.)'

table = table.where(eval(subset))


logg, feh, ewmg, gmi = loadtxt('synthetic.out', usecols=(2,3,4,6), unpack=True)

giants_gmi = []
giants_feh = []
giants_ewmg = []

dwarfs_gmi = []
dwarfs_feh = []
dwarfs_ewmg = []
for i in range(0, len(logg)):
    if logg[i] == 2.0 and feh[i] in [-0.5, -1.5, -2.5]:
        [typelist.append(dataval[i]) for typelist, dataval in zip([giants_gmi, giants_feh, giants_ewmg], [gmi, feh, ewmg])]

    elif logg[i] == 4.5 and feh[i] in [-0.5, -1.5, -2.5]:
        [typelist.append(dataval[i]) for typelist, dataval in zip([dwarfs_gmi, dwarfs_feh, dwarfs_ewmg], [gmi, feh, ewmg])]


fig = plt.figure()
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95)
ax = fig.add_subplot(111)

x = arange(0, 4, 0.1)
y = 3 * x + 0.5 # dividing range

ax.plot(x, y, '-', color='#333333', linewidth=3, zorder=3)

cmap = cmap_discretize(matplotlib.colors.LinearSegmentedColormap.from_list('mycmp', ['#2165BF', '#F2BA52', '#CE040F'], N=3), 3)
# do the dwarfs first

scat_dwarfs = ax.scatter(dwarfs_gmi, dwarfs_ewmg, s=70., cmap=cmap, c=dwarfs_feh, marker='^', facecolor='none', label=r'Dwarfs ($\log{g} = 4.5$)')

scat_giants = ax.scatter(giants_gmi, giants_ewmg, s=70., cmap=cmap, c=giants_feh, marker='v', facecolor='none', label=r'Giants ($\log{g} = 2.0$)')


#observed data
ax.scatter(table.data['(g-i)'], table.data['EW_{MgH}'], s=7., c='k', marker='o', label='Observations', zorder=2)


ax.set_xlabel(r'$(g - i)_o$', fontsize=fontsize)
ax.set_ylabel(r'EW Mg I (\AA{})', fontsize=fontsize)

def make_axes(parent, **kw):
    orientation = kw.setdefault('orientation', 'vertical')
    fraction = kw.pop('fraction', 0.15)
    shrink = kw.pop('shrink', 1.0)
    aspect = kw.pop('aspect', 20)
    pb = parent.get_position(original=True).frozen()
    if orientation == 'vertical':
        location = kw.pop('location', 1)
        pad = kw.pop('pad', 0.05)
        if location:
            x1 = 1.0-fraction
            pb1, pbx, pbcb = pb.splitx(x1-pad, x1)
            pbcb = pbcb.shrunk(1.0, shrink).anchored('C', pbcb)
            anchor = (0.0, 0.5)
            panchor = (1.0, 0.5)
        else:
            pbcb, pbx, pb1 = pb.splitx(fraction, fraction+pad)
            pbcb = pbcb.shrunk(1.0, shrink).anchored('C', pbcb)
            anchor = (1.0, 0.5)
            panchor = (0.0, 0.5)
    else:
        location = kw.pop('location', 0)
        pad = kw.pop('pad', 0.15)
        if location:
            y1 = 1.0-fraction
            pb1, pbx, pbcb = pb.splity(y1-pad, y1)
            pbcb = pbcb.shrunk(shrink, 1.0).anchored('C', pbcb)
            anchor = (0.5, 0.0)
            panchor = (0.5, 1.0)
        else:
            pbcb, pbx, pb1 = pb.splity(fraction, fraction+pad)
            pbcb = pbcb.shrunk(shrink, 1.0).anchored('C', pbcb)
            anchor = (0.5, 1.0)
            panchor = (0.5, 0.0)
        aspect = 1.0/aspect
    parent.set_position(pb1)
    parent.set_anchor(panchor)
    fig = parent.get_figure()
    cax = fig.add_axes(pbcb)
    cax.set_aspect(aspect, anchor=anchor, adjustable='box')
    return cax, kw


#ax.legend(loc=2)
#cax,kw = make_axes(ax, orientation='horizontal', location=1.)

labelsize=14

#cbar = plt.colorbar(scat_giants, ax=ax, cax=cax, **kw)
#cbar.set_ticks(arange(min(dwarfs_feh), max(dwarfs_feh) + 0.5, 0.5))
#cbar.set_ticks([-2.5 + 1/3., -1.5, -0.5 - 1/3.])
#cbar.set_ticklabels(['-2.5', '-1.5', '-0.5'])
#cbar.ax.xaxis.set_tick_params(labelsize=labelsize)
ax.set_xlim(0.6001, 1.80001)
ax.set_ylim(0, 15.0001)
ax.set_xticks([0.8, 1.0, 1.2, 1.4, 1.6, 1.8])
ax.set_xticklabels(map(str, [0.8, 1.0, 1.2, 1.4, 1.6, 1.8]))
ax.set_yticks([0, 2, 4, 6, 8, 10, 12, 14])
ax.set_yticklabels(map(str, [0, 2, 4, 6, 8, 10, 12, 14]))
ax.xaxis.set_tick_params(labelsize=labelsize)
ax.yaxis.set_tick_params(labelsize=labelsize)

#cax.set_position([0.25, 0.76, 0.5, 0.8])
#cax.set_title('[Fe/H]', fontsize=fontsize)


#plt.draw()
plt.savefig('dwarfgiants.eps')
#plt.savefig('dwarfgiants.png')
