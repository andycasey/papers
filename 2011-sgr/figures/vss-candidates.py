import matplotlib.pyplot as plt
import matplotlib.patches as patches
import atpy
from numpy import *
import matplotlib
import matplotlib.text as text
from pylab import *
from matplotlib.image import NonUniformImage, AxesImage
from scipy import interpolate
plt.rcParams['text.usetex'] = True
plt.rcParams['text.fontsize'] = 13.
plt.rcParams['font.size'] = 13.
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

# standard binning
xbound = [-300., 300]
xerr = 15.
yerr = 0.28
ybound = [-2.46, -0.3]

# offset by half
#xbound = [-307.5, 300]
#xerr = 15.
#yerr = 0.30
#ybound = [-2.60, -0.5]

vgsr = []
feh = []
feherr = []
vgsrerr = []

for item in table.data:
    if item['kgiants'] and (item['V_{gsr} (fxcor)'] > xbound[0]) and (item['V_{gsr} (fxcor)'] < xbound[1]) and (item['[Fe/H] (Battaglia)'] > ybound[0]) and (ybound[1] > item['[Fe/H] (Battaglia)']):
        feh.append(item['[Fe/H] (Battaglia)'])
        feherr.append(item['[Fe/H]_{err} (Battaglia)'])
        vgsr.append(item['V_{gsr} (fxcor)'])
        vgsrerr.append(item['V_{err} (fxcor)'])
        

plt.close('all')
myfig = plt.figure()

myax = myfig.add_subplot(121)

myax.errorbar(vgsr, feh, xerr=vgsrerr, yerr=feherr, fmt='o', color='k', markersize=4.5)
myax.set_xlim(xbound)
myax.set_ylim(ybound)

myax.set_xlabel(r'$V_{GSR}$ (km s$^{-1}$)')
myax.set_ylabel(r'[Fe/H]')
myax.set_ylim(ybound[0], ybound[1])
myax.set_xticks([-300, -200, -100, 0, 100, 200,300])
myax.set_xticklabels(map(str, ['', -200, -100, 0, 100, 200, '']))
myax.set_yticks([-2.5, -2.0, -1.5, -1.0, -0.5])
myax.set_yticklabels(map(str, ['', -2.0, -1.5, -1.0, -0.5]))


#ellipse = patches.Ellipse([130., -1.95], 40., 0.3, linestyle='dashed', facecolor='none', edgecolor='k')
#myax.add_artist(ellipse)

#text = myax.text(170., -2.2, "VSS K-Giant \nmembers",verticalalignment='top')
#myax.add_artist(text)


#rect = Rectangle((120, myax.get_ylim()[0]), 35, 4, facecolor='#dddddd', edgecolor='none')
#myax.add_patch(rect)

myax.set_xlim(-299, 299)
#arrow = patches.FancyArrowPatch([170., -2.2], [130., -2.1], arrowstyle='->', mutation_scale=15.)
#myax.add_patch(arrow)

#ax.errorbar(-270., -2.25, yerr=0.14, color='k')

plt.draw()
#fig.savefig('vss-candidates.eps')

x = arange(xbound[0], xbound[1], xerr)
y = arange(ybound[0], ybound[1], yerr)

bins, xdata, ydata = histogram2d(vgsr, feh, bins=[x, y])

bins = flipud(transpose(bins))
'''
fig = plt.figure()
ax = fig.add_subplot(111)

cmap = cmap_discretize(matplotlib.cm.RdBu_r, 1 + max(max(i) for i in bins))


#im = AxesImage(ax, interpolation='nearest', cmap=cmap, extent=(xbound[0], xbound[1], ybound[0], ybound[1]))


# Monte Carlo analysis


im = imshow(bins, interpolation='nearest', cmap=cmap, extent=(xdata[0], xdata[-1], ydata[0], ydata[-1]), aspect='auto')
ax.grid(True)

#ticks = list(ax.get_xticks()[1::]*10. -300)
#ticks.insert(0, ' ')

#ax.set_xticklabels(ticks)

#ticks = -ax.get_yticks()/max(ax.get_yticks()) * (ybound[1]-ybound[0])
#im.set_data(xdata[0:-1] + xerr/2., ydata[0:-1] + yerr/2., bins)
#im.set_data(bins)
#ax.images.append(im)
ax.set_xlim(xdata[0], xdata[-1])
ax.set_ylim(ydata[0], ydata[-1])

maxval = max(max(i) for i in bins)
cbar = plt.colorbar(im)
#cbar.set_clim(0, 1 + max(max(i) for i in bins))
cbar.set_ticks((arange(0, maxval + 1) + 0.5)/((maxval + 1)/maxval))


#cbar.set_ticklabels(map(int, arange(0, max(max(i) for i in bins))))
cbar.set_ticklabels(map(int, arange(0, maxval + 1)))

ax.scatter(vgsr, feh, marker='o', color='k', s=3)
'''

#plt.savefig('mc-grid2.eps')

print "RUNNING MONTE CARLO SIMULATION"
import random

x = array([-4.  , -3.4 , -3.  , -2.5 , -2.  , -1.7 , -1.5 , -1.1 , -0.89, -0.5 ,  0.  ,  0.5 ])
y = array([0, 25, 30, 90, 240, 270, 290, 307, 310, 390, 430, 540])

f = interpolate.splrep(y, x)

number_occurred = zeros((bins.shape))

for i in range(0, 10**4):
    print "Simulation #%s" % i
    sim_vgsr = zeros(len(vgsr))
    sim_feh = []
    
    
    for j in range(0, len(vgsr)):
        sim_vgsr[j] = random.gauss(0., 101.6) #Sirko et al 2004

    
    for j in range(0, round(len(vgsr)*0.02)):
        sim_feh.append(random.gauss(-0.89, 0.3))
        
    for j in range(0, round(len(vgsr)*0.12)):
        sim_feh.append(random.gauss(-0.5, 0.4))

    num = len(sim_vgsr) - len(sim_feh)
    for j in range(0, num):
        sim_feh.append(random.gauss(-1.7, 0.7))

    #raise
    sim_bins, sim_xdata, sim_ydata = histogram2d(sim_vgsr, np.array(sim_feh), bins=[xdata,ydata])
    sim_bins = transpose(sim_bins)
    
    difference = sim_bins - bins
    

    for k, v in enumerate(difference):
        for l, w in enumerate(v):
            if w < 0.:
                difference[k][l] = 0.
            
            elif w > 0.:
                difference[k][l] = 1.
                
    number_occurred += difference
                
        


number_occurred = log10(number_occurred + 1)
textsize = 15
labelsize = 14

#fig = figure()
#fig.subplots_adjust(left=, right=, top=, bottom=)

ax = myfig.add_subplot(122)
im = ax.imshow(number_occurred, interpolation='nearest', cmap=matplotlib.cm.binary_r, extent=(sim_xdata[0], sim_xdata[-1], sim_ydata[0], sim_ydata[-1]), aspect='auto')
ax.grid(True)
ax.set_xlim(xbound[0], xbound[1])
ax.set_ylim(ybound[0], ybound[1])
ax.set_xticks([-300, -200, -100, 0, 100, 200,300])
ax.set_xticklabels(map(str, ['', -200, -100, 0, 100, 200, '']))
ax.set_yticks([-2.5, -2.0, -1.5, -1.0, -0.5])
ax.set_yticklabels(map(str, ['', -2.0, -1.5, -1.0, -0.5]))
ax.set_yticklabels(['']*5)

#ax.yaxis.set_tick_params(labelsize=labelsize)
ax.xaxis.set_tick_params(labelsize=labelsize)
#ax.set_xlabel(r'$V_{GSR}$ (km s$^{-1}$)', fontsize=textsize)
#ax.set_xlabel(r'$\textrm{V}_{\textrm{\textsc{GSR}}}\;\;(\textrm{km s}^{-1})$', fontsize=textsize)
ax.set_xlabel(r'$V_{GSR}$ (km s$^{-1}$)')
#ax.set_ylabel(r'[Fe/H]', fontsize=textsize)
#ax.errorbar(vgsr, feh, xerr=vgsrerr, yerr=feherr, fmt='o', color='k', markersize=4.5)
ax.set_xlim(sim_xdata[1], sim_xdata[-1])
ax.set_ylim(sim_ydata[0], sim_ydata[-1])


cbar = plt.colorbar(im)
cbar.set_label(r'$\log{\left(1+N\right)}$', fontsize=textsize)
cbar.ax.set_yticks([0., 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6])
cbar.ax.set_yticklabels(map(str, [0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8, 3.2, 3.6]))
cbar.ax.yaxis.set_tick_params(labelsize=labelsize)



ellipse = patches.Ellipse([130., -2.05], 40., 0.3, linewidth=2, linestyle='solid', facecolor='none', edgecolor='k', zorder=100)
ax.add_artist(ellipse)
ellipse = patches.Ellipse((-90, -1.7), 130., 0.8, linewidth=2, linestyle='solid', facecolor='none', edgecolor='k', zorder=100)
ax.add_artist(ellipse)

ellipse = patches.Ellipse((220., -1.3), 50., 0.7, linewidth=2, linestyle='solid', facecolor='none', edgecolor='k', zorder=100)
ax.add_artist(ellipse)


ax.text(-252, -1.25, 'Feature A')
ax.text(140, -2.3, 'Feature B')
ax.text(160, -0.85, 'Feature C')

myax.set_ylim(-2.45, -0.5)
ax.set_ylim(-2.45, -0.5)
#myax.set_ylim(sim_ydata[0], sim_ydata[-1])#ybound[0], ybound[1])
#myax.set_yticks([-2.5, -2.0, -1.5, -1.0, -0.5])
#myax.set_yticklabels(map(str, [-2.5, -2.0, -1.5, -1.0, -0.5]))

plt.draw()

#plt.savefig('montecarlo.pdf')
#plt.savefig('montecarlo.png')
#plt.savefig('montecarlo.eps')
'''
        
    
fig = figure()
mins = [min(_) for _ in transpose(number_occurred)]
xdata_ = xdata[0:-1] + (xdata[1] - xdata[0])/2. #assumes equal spacing

plot(xdata_, mins)
xlabel(r'$V_{GSR}$ (km s$^{-1}$)')
ylabel(r'$\log{1+N}$')
#plt.savefig('num_occurred.eps')
'''
