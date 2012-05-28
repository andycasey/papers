

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import Image
from numpy import linspace, array
plt.rcParams['text.usetex'] = True
plot = Image.open('carbonstars.jpg')
dpi = plt.rcParams['figure.dpi']
figsize = plot.size[0]/float(dpi), plot.size[1]/float(dpi)


textsize = 16
labelsize = 15

xnum, ynum = 6, 4 #original bounds of the image
# flipping it around
xnum, ynum = 4, 6

plt.close('all')
fig = plt.figure(figsize=figsize)
fig.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95)

subplot = fig.add_subplot(111)
#subplot.set_frame_on(False)
#subplot.set_axis_off()
subplot.imshow(plot, origin='lower', aspect='auto')


xbound, ybound = subplot.get_xbound(), subplot.get_ybound()
xscale, yscale = (xbound[1]-xbound[0])/float(xnum), (ybound[1]-ybound[0])/float(ynum)

subplot.set_xticks(linspace(xbound[0], xbound[1], xnum + 1))
labels = ['']
labels.extend(range(1, xnum+1))

subplot.set_xticklabels(map(str, labels))

subplot.xaxis.set_tick_params(labelsize=labelsize)


subplot.set_yticks(linspace(ybound[0], ybound[1], ynum+1))
labels = ['']
labels.extend(range(1, ynum+1))
subplot.set_yticklabels(map(str, labels))
subplot.yaxis.set_tick_params(labelsize=labelsize)
#subplot.set_xlabel(r'\textit{u} $-$ \textit{g}', fontsize=textsize)
#subplot.set_ylabel(r'\textit{g} $-$ \textit{r}', fontsize=textsize)
# flipping it around
subplot.set_xlabel(r'$(g - r)_o$', fontsize=textsize)
subplot.set_ylabel(r'$(u - g)_o$', fontsize=textsize)


#Add our carbon stars

stars = [
            [1.096, 0.538],
            [1.668, 0.725],
            [1.789, 0.799],
            [1.428, 0.641],
            [1.719, 0.728]
        ]

# flipping it around
stars = [[b,a] for a, b in stars]

stars = array(stars).transpose()

subplot.scatter(stars[0]*xscale, stars[1]*yscale, s=50, c='k', marker='o', edgecolor='k', facecolor='orange', label='Carbon star contaminants (this work)')
markers = []
alpha = 1.0
markers.append(subplot.scatter(0, 0, s=50, c='k', alpha=alpha, marker='x', facecolor='none', label='SDSS FHLC dwarfs (Downes et al 2004)'))
markers.append(subplot.scatter(0, 0, s=50, c='k', alpha=alpha, marker='s', facecolor='none', label='SDSS FHLC giants (Downes et al 2004)'))
markers.append(subplot.scatter(0, 0, s=50, c='k', alpha=alpha, marker='^', facecolor='none', label='SDSS FHLCs - F/G type (Downes et al 2004)'))
markers.append(subplot.scatter(0, 0, s=50, c='k', alpha=alpha, marker='d', facecolor='none', label='SDSS FHLCs - uncertain (Downes et al 2004)'))
subplot.legend(loc=1)
plt.draw()
# hide markers

for marker in markers:
    marker.set_visible(False)
    
plt.draw()
fig.savefig('carbonstarcolor.eps')
fig.savefig('carbonstarcolor.pdf')
