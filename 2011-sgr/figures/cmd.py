#!/usr/bin/python


from matplotlib import pyplot as plt
from matplotlib import rc as rc
from matplotlib.patches import Polygon
import atpy, glob
from numpy import loadtxt, log10

g_o = []
gmi_o = []


data = []

SDSS_files = glob.glob('SDSS_*.xml')

for SDSS_file in SDSS_files:
    table = atpy.Table(SDSS_file, verbose=False)
    
    # g_o
    [g_o.append(item) for item in table.data['g']]
    
    # (g - i)_o
    colour = table.data['g'] - table.data['i']
    [gmi_o.append(item) for item in colour]
    

# Plot it up
fig = plt.figure()
fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)
ax = fig.add_subplot(111)

fontsize = 14
labelsize = 14.5
ax.scatter(gmi_o, g_o, color='k', s=1, alpha=0.10)

# load up the isochrone data
isochrone_m1p5 = loadtxt('dartmouth_12gr_fehm15_alpha0.isochrone', usecols=(6,8))
isochrone_m0p5 = loadtxt('dartmouth_12gr_fehm05_alpha0.isochrone', usecols=(6,8))

d = 40 # kpc
gmi = isochrone_m0p5[:, 0] - isochrone_m0p5[:, 1]
g = -5 + 5*log10(d*10**3) + isochrone_m0p5[:,0]

ax.scatter(gmi, g, marker='o', edgecolor='#CE040F', facecolor='None', s=20, zorder=5, label='Dartmouth isochrone (12 Gyr, $r_\odot = 40$ kpc, [Fe/H] = -0.5)')
#ax.text(1.5, 16, '12 Gyr Isochrone at $r_\odot = 45$ kpc with [Fe/H] $= -1.5$', color='g', horizontalalignment='left', fontsize=labelsize)

gmi = isochrone_m1p5[:, 0] - isochrone_m1p5[:, 1]
g = - 5 + 5*log10(d*10**3) + isochrone_m1p5[:,0]

ax.scatter(gmi, g, marker='o', edgecolor='#2165BF', facecolor='None', s=20, zorder=5, label='Dartmouth isochrone (12 Gyr, $r_\odot = 40$ kpc, [Fe/H] = -1.5)')
#ax.text(1.75, 18., '12 Gyr Isochrone at $r_\odot = 45$ kpc with [Fe/H] $= -0.5$', color='b', horizontalalignment='left', fontsize=labelsize)


# Set the bounds
ax.set_xbound([-0.5, 3.5])
ax.set_ybound([13, 21])


# Reverse the y-axis
ax.set_ylim(ax.get_ylim()[::-1])

# Label axes
rc('text', usetex=True)


ax.yaxis.set_label_text(r'\textit{g}', fontsize=labelsize)
ax.xaxis.set_label_text(r'\textit{g} $-$ \textit{i}', fontsize=labelsize)
ax.set_yticks([21, 20, 19, 18, 17, 16, 15, 14, 13])                                                                                                                                                                                 
ax.set_yticklabels(map(str, [21, 20, 19, 18, 17, 16, 15, 14, 13]))
ax.set_xticks([-0.5, 0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
ax.set_xticklabels(map(str, ['', 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]))
ax.xaxis.set_tick_params(labelsize=labelsize)
ax.yaxis.set_tick_params(labelsize=labelsize)


# Draw selection box

# 0.6 < (g-i)_o < 1.7
# 15 < i_o < 18
# -15(g-i) + 27 < g_o < -3.75(g-i)_o

# This was drawn in TOPCAT and the following x, y bound points were found for a polygon that would
# describe these equation criterion.

pol = Polygon([[0.8, 15.8], [1.6, 16.8], [1.0, 19.1], [0.6, 18.7], [0.6, 17.9]], closed=True, facecolor='none', fill=False, linewidth=2)
ax.add_patch(pol)
ax.legend()
plt.draw()
#plt.savefig('cmd.pdf')
#plt.savefig('cmd.png')

