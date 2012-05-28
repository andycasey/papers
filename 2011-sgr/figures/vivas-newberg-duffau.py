

import atpy
import ephem
import matplotlib.pyplot as plt
fontsize=13
vivas = atpy.Table('Vivas2008.xml', verbose=False)
sirko = atpy.Table('Sirko2004a.vot', verbose=False)


RRLs = vivas.where((6.5 < vivas.Dis) & (vivas.Dis < 9.5) & (-80 < vivas.Vgsr) & (vivas.Vgsr < -10))
BHBs = sirko.where((6.5 < sirko.Dist) & (sirko.Dist < 10.5) & (-80 < sirko.Vgal) & (sirko.Vgal < -10) & (170 < sirko.RAJ2000) & (sirko.RAJ2000 < 215) & (-4 < sirko.DEJ2000) & (sirko.DEJ2000 < 4))

plt.close('all')
fig = plt.figure(figsize=(15.0, 4.8))
ax = fig.add_subplot(111)
plt.subplots_adjust(left=0.05, bottom=0.05, right=0.97, top=0.97, wspace=0.20, hspace=0.20)
plt.rcParams['text.usetex'] = True
#'font.size' : 10,
#          'axes.labelsize' : 10,
#          'font.size' : 10,
#          'text.fontsize' : 10,
#          'legend.fontsize': 10,
#          'xtick.labelsize' : 8,
#          'ytick.labelsize' : 8,
plt.rcParams['font.size'] = 11
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['text.fontsize'] = 12
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
#fig.set_figheight(6.8125)
#fig.set_figwidth(16.712499999999999)


# Newberg 2007 SEGUE plates


newbergPlates = [[300, 55], [288, 62]]
area = 7.0 # SEGUE overview document on ADS
radius = (area/pi)**0.5

for newbergPlate in newbergPlates:

    l, b = map(radians, newbergPlate)
     
    coords = ephem.Galactic(l, b, epoch='1958')
    ra, dec = map(degrees, coords.to_radec())
    
    circle = plt.Circle((ra, dec), radius=radius, ec='#000000', fc='none', linestyle='dashed')

    ax.add_patch(circle)
    
    
# BHB Sirko stars (used by Vivas)


RA_idx = BHBs.columns.keys.index('RAJ2000')
Dec_idx = BHBs.columns.keys.index('DEJ2000')

    
ax.scatter(BHBs.data['RAJ2000'], BHBs.data['DEJ2000'], marker='^', facecolor='none')

# RRL Vivas stars (used by Vivas, obviously
ax.scatter(RRLs.data['RAdeg'], RRLs.data['DEdeg'], marker='^', facecolor='none', label='Vivas et al. 2008, $-80 < V_{GSR} < -10$ $\mbox{km s}^{-1}$')
    
    
ax.set_xlim(165, 215)
ax.set_ylim(-4, 4)
ax.axis('scaled')

observations = atpy.Table('kgiants.xml', verbose=False)
subset1 = '(observations.K_Giants == True) & (-54. < observations.Vgsr) & (observations.Vgsr < -34.)'
subset2 = '(observations.K_Giants == True) & (-91. < observations.Vgsr) & (observations.Vgsr < -71.)'


# Generate the subset

subset1 = observations.where(eval(subset1))
subset2 = observations.where(eval(subset2))

ax.scatter(subset1.data['RA'], subset1.data['DEC'], marker='s', c='g', label=r'This work, $-54 < V_{GSR} < -34$ $\mbox{km s}^{-1}$')
ax.scatter(subset2.data['RA'], subset2.data['DEC'], marker='o', c='b', label=r'This work, $-91 < V_{GSR} < -71$ $\mbox{km s}^{-1}$')


# Our observed paltes
fields = [
            [180., 0.,],
            [185., -1.],
            [190., -2.],
            [194., -2.7]
            ]

for field in fields:
    ra, dec = field
    circle = plt.Circle((ra, dec), radius=1., ec='#000000', fc='none')
    ax.add_patch(circle)

# Duffau et al 2006 region
rect = plt.Rectangle((180, -4), 30, 2, linestyle='dashdot', ec="#000000", fill=False)
ax.add_patch(rect)




ax.set_xlabel(r'RA (deg)')
ax.set_ylabel(r'Dec (deg)')
ax.legend(loc=4)
plt.draw()

"""
fig = plt.figure()

subset3 = '(observations.K_Giants == True) & (-80. < observations.Vgsr) & (observations.Vgsr < -10.)'

subset3 = observations.where(eval(subset3)).data['[Fe/H] (Battaglia)']
subset3 = subset3[subset3<0]
subset3 = subset3[-3<subset3]

subset2 = subset2.data['[Fe/H] (Battaglia)']
subset2 = subset2[subset2<0]
subset2 = subset2[-3<subset2]

subset1 = subset1.data['[Fe/H] (Battaglia)']
subset1 = subset1[subset1<0]
subset1 = subset1[-3<subset1]

ax = fig.add_subplot(111)
ax.bar(arange(-2.9, -0.5, 0.3), subset1, color='r')
ax.bar(arange(-2.9, -0.5, 0.3), subset2, bottom=subset1, color='g')
#ax.hist((subset3, subset2, subset1), bins=arange(-2.9, -0.5, 0.3), facecolor=['w', 0.3, 0.6], rwidth=1.0, histtype='barstacked')
#ax.hist(subset3, bins=arange(-2.9, -0.5, 0.3), edgecolor='#000000', facecolor='none', linestyle='solid', label=r'$-80 < V_{GSR} < -10$ $\mbox{km s}^{-1}$')
#ax.hist(subset2, bins=arange(-2.9, -0.5, 0.3), edgecolor='#000000', facecolor='none', linestyle='dashed', label=r'$-86 < V_{GSR} < -66$ $\mbox{km s}^{-1}$')
#ax.hist(subset1, bins=arange(-2.9, -0.5, 0.3), edgecolor='#000000', facecolor='none', linestyle='dotted', label=r'$-59 < V_{GSR} < -39$ $\mbox{km s}^{-1}$')
#ax.legend()
    
   
"""

fig.text(0.305, 0.555, 'Field 1')#, fontsize=fontsize)
fig.text(0.38, 0.52, 'Field 2')#, fontsize=fontsize)
fig.text(0.455, 0.47, 'Field 3')#, fontsize=fontsize)
fig.text(0.52, 0.45, 'Field 4')#, fontsize=fontsize)
fig.text(0.47, 0.625, 'Newberg et al. 2007 plate at $(l, b) = (288^\circ, 62^\circ)$')#, fontsize=fontsize)
fig.text(0.245, 0.29, 'Newberg et al. 2007 plate at $(l, b) = (300^\circ, 55^\circ)$')#, fontsize=fontsize)
fig.text(0.70, 0.52, 'Duffau et al. 2006 region')
plt.draw()
plt.savefig('vivas-newberg-duffau.eps')
