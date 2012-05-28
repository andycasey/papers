

import matplotlib.pyplot as plt

import atpy
from numpy import *


column = 'V_{gsr} (fxcor)'
table = atpy.Table('multiple_observations.xml', verbose=False)



column = 'Derived error ' + column

plt.close('all')
fig = plt.figure()
"""
ax1 = fig.add_subplot(121)

binsize = 1
groupSizes = table.data['GroupSize']
for i, groupSize in enumerate(table.data['GroupSize']):
    if groupSize < 0:
        groupSizes[i] = 1
    else: groupSizes[i] = int(groupSize)
        
ax1.hist(groupSizes, bins=max(groupSizes), range=[1, max(groupSizes) + 1], align='mid', rwidth=0.90, color='#666666')
ax1.set_xlabel(r'Number of repeated observations')
ax1.set_ylabel(r'Number of stars')
    
ax1.set_xticklabels(['', '1', '', '2', '', '3', '', '4', ''])

"""

ax = fig.add_subplot(111)
binsize = 2.5 # km s
range = [0, 20]

bins = range[1]/binsize
bins = ax.hist(table.data[column], bins=bins, range=range, align='mid', rwidth=0.90, color='#666666')
ax.set_xlabel(r'Internal velocity error (km s$^{-1}$)')
ax.set_ylabel('Number of observations')

x = arange(0, 20, 0.01)
variance = std(table.data[column])

y = 1/sqrt(2*pi*variance**2) * exp((-x**2) / (2*variance**2))
 
y = bins[0][0] * y / max(y)
ax.plot(x, y, 'k--')


from matplotlib.patches import FancyArrowPatch

for i, _ in enumerate(y):
    if bins[0][0]/2. > _:
        break

i = i - 1

    
arrow = FancyArrowPatch([0, bins[0][0]/2.], [x[i], bins[0][0]/2.], arrowstyle='<|-|>', mutation_scale=15., facecolor='k', edgecolor='k')
ax.add_patch(arrow)

plt.figtext((x[i]*2.)/ax.get_xbound()[1] - 0.07, (bins[0][0]/2.)/ax.get_ybound()[1], r'HWHM = %3.2f km s$^{-1}$' % x[i], color='k', horizontalalignment='left')

plt.draw()

fig.savefig('vel_error.eps')