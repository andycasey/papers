#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import atpy
import pdb
from numpy import *

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib.patches as mpatches
plt.rcParams['text.usetex'] = True

modelFiles = ['prolate.dat', 'oblate.dat', 'spherical.dat', 'SgrTriax_DYN.dat']
#modelFiles = ['SgrTriax_DYN.dat']
textsize = 15
labelsize = 14
observations = atpy.Table('combined.xml', verbose=False)

observedSpatialLimits = ['Sgr_lambda', 'Sgr_b']
modelSpatialLimits = ['lambda', 'beta']
spatialTolerances = [0.0, 0.0]
xbounds = [250, 345]
ybounds = [-30, 30]
distbounds = [0, 60]
inc = 15

observedDataParameter = 'Vgsr'
observedErrorParameter = 'Verr'
modelParameter = 'vgsr'
modelDistanceParameter = 'dist'
modelMarkerParameter = 'Pcol'

subsetLogic = '(observations.kgiants == True) & (-250. < observations.Vgsr) & (observations.Vgsr < 250.)'

# Number of data points in x-dimension
resolution = 5000


# Generate the subset

observationsSubset = observations.where(eval(subsetLogic)).data
observedData = observationsSubset[observedDataParameter]
observedError = observationsSubset[observedErrorParameter]

# Calculate bandwidth and ranges

# Bandwidth estimated using a Epanechnikov kernel
epanechnikov = ((40*(np.pi)**0.5)/len(observationsSubset[observedDataParameter]))**(0.2)

bandwidth = 10.0

offset = 100.0
x = np.linspace(min(observedData) - offset, offset + max(observedData), resolution)
#x = np.linspace(min(xBound), max(xBound), resolution)
#print "min, max", min(observedData), max(observedData)
#pdb.set_trace()

observedHistogram = np.zeros(len(x))

spatialLimits = []
for observedSpatialLimit in observedSpatialLimits:
    spatialLimits.append([min(observationsSubset[observedSpatialLimit]), max(observationsSubset[observedSpatialLimit])])

print spatialLimits

# Generate the observed generalized histogram

for i, (observationMean, observationalError) in enumerate(zip(observationsSubset[observedDataParameter], observationsSubset[observedErrorParameter])):
    observedHistogram += exp((-(x-observationMean)**2)/(2*bandwidth**2))/(bandwidth * sqrt(2*np.pi))
    
# Normalise and smooth
observedHistogram /= len(observedHistogram)

# Generate the model generalized histograms
modelHeaders =  {
                    'SgrTriax_DYN.dat' : 'lambda beta ra dec l b xgc ygc zgc xsun ysun zsun x4 y4 z4 u v w dist vgsr mul mub mua mud Pcol Lmflag',
                    'oblate.dat'    :'lambda beta l b xgc ygc zgc xsun ysun zsun u v w dist vgsr Pcol',
                    'prolate.dat'   :'lambda beta l b xgc ygc zgc xsun ysun zsun u v w dist vgsr Pcol',
                    'spherical.dat' :'lambda beta l b xgc ygc zgc xsun ysun zsun u v w dist vgsr Pcol'
                }

colours = ['#F2BA52', '#CE040F', '#CE040F', '#2165BF', '#2165BF', 'green', 'green']


xBound = [min(observedData), max(observedData)]
    
plt.close('all')

field_rad = (2./3.14159)**0.5
fields = [(257.25, 16.30), (262.32, 14.97), (267.36, 13.50), (271.3, 12.21)]

field_xlims = [min([field[0] for field in fields]) - field_rad, max([field[0] for field in fields]) + field_rad]
field_ylims = [min([field[1] for field in fields]) - field_rad, max([field[1] for field in fields]) + field_rad]



#fig, (triaxial, prolate, spherical, oblate) = plt.subplots(4, 1, sharex=True, sharey=False)
fig = plt.figure(figsize=(4.5,9))
fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.1, wspace=0.05)
triaxial = plt.subplot(411)
prolate = plt.subplot(412)
spherical = plt.subplot(413)
oblate = plt.subplot(414)

model = loadtxt('SgrTriax_DYN.dat', usecols=(0, 1, 25, 18))
joined = unique([' '.join(map(str, [round(j, 2) for j in i])) for i in [_ for _ in model]])
#joined = [' '.join(map(str, _)) for _ in model]
joined = filter(lambda x: 3 >= float(x.split()[2]) and float(x.split()[2]) >= 0, joined)
joined = filter(lambda x: xbounds[1] >= float(x.split()[0]) and float(x.split()[0]) >= xbounds[0], joined) # second x.plit() was [1] instead of [0]
joined = filter(lambda x: ybounds[1] >= float(x.split()[1]) and float(x.split()[1]) >= ybounds[0], joined)
joined = filter(lambda x: distbounds[1] >= float(x.split()[3]) and float(x.split()[3]) >= distbounds[0], joined)


model = [map(float, _.split()) for _ in joined]

for i, point in enumerate(model):

    triaxial.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))

triaxial.set_xbound(xbounds)
triaxial.set_xticklabels([], alpha=0.0)
triaxial.set_ybound(ybounds)
triaxial.set_yticks(np.arange(ybounds[0], ybounds[1] + inc, inc))
triaxial.set_yticklabels(map(str, np.arange(ybounds[0], ybounds[1] + inc, inc)))
triaxial.yaxis.set_tick_params(labelsize=labelsize)


[triaxial.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]

axins = zoomed_inset_axes(triaxial, 4, loc=5) # zoom = 6

for i, point in enumerate(model):
    axins.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))

[axins.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]

axins.xaxis.set_ticks([255,260, 265, 270])
axins.yaxis.set_ticks([12, 14, 16])
axins.xaxis.set_ticklabels(map(str, [255, 260, 265, 270]))
axins.yaxis.set_ticklabels(map(str, [12, 14, 16]))
axins.xaxis.set_tick_params(labelsize=labelsize-1)
axins.yaxis.set_tick_params(labelsize=labelsize-1)

axins.set_xlim(field_xlims)
axins.set_ylim(field_ylims)

mark_inset(triaxial, axins, loc1=1, loc2=3, fc='none', ec="#000000", zorder=200)
#raise



model = loadtxt('prolate.dat', usecols=(0, 1, 15, 13))
joined = unique([' '.join(map(str, [float(j) for j in i])) for i in [_ for _ in model]])

joined = filter(lambda x: 3 >= float(x.split()[2]) and float(x.split()[2]) >= 0, joined)
joined = filter(lambda x: xbounds[1] >= float(x.split()[0]) and float(x.split()[0]) >= xbounds[0], joined)
joined = filter(lambda x: ybounds[1] >= float(x.split()[1]) and float(x.split()[1]) >= ybounds[0], joined)
joined = filter(lambda x: distbounds[1] >= float(x.split()[3]) and float(x.split()[3]) >= distbounds[0], joined)

model = [map(float, _.split()) for _ in joined]


for i, point in enumerate(model):
    
    prolate.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))

prolate.set_xbound(xbounds)
prolate.set_xticklabels([], alpha=0.0)
prolate.set_ybound(ybounds)
prolate.set_yticks(np.arange(ybounds[0], ybounds[1] + inc, inc))
prolate.set_yticklabels(map(str, np.arange(ybounds[0], ybounds[1] + inc, inc)))
prolate.yaxis.set_tick_params(labelsize=labelsize)


[prolate.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]

axins = zoomed_inset_axes(prolate, 4, loc=5) # zoom = 6

for i, point in enumerate(model):
    axins.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))

[axins.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]


axins.xaxis.set_ticks([255,260, 265, 270])
axins.yaxis.set_ticks([12, 14, 16])
axins.xaxis.set_ticklabels(map(str, [255, 260, 265, 270]))
axins.yaxis.set_ticklabels(map(str, [12, 14, 16]))
axins.xaxis.set_tick_params(labelsize=labelsize-1)
axins.yaxis.set_tick_params(labelsize=labelsize-1)


axins.set_xlim(field_xlims)
axins.set_ylim(field_ylims)

#axins.xaxis.set_visible(False)
#axins.yaxis.set_visible(False)

mark_inset(prolate, axins, loc1=1, loc2=3, fc='none', ec="#000000", zorder=200)

        

model = loadtxt('spherical.dat', usecols=(0, 1, 15, 13))
joined = unique([' '.join(map(str, [float(j) for j in i])) for i in [_ for _ in model]])
#joined = [' '.join(map(str, _)) for _ in model]
joined = filter(lambda x: 3 >= float(x.split()[2]) and float(x.split()[2]) >= 0, joined)
joined = filter(lambda x: xbounds[1] >= float(x.split()[0]) and float(x.split()[0]) >= xbounds[0], joined)
joined = filter(lambda x: ybounds[1] >= float(x.split()[1]) and float(x.split()[1]) >= ybounds[0], joined)
joined = filter(lambda x: distbounds[1] >= float(x.split()[3]) and float(x.split()[3]) >= distbounds[0], joined)

model = [map(float, _.split()) for _ in joined]

for i, point in enumerate(model):

    spherical.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))


spherical.set_xbound(xbounds)
spherical.set_xticklabels([], alpha=0.0)
spherical.set_ybound(ybounds)
spherical.set_yticks(np.arange(ybounds[0], ybounds[1] + inc, inc))
spherical.set_yticklabels(map(str, np.arange(ybounds[0], ybounds[1] + inc, inc)))
spherical.yaxis.set_tick_params(labelsize=labelsize)



[spherical.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]

axins = zoomed_inset_axes(spherical, 4, loc=5) # zoom = 6

for i, point in enumerate(model):
    axins.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))

[axins.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]


axins.xaxis.set_ticks([255,260, 265, 270])
axins.yaxis.set_ticks([12, 14, 16])
axins.xaxis.set_ticklabels(map(str, [255, 260, 265, 270]))
axins.yaxis.set_ticklabels(map(str, [12, 14, 16]))
axins.xaxis.set_tick_params(labelsize=labelsize-1)
axins.yaxis.set_tick_params(labelsize=labelsize-1)


axins.set_xlim(field_xlims)
axins.set_ylim(field_ylims)


#axins.xaxis.set_visible(False)
#axins.yaxis.set_visible(False)

mark_inset(spherical, axins, loc1=1, loc2=3, fc='none', ec="#000000", zorder=200)




model = loadtxt('oblate.dat', usecols=(0, 1, 15, 13))
joined = unique([' '.join(map(str, [float(j) for j in i])) for i in [_ for _ in model]])
#joined = [' '.join(map(str, _)) for _ in model]
joined = filter(lambda x: 3 >= float(x.split()[2]) and float(x.split()[2]) >= 0, joined)
joined = filter(lambda x: xbounds[1] >= float(x.split()[0]) and float(x.split()[0]) >= xbounds[0], joined)
joined = filter(lambda x: ybounds[1] >= float(x.split()[1]) and float(x.split()[1]) >= ybounds[0], joined)
joined = filter(lambda x: distbounds[1] >= float(x.split()[3]) and float(x.split()[3]) >= distbounds[0], joined)

model = [map(float, _.split()) for _ in joined]

for i, point in enumerate(model):

    oblate.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))


oblate.set_xbound(xbounds)
oblate.set_xticks(np.arange(260, 360, 20))
oblate.set_xticklabels(map(str, np.arange(260, 360, 20)))

oblate.set_ybound(ybounds)
oblate.set_yticks(np.arange(ybounds[0], ybounds[1] + inc, inc))
oblate.set_yticklabels(map(str, np.arange(ybounds[0], ybounds[1] + inc, inc)))
oblate.yaxis.set_tick_params(labelsize=labelsize)


[oblate.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]

axins = zoomed_inset_axes(oblate, 4, loc=5) # zoom = 6

for i, point in enumerate(model):
    axins.scatter(point[0], point[1], s=1, color=colours[int(point[2])], zorder=(100 - int(point[3])))

[axins.add_patch(mpatches.Circle(field, field_rad, edgecolor='#000000', facecolor='none', zorder=150)) for field in fields]

axins.xaxis.set_ticks([255,260, 265, 270])
axins.yaxis.set_ticks([12, 14, 16])
axins.xaxis.set_ticklabels(map(str, [255, 260, 265, 270]))
axins.yaxis.set_ticklabels(map(str, [12, 14, 16]))
axins.xaxis.set_tick_params(labelsize=labelsize-1)
axins.yaxis.set_tick_params(labelsize=labelsize-1)


axins.set_xlim(field_xlims)
axins.set_ylim(field_ylims)


#axins.xaxis.set_visible(False)
#axins.yaxis.set_visible(False)

mark_inset(oblate, axins, loc1=1, loc2=3, fc='none', ec="#000000", zorder=200)

        
    
oblate.set_xlabel(r'$\Lambda_\odot$ (deg)', fontsize=textsize)
oblate.set_ylabel(r'$B_\odot$ (deg)',fontsize=textsize)
spherical.set_ylabel(r'$B_\odot$ (deg)',fontsize=textsize)
triaxial.set_ylabel(r'$B_\odot$ (deg)',fontsize=textsize)
prolate.set_ylabel(r'$B_\odot$ (deg)',fontsize=textsize)

xlim, ylim = spatialLimits




#rect = mpatches.Rectangle((xlim[0], ylim[0]), xlim[1] - xlim[0], ylim[1] - ylim[0], edgecolor='#000000', facecolor='none', zorder=100)

#[triaxial.add_patch(fieldN) for fieldN in [fieldA, fieldB, fieldC, fieldD]]
#[prolate.add_patch(fieldN) for fieldN in [fieldA, fieldB, fieldC, fieldD]]
#[spherical.add_patch(fieldN) for fieldN in [fieldA, fieldB, fieldC, fieldD]]
#[oblate.add_patch(fieldN) for fieldN in [fieldA, fieldB, fieldC, fieldD]]

triaxial.text(339, 21, r'Tri-axial (LM10)', horizontalalignment='right',fontsize=labelsize, zorder=150)
prolate.text(339, 21, r'Prolate (LJM05)', horizontalalignment='right',fontsize=labelsize,zorder=150)
spherical.text(339, 21, r'Spherical (LJM05)', horizontalalignment='right', fontsize=labelsize,zorder=150)
oblate.text(339, 21, r'Oblate (LJM05)', horizontalalignment='right',fontsize=labelsize,zorder=150)

plt.draw()
import os
os.system('rm -f law_test_spatial.*')
plt.savefig('law_test_spatial.pdf')
plt.savefig('law_test_spatial.png')

"""
triaxial.add_patch(fieldA)
triaxial.add_patch(fieldB)
triaxial.add_patch(fieldC)
triaxial.add_patch(fieldD)



prolate.add_patch(fieldA)
prolate.add_patch(fieldB)
prolate.add_patch(fieldC)
prolate.add_patch(fieldD)

spherical.add_patch(fieldA)
spherical.add_patch(fieldB)
spherical.add_patch(fieldC)
spherical.add_patch(fieldD)

oblate.add_patch(fieldA)
oblate.add_patch(fieldB)
oblate.add_patch(fieldC)
oblate.add_patch(fieldD)
"""


