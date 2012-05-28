#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import atpy
from numpy import *
from matplotlib.patches import Rectangle

useColor = False

modelFiles = ['prolate.dat', 'oblate.dat', 'spherical.dat', 'SgrTriax_DYN.dat']

labelsize=10
fontsize=11

observations = atpy.Table('combined.xml', verbose=False)

observedSpatialLimits = ['Sgr_lambda', 'Sgr_b']
modelSpatialLimits = ['lambda', 'beta']
spatialTolerances = [0.0, 0.0]

observedDataParameter = 'Vgsr'
observedErrorParameter = 'Verr'
modelParameter = 'vgsr'
modelDistanceParameter = 'dist'
modelMarkerParameter = 'Pcol'

subsetLogic = '(observations.kgiants == True) & (-250. < observations.Vgsr) & (observations.Vgsr < 250.)'

# Number of data points in x-dimension
resolution = 500

xBound = [-400, 400]
distBound = [0, 60]

# Generate the subset

observationsSubset = observations.where(eval(subsetLogic)).data
observedData = observationsSubset[observedDataParameter]
observedError = observationsSubset[observedErrorParameter]

# Calculate bandwidth and ranges

# Bandwidth estimated using a Epanechnikov kernel
epanechnikov = ((40*(np.pi)**0.5)/len(observationsSubset[observedDataParameter]))**(0.2)

bandwidth = 10.0
#bandwidth = np.mean(observationsSubset[observedErrorParameter]) * epanechnikov

offset = 100.0
x = np.linspace(min(observedData) - offset, offset + max(observedData), resolution)
x = np.linspace(min(xBound), max(xBound), resolution)


#pdb.set_trace()

observedHistogram = zeros(len(x))

spatialLimits = []
for observedSpatialLimit in observedSpatialLimits:
    spatialLimits.append([min(observationsSubset[observedSpatialLimit]), max(observationsSubset[observedSpatialLimit])])

# Generate the observed generalibzed histogram


for i, (observationMean, observationalError) in enumerate(zip(observationsSubset[observedDataParameter], observationsSubset[observedErrorParameter])):
#for i, (num, observationMean) in enumerate(zip(num, bins)):
    u = (x - observationMean)/bandwidth
    
    observedHistogram += 1/sqrt(2*np.pi * bandwidth**2) * exp(-0.5 * u**2)
   

    
# Normalise and smooth
observedHistogram /= len(observationsSubset)

# Generate the model generalized histograms
modelHeaders =  {
                    'SgrTriax_DYN.dat' : 'lambda beta ra dec l b xgc ygc zgc xsun ysun zsun x4 y4 z4 u v w dist vgsr mul mub mua mud Pcol Lmflag',
                    'oblate.dat'    :'lambda beta l b xgc ygc zgc xsun ysun zsun u v w dist vgsr Pcol',
                    'prolate.dat'   :'lambda beta l b xgc ygc zgc xsun ysun zsun u v w dist vgsr Pcol',
                    'spherical.dat' :'lambda beta l b xgc ygc zgc xsun ysun zsun u v w dist vgsr Pcol'
                }

models = {}
modelHistograms = {}

for modelFile in modelFiles:

    modelPoints = np.loadtxt(modelFile)
    headers = modelHeaders[modelFile].split()
    spatialIndices = [headers.index(item) for item in modelSpatialLimits]
    
    modelData = []
    modelHistogram = zeros(len(x))

    for modelPoint in modelPoints:
        
        # Check to see if this is in the same spatial bounds as our observed data
        
        spatialIndex = spatialIndices[0]
        observedSpatialLimit = observedSpatialLimits[0]
        
        spatialLimitMin, spatialLimitMax = spatialLimits[0]
       
        
        if ((spatialLimitMax + spatialTolerances[0]) >= modelPoint[spatialIndex]) and (modelPoint[spatialIndex] >= (spatialLimitMin - spatialTolerances[0])):
            
            spatialIndex = spatialIndices[1]
            observedSpatialLimit = observedSpatialLimits[1]
            
            spatialLimitMin, spatialLimitMax = spatialLimits[1]
            
            if ((spatialTolerances[1] + spatialLimitMax) >= modelPoint[spatialIndex]) and (modelPoint[spatialIndex] >= (spatialLimitMin - spatialTolerances[1])) and (3 >= modelPoint[headers.index(modelMarkerParameter)]):
                modelData.append([modelPoint[headers.index(modelParameter)], modelPoint[headers.index(modelDistanceParameter)], modelPoint[headers.index(modelMarkerParameter)]])
                u = (x - modelPoint[headers.index(modelParameter)]) / bandwidth
                modelHistogram += 1/sqrt(2*np.pi * bandwidth**2) * exp(-0.5 * u**2)
                
                
   
    modelHistogram /= len(modelData)
    #modelHistogram /= len(observationsSubset)
                
    modelHistograms[modelFile] = modelHistogram
    
    models[modelFile] = modelData
    
    
            


# Begin plotting 

# Tri-axial -- generalised histogram
markerMap = ['o', '+', 'o', 's', 'd', '^', 'v', '>', '<']
plt.close('all')

#fig, ((triAxialGH, triAxialVD), (prolateGH, prolateVD), (sphericalGH, sphericalVD), (oblateGH, oblateVD)) = plt.subplots(4, 2, sharex=False, sharey=False)
fig, ((triAxialGH, triAxialVD), (sphericalGH, sphericalVD), (oblateGH, oblateVD)) = plt.subplots(3, 2, sharex=False, sharey=False)
fig.subplots_adjust(left=0.08, right=0.92, bottom=0.08, top=0.95, wspace=0.03, hspace=0.05)


triAxialGH.plot(x, observedHistogram, 'k:')
triAxialGH.plot(x, modelHistograms['SgrTriax_DYN.dat'], 'k-')

triAxialGH.set_xticklabels([], alpha=0.0)
triAxialGH.set_xbound(xBound[0], xBound[1])
triAxialGH.set_ybound(0, 16.*10**-3)
triAxialGH.set_yticks([0, 0.004, 0.008, 0.012, .016])
#triAxialGH.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
triAxialGH.set_yticklabels([0, 0.4, 0.8, 1.2, 1.6])
triAxialGH.text(-400, 1.78*10**-2, r"$\times{}10^{-2}$", horizontalalignment='left', verticalalignment='top', fontsize=labelsize)

triAxialGH.tick_params(labelsize=labelsize)


# Tri-axial -- kinematics and distance


triAxialVD.get_xaxis().set_ticks([])


triAxialVD.set_ybound([distBound[0], distBound[1]])

ay = triAxialVD.twinx()
ay.set_yticks(triAxialVD.get_yticks())
ay.set_yticks([0, 10, 20, 30, 40, 50])
ay.set_yticklabels(map(str, [0, 10, 20, 30, 40, 50]))
ay.yaxis.set_tick_params(labelsize=labelsize)

triAxialVD.set_yticklabels([], alpha=0.0)
triAxialVD.set_xticklabels([], alpha=0.0)
ay.tick_params(labelsize=labelsize)

for point in models['SgrTriax_DYN.dat']:
    
    vgsr, distance, lmflag = point
    lmflag = int(lmflag)

    triAxialVD.scatter([vgsr], [distance], marker=markerMap[lmflag], facecolor='w')

triAxialVD.set_xbound(xBound)
triAxialVD.set_ybound(distBound)
"""
# Prolate -- generalised histogram

prolateGH.plot(x, observedHistogram, 'k:')
prolateGH.plot(x, modelHistograms['prolate.dat'], 'k-')


prolateGH.set_xticklabels([], alpha=0.0)
prolateGH.set_xbound(xBound)
prolateGH.set_ybound(0, 16.*10**-3)
prolateGH.set_yticks([0, 0.004, 0.008, 0.012])

prolateGH.set_yticklabels([0, 0.4, 0.8, 1.2])

prolateGH.tick_params(labelsize=labelsize)


# Prolate -- kinematics and distance

prolateVD.get_xaxis().set_ticks([])

prolateVD.set_xbound(xBound)
prolateVD.set_ybound(distBound)

ay = prolateVD.twinx()
ay.set_yticks(prolateVD.get_yticks())
ay.set_yticks([0, 10, 20, 30, 40, 50])

prolateVD.set_yticklabels([], alpha=0.0)
prolateVD.set_xticklabels([], alpha=0.0)
ay.tick_params(labelsize=labelsize)

for point in models['prolate.dat']:

    vgsr, distance, lmflag = point
    lmflag = int(lmflag)
    
    prolateVD.scatter([vgsr], [distance], marker=markerMap[lmflag], facecolor='w')



prolateVD.set_xbound(xBound)
prolateVD.set_ybound(distBound)
"""
# Spherical -- generalised histogram

sphericalGH.plot(x, observedHistogram, 'k:')
sphericalGH.plot(x, modelHistograms['spherical.dat'], 'k-')


sphericalGH.get_xaxis().set_ticks([])
sphericalGH.set_xticklabels([], alpha=0.0)
sphericalGH.set_xbound(xBound)
sphericalGH.set_ybound(0, .8*10**-2)
sphericalGH.set_yticks([0, 0.002, 0.004, 0.006])
sphericalGH.set_yticklabels([0.0, 0.2, 0.4, 0.6])
sphericalGH.tick_params(labelsize=labelsize)


# Spherical -- kinematics and distance

sphericalVD.get_xaxis().set_ticks([])

sphericalVD.set_xbound(xBound)
sphericalVD.set_ybound(distBound)
ay = sphericalVD.twinx()
ay.set_yticks(sphericalVD.get_yticks())
ay.set_yticks([0, 10, 20, 30, 40, 50])
ay.set_yticklabels(map(str, [0, 10, 20, 30, 40, 50]))
ay.yaxis.set_tick_params(labelsize=labelsize)

sphericalVD.set_yticklabels([], alpha=0.0)
sphericalVD.set_xticklabels([], alpha=0.0)

for point in models['spherical.dat']:

    vgsr, distance, lmflag = point
    lmflag = int(lmflag)
    

    sphericalVD.scatter([vgsr], [distance], marker=markerMap[lmflag], facecolor='w')


sphericalVD.set_xbound(xBound)
sphericalVD.set_ybound(distBound)

# Oblate -- generalised histogram

oblateGH.plot(x, observedHistogram, 'k:')
oblateGH.plot(x, modelHistograms['oblate.dat'], 'k-')

oblateGH.set_xbound(xBound)
oblateGH.set_ybound(0, 0.8*10**-2)
oblateGH.set_yticks([0, 0.002, 0.004, 0.006])
oblateGH.set_yticklabels(['0.0', '0.2', '0.4', '0.6'])
oblateGH.set_xticks([-300, -200, -100, 0, 100, 200, 300])
oblateGH.set_xticklabels(map(str, [-300, -200, -100, 0, 100, 200, 300]))
oblateGH.xaxis.set_tick_params(labelsize=labelsize)
oblateGH.tick_params(labelsize=labelsize)

# Oblate -- kinematics and distance

oblateVD.set_xbound(xBound)
oblateVD.set_ybound(distBound)

ay = oblateVD.twinx()
ay.set_yticks(oblateVD.get_yticks())
ay.set_yticks([0., 10., 20., 30., 40., 50.])
ay.set_yticklabels(map(str, [0, 10, 20, 30, 40, 50]))
ay.yaxis.set_tick_params(labelsize=labelsize)

oblateVD.set_yticklabels([], alpha=0.0)

for point in models['oblate.dat']:
    
    vgsr, distance, lmflag = point
    lmflag = int(lmflag)
    
    oblateVD.scatter([vgsr], [distance], marker=markerMap[lmflag], facecolor='w')


oblateVD.set_xbound(xBound)
oblateVD.set_xticks([-300, -200, -100, 0, 100, 200, 300])
oblateVD.set_xticklabels(map(str, [-300, -200, -100, 0, 100, 200, 300]))
oblateVD.xaxis.set_tick_params(labelsize=labelsize)
oblateVD.set_ybound(distBound)
oblateVD.tick_params(labelsize=labelsize)


# Adjust the plot space    
#plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, hspace=0.30, wspace=0.12)
#plt.subplots_adjust(left=0.10, bottom=0.10, right=0.90, top=0.90, hspace=0.03, wspace=0.03)


# Add axes
from matplotlib import rc
rc('text', usetex=True)
oblateVD.set_xlabel(r'$\textrm{V}_{\textrm{\textsc{GSR}}}\;\;(\textrm{km s}^{-1})$', fontsize=fontsize)
oblateGH.set_xlabel(r'$\textrm{V}_{\textrm{\textsc{GSR}}}\;\;(\textrm{km s}^{-1})$', fontsize=fontsize)
sphericalGH.yaxis.set_label_text(r'$\sum{P\left(\textit{V}_{\textrm{\textsc{GSR}}}\right)}$', fontsize=fontsize)
#sphericalGH.yaxis.set_label_coords(-0.13, 1.15)

sphericalVD.yaxis.set_label_position('right')
sphericalVD.yaxis.set_label_coords(1.1, 0.5)
#sphericalVD.yaxis.set_label_coords(1.13, 1.15)
sphericalVD.yaxis.set_label_text(r'$D\;\;(\textrm{kpc})$', fontsize=fontsize)
sphericalVD.yaxis.set_tick_params(labelsize=labelsize)

plt.rcParams['font.size'] = 11.
#triAxialGH.text(xBound[1]*0.92, 0.9*2*10**-3, r"Tri-axial (LM 10)\n$\\q_{(2,1,z)} = (1.00, 1.36, 1.38)$", horizontalalignment='right', verticalalignment='top')
#oblateGH.text(xBound[1]*0.92, 0.9*2*10**-3, r"Oblate (LJM 05)\n $q = 0.90$", horizontalalignment='right', verticalalignment='top')
#prolateGH.text(xBound[1]*0.92, 0.9*2*10**-3, r"Prolate (LJM 05)\n $q = 1.25$", horizontalalignment='right', verticalalignment='top')
#sphericalGH.text(xBound[1]*0.92, 0.9*2*10**-3, r"Spherical (LJM 05)\n $q = 1.00$", horizontalalignment='right', verticalalignment='top')
triAxialVD.text(xBound[1]*0.92, 52, r"Tri-axial (LM10)", horizontalalignment='right', verticalalignment='top')
triAxialVD.text(xBound[1]*0.92, 44, r"%s particles" % (len(models['SgrTriax_DYN.dat'])), horizontalalignment='right', verticalalignment='top')
oblateVD.text(xBound[1]*0.92, 52, r"Oblate (LJM05)", horizontalalignment='right', verticalalignment='top')
oblateVD.text(xBound[1]*0.92, 44, r"%s particles" % len(models['oblate.dat']), horizontalalignment='right', verticalalignment='top')
#prolateVD.text(xBound[1]*0.92, 0.9*60, r"Prolate (LJM 05)", horizontalalignment='right', verticalalignment='top')
sphericalVD.text(xBound[1]*0.92, 52, r"Spherical (LJM05)", horizontalalignment='right', verticalalignment='top')
sphericalVD.text(xBound[1]*0.92, 44, r"%s particles" % len(models['spherical.dat']), horizontalalignment='right', verticalalignment='top')


if useColor:
    plt.savefig('lawvelcompare-color.png')
    plt.savefig('lawvelcompare-color.eps')
    
else:
    plt.savefig('lawvelcompare-bw.eps')
    plt.savefig('lawvelcompare-bw.png')

figGH = plt.figure(2)
figGH.subplots_adjust(left=0.08, bottom=0.11, right=0.95, top=0.95, wspace=0.2, hspace=0.2)
ax = figGH.add_subplot(111)
labelsize = 14.5
fontsize = 13

ax.plot(x, observedHistogram, 'k-', label='Observed distribution')
ax.set_xbound(xBound)

haloMean = 0.0
haloStd = 101.6
y = exp(-((x-haloMean)**2)/(2*haloStd**2))*(1/sqrt(2*pi*haloStd**2))

# Exclude Q sigma data points to reduce observationsSubset length
rect = Rectangle((-135, 0), 100, ax.get_ylim()[1], facecolor='#dddddd', edgecolor='none')
ax.add_patch(rect)

rect = Rectangle((120, 0), 30, ax.get_ylim()[1], facecolor='#dddddd', edgecolor='none')
ax.add_patch(rect)

rect = Rectangle((190, 0), 60, ax.get_ylim()[1], facecolor='#dddddd', edgecolor='none')
ax.add_patch(rect)

#figGH.text(-130., 0.05, 'Feature A')
#ax.text(140., 0.05, 'Feature B')
figGH.text(0.375, 0.87, 'Feature A', fontsize=fontsize)
figGH.text(0.61, 0.43, 'Feature B', fontsize=fontsize)
figGH.text(0.71, 0.34, 'Feature C', fontsize=fontsize)

significance = (observedHistogram - y) * len(observationsSubset)/sqrt(y)

excessSigma = 2.5
significant = filter(lambda alpha: excessSigma > alpha > -excessSigma, significance)
fractionHalo = len(significant) / float(len(observedHistogram))

print "Using %s sigma as excess, there is %s percent halo" % (excessSigma, fractionHalo*100.,)
# Generate the N sigma lines

sigmas = []
#sigmas.append(1) # 67% confidence level
#sigmas.append(2) # 95% confidence level
sigmas.append(3) # 99.7% confidence level

#gpSigmas = [(g*sampleNumber + k*sqrt(g*sampleNumber**2))/sampleNumber**2 for k in sigmas]
gpSigmas = [y + k*sqrt(y)/len(observationsSubset) for k in sigmas]
gmSigmas = [y - k*sqrt(y)/len(observationsSubset) for k in sigmas]

# Some minus (unrealistic) values may result, let's minimize the floor to zero
gmSigmas = [[max(0, gmSigmas[i][j]) for j in range(0, len(gmSigmas[i]))] for i in range(0, len(gmSigmas))]

# Normalise the gaussian distributions
y *= fractionHalo
gmSigmas = [array(k)*fractionHalo for k in gmSigmas]
gpSigmas = [array(k)*fractionHalo for k in gpSigmas]



ax.plot(x,y, 'k--', label='Canonical halo')
for gp, gm in zip(gpSigmas, gmSigmas):
    ax.plot(x, gp, 'k-.', label=r'$\pm{}3\sigma$ level')
    ax.plot(x, gm, 'k-.')
    
ax.set_xticks([-300, -200, -100, 0, 100, 200, 300])
ax.set_xticklabels(map(str, [-300, -200, -100, 0, 100, 200, 300]))
ax.set_yticks([0.0, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006])
ax.set_yticklabels(map(str, [0, 1, 2, 3, 4, 5, 6]))
ax.text(-390, 0.063, r"$\times{}10^{-3}$", horizontalalignment='left', verticalalignment='top', fontsize=labelsize)
ax.set_xlabel(r'$\textrm{V}_{\textrm{\textsc{GSR}}}\;\;(\textrm{km s}^{-1})$', fontsize=labelsize)
ax.set_ylabel(r'$\sum{P\left(\textit{V}\,\right)}$', fontsize=labelsize)
ax.xaxis.set_tick_params(labelsize=labelsize)
ax.yaxis.set_tick_params(labelsize=labelsize)
ax.legend(loc=1)
plt.draw()
plt.savefig('velocity-halo.eps')
# BESANCON MODEL

raise InputError

bes = plt.figure()
ax = bes.add_subplot(111)

ax.plot(x, observedHistogram, 'k-')
ax.set_xbound(xBound)

besanconData = loadtxt('besancon.dat')
lbound = [274.40, 307.19]
bbound = [59.17, 61.84]

b = (bbound[1] - bbound[0])/2 + bbound[0]
l = (lbound[1] - lbound[0])/2 + lbound[0]

vgsrOffset = 220.0*sin(radians(l))*cos(radians(b)) + 16.5*(sin(radians(b))*sin(radians(25)) + cos(radians(b))*cos(radians(25))*cos(radians(l-53)))
print "vgsr offset", vgsrOffset
vgsrBounds = [220.0*sin(radians(l))*cos(radians(b)) + 16.5*(sin(radians(b))*sin(radians(25)) + cos(radians(b))*cos(radians(25))*cos(radians(l-53))) for l, b in zip(lbound, bbound)]
print vgsrBounds


vr = besanconData[:,15]
vgsr = vr + vgsrOffset

besanconHistogram = zeros(len(observedHistogram))

for v in vgsr:
#for i, (num, observationMean) in enumerate(zip(num, bins)):
    u = (x - v)/bandwidth
    
    besanconHistogram += 1/sqrt(2*np.pi * bandwidth**2) * exp(-0.5 * u**2)
    
besanconHistogram /= len(vgsr)

ax.plot(x, besanconHistogram, 'k--')

bpSigmas = [besanconHistogram + k*sqrt(besanconHistogram)/len(besanconHistogram) for k in sigmas]
bmSigmas = [besanconHistogram - k*sqrt(besanconHistogram)/len(besanconHistogram) for k in sigmas]



# Some minus (unrealistic) values may result, let's minimize the floor to zero
bmSigmas = [[max(0, bmSigmas[i][j]) for j in range(0, len(bmSigmas[i]))] for i in range(0, len(bmSigmas))]

# Normalise the gaussian distributions


#for bp, bm in zip(bpSigmas, bmSigmas):
#    ax.plot(x, bp, 'k-.')
#    ax.plot(x, bm, 'k-.')


ax.set_xticks([-300, -200, -100, 0, 100, 200, 300])
ax.set_xlabel(r'$\textrm{V}_{\textrm{\textsc{GSR}}}\;\;(\textrm{km s}^{-1})$', fontsize=fontsize)
ax.set_ylabel(r'$\sum{P\left(\textit{V}\,\right)}$', fontsize=fontsize)

# Determine the chi-squared value
sd = 6.7
fig = plt.figure()
for model in ['oblate.dat', 'spherical.dat', 'SgrTriax_DYN.dat']:
    E, bins, e = plt.hist(transpose(models[model])[0], bins=ceil(-xBound[0] + xBound[1])/bandwidth, range=xBound)
    O, bins, e = plt.hist(observationsSubset[observedDataParameter], bins=ceil(-xBound[0] + xBound[1])/bandwidth, range=xBound)
    
    #chiSq = sum(((observedHistogram[0:470] - 471*modelHistograms[model][0:470])**2)/(471*modelHistograms[model][0:470]))
    #print O, E
    #chiSq = sum(((O - E)**2)/E)
    
    chiSq = sum((O - E)**2)/float(len(O) - 3)

    print model, chiSq
    #print model, chiSq, R

sigfig = plt.figure()
ax = sigfig.add_subplot(111)

sig = (observedHistogram - y)/gpSigmas[0]
for i, sig_ in enumerate(sig):
    if sig_ > 0:
        sig[i] = sig_ + 3.
    else: sig[i] = sig_ -3.0


plot(x,sig)
