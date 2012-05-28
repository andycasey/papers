
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from numpy import *
import atpy


filename = 'combined.xml'
column = '[Fe/H] (Battaglia)'

table = atpy.Table(filename, verbose=False)

histData1 = []
histData2 = []
histData3 = []
histData4 = []
histData5 = []

for item in table.data:
    #if(item[column] > -3.0) and (item[column] < -0.5) and item['K-Giants'] and (-350. < item['V_{gsr} (fxcor)']) and (item['V_{gsr} (fxcor)'] < 0.):
    if (item[column] > -3.0) and (item[column] < -0.5) and item['K-Giants'] and (-140. < item['V_{gsr} (fxcor)']) and (item['V_{gsr} (fxcor)'] < -30.):
    #if (item[column] > -3.0) and (item[column] < -0.5) and item['K-Giants'] and (item['V_{gsr} (fxcor)'] > -300.):# (item['V_{gsr} (fxcor)'] < 0.0):
    #iif (item['K-Giants']):
        if (-120. < item['V_{gsr} (fxcor)']) and (item['V_{gsr} (fxcor)'] < -30.):
            histData1.append(item[column])
        else:
            histData5.append(item[column])
            
        if (-54. < item['V_{gsr} (fxcor)']) and (item['V_{gsr} (fxcor)'] < -34.):
            histData2.append(item[column])
            
        if (-91. < item['V_{gsr} (fxcor)']) and (item['V_{gsr} (fxcor)'] < -71.):
            histData3.append(item[column])
            
        
        histData4.append(item[column])
            
            
rwidth = 0.90
edgecolor = '#000000'
facecolor = 'none'
bins = arange(-2.6, -0.8, 0.20)


plt.close('all')

fig = plt.figure()

ax1 = fig.add_subplot(211)

a1, b, c = ax1.hist(histData2, bins=bins, rwidth=rwidth, edgecolor=edgecolor, facecolor=facecolor)



ax2 = fig.add_subplot(212)
a2, b, c = ax2.hist(histData3, bins=bins, rwidth=rwidth, edgecolor=edgecolor, facecolor=facecolor)

ax2.set_xlabel('[Fe/H]')
ax1.set_ylabel('Number')
ax2.set_ylabel('Number')
plt.setp(ax1.get_xticklabels(), visible=False)

[axis.set_ylim(0, 1 + max(max(_) for _ in [a1, a2])) for axis in [ax1, ax2]]

ax1.text(-1.4, 7, r'$V_{GSR} = -44$ km s$^{-1}$ peak')
ax1.text(-1.4, 6, r'$-54 < V_{GSR} < -34$ km s$^{-1}$')

ax2.text(-1.4, 7, r'$V_{GSR} = -81$ km s$^{-1}$ peak')
ax2.text(-1.4, 6, r'$-91 < V_{GSR} < -71$ km s$^{-1}$')

ax1.get_xticklabels()[0].set_visible(False)
ax1.get_xticklabels()[-1].set_visible(False)
ax2.get_xticklabels()[0].set_visible(False)
ax2.get_xticklabels()[-1].set_visible(False)


plt.savefig('feh-peak-bins.eps')


fig = plt.figure()
"""
ax5 = fig.add_subplot(311)

a5, b, c = ax5.hist(histData1, bins=bins, rwidth=rwidth, edgecolor=edgecolor, facecolor=facecolor)
ax5.set_ylabel('Number')


plt.setp(ax5.get_xticklabels(), visible=False)

ax6 = fig.add_subplot(312)

a6, b, c = ax6.hist(histData5, bins=bins, rwidth=rwidth, edgecolor=edgecolor, facecolor=facecolor)


plt.setp(ax6.get_xticklabels(), visible=False)
ax6.set_ylabel('Number')

ax5 = fig.add_subplot(211)

a5, b, c = ax5.hist(histData1, bins=bins, rwidth=rwidth, edgecolor=edgecolor, facecolor=facecolor)
ax5.set_ylabel('Number')


plt.setp(ax5.get_xticklabels(), visible=False)
"""
labelsize=14
textsize=15
facecolor='#cccccc'
ax7 = fig.add_subplot(111)
fig.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.95)
a7, b, c = ax7.hist(histData4, bins=np.arange(-2.4, -0.4, 0.2),rwidth=rwidth, edgecolor=edgecolor, facecolor=facecolor)

ax7.set_xbound(-2.5, -0.5)
ax7.set_xticks(np.arange(-2.4, 0., 0.2))
ax7.set_xticklabels(map(str, ['', -2.2, '', -1.8, '', -1.4, '', -1.0, '', -0.6, '']))
ax7.xaxis.set_tick_params(labelsize=labelsize)

ax7.set_ybound(0, 25)
ax7.set_yticks([0, 5, 10, 15, 20, 25])
ax7.set_yticklabels(map(str, ['', 5, 10, 15, 20, 25]))
ax7.yaxis.set_tick_params(labelsize=labelsize)

ax7.set_ylabel('Number',fontsize=textsize)

ax7.text(-0.5, 22, "Sgr K-Giant Sample", horizontalalignment='right', fontsize=labelsize)
ax7.text(-0.5, 20, "$\\textrm{V}_{\\textrm{GSR}} < 0$ km s$^{-1}$", horizontalalignment='right', fontsize=13)

[axis.set_ylim(0, 3 + max(max(_) for _ in [a7])) for axis in [ax7]]


ax7.set_xlabel('[Fe/H]', fontsize=textsize)
plt.savefig('fehhist.eps')

