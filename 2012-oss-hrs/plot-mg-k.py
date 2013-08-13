import numpy as np


mucciarelli_headers = 'ID Teff logg vt [Fe/H] [Fe/H]_err [Mg/Fe] [Mg/Fe]_err [K/Fe] [K/Fe]_err [Ti/Fe] [Ti/Fe]_err [Ca/Fe] [Ca/Fe]_err'

mg_fe, mg_fe_err, k_fe, k_fe_err = range(4)
mucciarelli_columns = ['[Mg/Fe]', '[Mg/Fe]_err', '[K/Fe]', '[K/Fe]_err']

mucciarelli_column_indices = [mucciarelli_headers.split().index(item) for item in mucciarelli_columns]
mucciarelli = np.loadtxt('ngc2419-mucciarelli.txt', usecols=mucciarelli_column_indices)

fig = plt.figure()
fig.subplots_adjust(left=0.10, right=0.95, bottom=0.07, wspace=0, hspace=0.14, top=0.95)
ax = fig.add_subplot(111)

ax.errorbar(mucciarelli[:, mg_fe], mucciarelli[:, k_fe], xerr=mucciarelli[:, mg_fe_err], yerr=mucciarelli[:, k_fe_err], fmt=None, ecolor='k')
ax.scatter(mucciarelli[:, mg_fe], mucciarelli[:, k_fe], facecolor='k', edgecolor='k', label='Mucciarelli et al. (DEIMOS; 2012)')

# Cohen data
cohen = np.loadtxt('ngc2419-cohen.txt', usecols=(1,2,3,4, ))

ax.errorbar(cohen[:, mg_fe], cohen[:, k_fe], xerr=cohen[:, mg_fe_err], yerr=cohen[:, k_fe_err], fmt=None, ecolor='g')
ax.scatter(cohen[:, mg_fe], cohen[:, k_fe], marker='s', edgecolor='g', facecolor='g', label='Cohen et al. (HIRES; 2011, 2012)')


ax.set_ylabel('[K/Fe]')
ax.set_xlabel('[Mg/Fe]')

ax.legend(loc=3)

