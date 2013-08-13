import matplotlib.pyplot as plt
import numpy as np

aquarius_colors = ['c', 'g', 'm', 'r', 'b']

# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Sc_Fe = []
e_Sc_Fe = []

V_Fe = []
e_V_Fe = []

Cr_Fe = []
e_Cr_Fe = []

Mn_Fe = []
e_Mn_Fe = []


Co_Fe = []
e_Co_Fe = []

Ni_Fe = []
e_Ni_Fe = []

Cu_Fe = []
e_Cu_Fe = []


Fe_H = []
e_Fe_H = []



observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -6, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Sc') + 1 # (Use Sc II)
	_Sc_Fe, _e_Sc_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('V')
	_V_Fe, _e_V_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Cr')
	_Cr_Fe, _e_Cr_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Mn')
	_Mn_Fe, _e_Mn_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Co')
	_Co_Fe, _e_Co_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Ni')
	_Ni_Fe, _e_Ni_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Cu')
	_Cu_Fe, _e_Cu_Fe = data[idx, 2:]

	idx = list(data[:, 0]).index('Fe')
	_Fe_H, _e_Fe_H = data[idx, 1], data[idx, -1]

	Sc_Fe.append(_Sc_Fe)
	e_Sc_Fe.append(_e_Sc_Fe)

	V_Fe.append(_V_Fe)
	e_V_Fe.append(_e_V_Fe)

	Cr_Fe.append(_Cr_Fe)
	e_Cr_Fe.append(_e_Cr_Fe)

	Mn_Fe.append(_Mn_Fe)
	e_Mn_Fe.append(_e_Mn_Fe)

	Co_Fe.append(_Co_Fe)
	e_Co_Fe.append(_e_Co_Fe)

	Ni_Fe.append(_Ni_Fe)
	e_Ni_Fe.append(_e_Ni_Fe)

	Cu_Fe.append(_Cu_Fe)
	e_Cu_Fe.append(_e_Cu_Fe)

	Fe_H.append(_Fe_H)
	e_Fe_H.append(_e_Fe_H)



# Clean up Mg_Fe, Fe_H, e_Mg_Fe, e_Fe_H
Sc_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Sc_Fe])
V_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in V_Fe])
Cr_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Cr_Fe])
Mn_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Mn_Fe])
Co_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Co_Fe])
Ni_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Ni_Fe])
Cu_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Cu_Fe])
Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Fe_H])

e_Sc_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Sc_Fe])
e_V_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_V_Fe])
e_Cr_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Cr_Fe])
e_Mn_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Mn_Fe])
e_Co_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Co_Fe])
e_Ni_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Ni_Fe])
e_Cu_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Cu_Fe])
e_Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Fe_H])


Sc_Fe, Cr_Fe, V_Fe, Mn_Fe, Co_Fe, Ni_Fe, Cu_Fe, Fe_H = map(np.array, [Sc_Fe, Cr_Fe, V_Fe, Mn_Fe, Co_Fe, Ni_Fe, Cu_Fe, Fe_H])
e_Sc_Fe, e_Cr_Fe, e_V_Fe, e_Mn_Fe, e_Co_Fe, e_Ni_Fe, e_Cu_Fe, e_Fe_H = map(np.array, [e_Sc_Fe, e_Cr_Fe, e_V_Fe, e_Mn_Fe, e_Co_Fe, e_Ni_Fe, e_Cu_Fe, e_Fe_H])


print "[Sc II/Fe] = ", Sc_Fe
print "[V/Fe] = ", V_Fe
print "[Cr/Fe] = ", Cr_Fe
print "[Mn/Fe] = ", Mn_Fe
print "[Co/Fe] = ", Co_Fe
print "[Ni/Fe] = ", Ni_Fe
print "[Cu/Fe] = ", Cu_Fe


# Bounds
xlim = np.array([-1.7, -0.5])
ylim_Sc = [-0.3, 0.3]
ylim_V = [-0.2, 0.5]
ylim_Cr = [-0.5, 0.1]
ylim_Mn = [-0.7, 0.05]

ylim_Co = [-0.2, 0.15]
ylim_Ni = [-0.3, 0.3]
ylim_Cu = [-1.0, 0.4]

num_plots = 6

fig = plt.figure(figsize=(8, 10))
fig.subplots_adjust(left=0.10,bottom=0.10,right=0.99, top=0.99, hspace=0.0)

# Scantium
ax = fig.add_subplot(num_plots, 1, 1)

ax.errorbar(Fe_H, Sc_Fe, xerr=e_Fe_H, yerr=e_Sc_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, Sc_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Sc_Fe)[0]

print "Sc II: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Sc)
ax.set_xticklabels([])
ax.set_ylabel('[Sc/Fe]')
ax.set_yticks([-0.2, 0.0, 0.2, 0.4])



# Vandium
ax = fig.add_subplot(num_plots, 1, 2)

ax.errorbar(Fe_H, V_Fe, xerr=e_Fe_H, yerr=e_V_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, V_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, V_Fe)[0]

print "V: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_V)
ax.set_xticklabels([])
ax.set_ylabel('[V/Fe]')
ax.set_yticks([-0.2, 0.0, 0.2, 0.4])


# Cromium
ax = fig.add_subplot(num_plots, 1, 3)

ax.errorbar(Fe_H, Cr_Fe, xerr=e_Fe_H, yerr=e_Cr_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, Cr_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Cr_Fe)[0]

print "Cr: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Cr)

ax.set_ylabel('[Cr/Fe]')
ax.set_xlabel('[Fe/H]')
ax.set_yticks([-0.4, -0.2, 0.0])


# Manganese
ax = fig.add_subplot(num_plots, 1, 4)

ax.errorbar(Fe_H, Mn_Fe, xerr=e_Fe_H, yerr=e_Mn_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, Mn_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Mn_Fe)[0]

print "Mn: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Mn)

ax.set_ylabel('[Mn/Fe]')
ax.set_xlabel('[Fe/H]')
ax.set_yticks([-0.6, -0.3, 0.0, 0.3])


# Cobalt
ax = fig.add_subplot(num_plots, 1, 5)

ax.errorbar(Fe_H, Co_Fe, xerr=e_Fe_H, yerr=e_Co_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, Co_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)


A = np.vstack([Fe_H[np.isfinite(Co_Fe)], np.ones(np.sum([np.isfinite(Co_Fe)]))]).T
m, c = np.linalg.lstsq(A, Co_Fe[np.isfinite(Co_Fe)])[0]

print "Co: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Co)

ax.set_ylabel('[Co/Fe]')
ax.set_xlabel('[Fe/H]')
ax.set_yticks([-0.20, -0.1, 0.0, 0.1])

'''
# Nickel
ax = fig.add_subplot(num_plots, 1, 6)

ax.errorbar(Fe_H, Ni_Fe, xerr=e_Fe_H, yerr=e_Ni_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, Ni_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Ni_Fe)[0]

print "Ni: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Ni)

ax.set_ylabel('[Ni/Fe]')
ax.set_xlabel('[Fe/H]')
'''

#Copper
ax = fig.add_subplot(num_plots, 1, 6)

ax.errorbar(Fe_H, Cu_Fe, xerr=e_Fe_H, yerr=e_Cu_Fe, fmt=None, ecolor='k', zorder=-32)
ax.scatter(Fe_H, Cu_Fe, color=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Cu_Fe)[0]

print "Cu: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, '-', c='#666666', lw=1.2)
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Cu)

ax.set_ylabel('[Cu/Fe]')
ax.set_xlabel('[Fe/H]')
ax.set_yticks([-1.0, -0.5, 0.0])



plt.savefig('aquarius-fe-peak.eps')
plt.savefig('aquarius-fe-peak.pdf')
