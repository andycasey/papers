import matplotlib.pyplot as plt
import numpy as np

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

	Fe_H.append(_Fe_H)
	e_Fe_H.append(_e_Fe_H)



# Clean up Mg_Fe, Fe_H, e_Mg_Fe, e_Fe_H
Sc_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Sc_Fe])
V_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in V_Fe])
Cr_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Cr_Fe])
Mn_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Mn_Fe])
Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Fe_H])

e_Sc_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Sc_Fe])
e_V_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_V_Fe])
e_Cr_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Cr_Fe])
e_Mn_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Mn_Fe])
e_Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Fe_H])

Sc_Fe, Cr_Fe, V_Fe, Mn_Fe, Fe_H = map(np.array, [Sc_Fe, Cr_Fe, V_Fe, Mn_Fe, Fe_H])
e_Sc_Fe, e_Cr_Fe, e_V_Fe, e_Mn_Fe, e_Fe_H = map(np.array, [e_Sc_Fe, e_Cr_Fe, e_V_Fe, e_Mn_Fe, e_Fe_H])


print "[Sc II/Fe] = ", Sc_Fe
print "[V/Fe] = ", V_Fe
print "[Cr/Fe] = ", Cr_Fe
print "[Mn/Fe] = ", Mn_Fe

# Bounds
xlim = np.array([-1.7, -0.5])
ylim_Sc = [-0.1, 0.4]
ylim_V = [-0.1, 0.5]
ylim_Cr = [-0.5, 0.1]
ylim_Mn = [-0.7, 0.6]


fig = plt.figure()
fig.subplots_adjust(left=0.10,bottom=0.10,right=0.95, top=0.95, hspace=0.0)

# Scantium
ax = fig.add_subplot(411)
ax.scatter(Fe_H, Sc_Fe, color='k')
ax.errorbar(Fe_H, Sc_Fe, xerr=e_Fe_H, yerr=e_Sc_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Sc_Fe)[0]

print "Sc II: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Sc)
ax.set_xticklabels([])
ax.set_ylabel('[Sc/Fe]')
ax.set_yticks([0.0, 0.2, 0.4])



# Vandium
ax = fig.add_subplot(412)
ax.scatter(Fe_H, V_Fe, color='k')
ax.errorbar(Fe_H, V_Fe, xerr=e_Fe_H, yerr=e_V_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, V_Fe)[0]

print "V: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_V)
ax.set_xticklabels([])
ax.set_ylabel('[V/Fe]')
ax.set_yticks([0.0, 0.2, 0.4])


# Cromium
ax = fig.add_subplot(413)
ax.scatter(Fe_H, Cr_Fe, color='k')
ax.errorbar(Fe_H, Cr_Fe, xerr=e_Fe_H, yerr=e_Cr_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Cr_Fe)[0]

print "Cr: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Cr)

ax.set_ylabel('[Cr/Fe]')
ax.set_xlabel('[Fe/H]')
ax.set_yticks([-0.4, -0.2, 0.0])


# Manganese
ax = fig.add_subplot(414)
ax.scatter(Fe_H, Mn_Fe, color='k')
ax.errorbar(Fe_H, Mn_Fe, xerr=e_Fe_H, yerr=e_Mn_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Mn_Fe)[0]

print "Mn: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_Mn)

ax.set_ylabel('[Mn/Fe]')
ax.set_xlabel('[Fe/H]')
ax.set_yticks([-0.60, -0.4, -0.20, 0.0, 0.2, 0.4])




plt.savefig('aquarius-fe-peak.eps')
plt.savefig('aquarius-fe-peak.pdf')
