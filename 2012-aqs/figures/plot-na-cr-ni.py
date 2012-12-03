import matplotlib.pyplot as plt
import numpy as np

# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Na_Fe = []
e_Na_Fe = []

Cr_Fe = []
e_Cr_Fe = []

Ni_Fe = []
e_Ni_Fe = []

Fe_H = []
e_Fe_H = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -6, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Na')
	_Na_Fe, _e_Na_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Cr')
	_Cr_Fe, _e_Cr_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Ni')
	_Ni_Fe, _e_Ni_Fe = data[idx, 2:]

	idx = list(data[:, 0]).index('Fe')
	_Fe_H, _e_Fe_H = data[idx, 1], data[idx, -1]

	Na_Fe.append(_Na_Fe)
	e_Na_Fe.append(_e_Na_Fe)

	Cr_Fe.append(_Cr_Fe)
	e_Cr_Fe.append(_e_Cr_Fe)

	Ni_Fe.append(_Ni_Fe)
	e_Ni_Fe.append(_e_Ni_Fe)

	Fe_H.append(_Fe_H)
	e_Fe_H.append(_e_Fe_H)



# Clean up Mg_Fe, Fe_H, e_Mg_Fe, e_Fe_H
Na_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Na_Fe])
Cr_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Cr_Fe])
Ni_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Ni_Fe])
Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Fe_H])

e_Na_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Na_Fe])
e_Cr_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Cr_Fe])
e_Ni_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Ni_Fe])
e_Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Fe_H])

Na_Fe, Ni_Fe, Cr_Fe, Fe_H = map(np.array, [Na_Fe, Ni_Fe, Cr_Fe, Fe_H])
e_Na_Fe, e_Ni_Fe, e_Cr_Fe, e_Fe_H = map(np.array, [e_Na_Fe, e_Ni_Fe, e_Cr_Fe, e_Fe_H])


print "[Na/Fe] = ", Na_Fe
print "[Cr/Fe] = ", Cr_Fe
print "[Ni/Fe] = ", Ni_Fe

# Bounds
xlim = np.array([-1.8, -0.4])
ylim_na = [-0.3, 0.4]
ylim_cr = [-0.4, 0.1]
ylim_ni = [-0.2, 0.2]


fig = plt.figure()

# Sodium

ax = fig.add_subplot(311)
ax.scatter(Fe_H, Na_Fe, color='k')
ax.errorbar(Fe_H, Na_Fe, xerr=e_Fe_H, yerr=e_Na_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Na_Fe)[0]

print "Na: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_na)
ax.set_xticklabels([])
ax.set_ylabel('[Na/Fe]')

# Cromium

ax = fig.add_subplot(312)
ax.scatter(Fe_H, Cr_Fe, color='k')
ax.errorbar(Fe_H, Cr_Fe, xerr=e_Fe_H, yerr=e_Cr_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Cr_Fe)[0]

print "Cr: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_cr)
ax.set_xticklabels([])
ax.set_ylabel('[Cr/Fe]')

# Nickel

ax = fig.add_subplot(313)
ax.scatter(Fe_H, Ni_Fe, color='k')
ax.errorbar(Fe_H, Ni_Fe, xerr=e_Fe_H, yerr=e_Ni_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Ni_Fe)[0]

print "Ni: y = %1.2fx + %1.2f" % (m, c, )

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_ylim(ylim_ni)

ax.set_ylabel('[Ni/Fe]')
ax.set_xlabel('[Fe/H]')





plt.savefig('aquarius-na-cr-ni-fe.eps')
plt.savefig('aquarius-na-cr-ni-fe.pdf')
