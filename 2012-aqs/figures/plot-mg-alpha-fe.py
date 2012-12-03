import matplotlib.pyplot as plt
import numpy as np

# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Si_Fe = []
e_Si_Fe = []

Ca_Fe = []
e_Ca_Fe = []

Ti_Fe = []
e_Ti_Fe = []

Mg_Fe = []
e_Mg_Fe = []

Fe_H = []
e_Fe_H = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -6, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Si')
	_Si_Fe, _e_Si_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Ca')
	_Ca_Fe, _e_Ca_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Mg')
	_Mg_Fe, _e_Mg_Fe = data[idx, 2:]

	idx = list(data[:,0]).index('Ti') + 1 # Get Ti II because it's more reliable
	_Ti_Fe, _e_Ti_Fe = data[idx, 2:]

	idx = list(data[:, 0]).index('Fe')
	_Fe_H, _e_Fe_H = data[idx, 1], data[idx, -1]

	Si_Fe.append(_Si_Fe)
	e_Si_Fe.append(_e_Si_Fe)

	Ca_Fe.append(_Ca_Fe)
	e_Ca_Fe.append(_e_Ca_Fe)

	Mg_Fe.append(_Mg_Fe)
	e_Mg_Fe.append(_e_Mg_Fe)

	Ti_Fe.append(_Ti_Fe)
	e_Ti_Fe.append(_e_Ti_Fe)

	Fe_H.append(_Fe_H)
	e_Fe_H.append(_e_Fe_H)



# Clean up Mg_Fe, Fe_H, e_Mg_Fe, e_Fe_H
Si_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Si_Fe])
Ca_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Ca_Fe])
Ti_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Ti_Fe])
Mg_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Mg_Fe])
Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in Fe_H])

e_Si_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Si_Fe])
e_Ca_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Ca_Fe])
e_Ti_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Ti_Fe])
e_Mg_Fe = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Mg_Fe])
e_Fe_H = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in e_Fe_H])

Si_Fe, Ti_Fe, Ca_Fe, Mg_Fe, Fe_H = map(np.array, [Si_Fe, Ti_Fe, Ca_Fe, Mg_Fe, Fe_H])
e_Si_Fe, e_Ti_Fe, e_Ca_Fe, e_Mg_Fe, e_Fe_H = map(np.array, [e_Si_Fe, e_Ti_Fe, e_Ca_Fe, e_Mg_Fe, e_Fe_H])

# Bounds
xlim = np.array([-1.8, -0.4])
ylim = np.array([-0.2, 0.85])

fig = plt.figure()

# Magnesium

ax = fig.add_subplot(211)
ax.scatter(Fe_H, Mg_Fe, color='k')
ax.errorbar(Fe_H, Mg_Fe, xerr=e_Fe_H, yerr=e_Mg_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, Mg_Fe)[0]

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_xlim(xlim)
ax.set_xticklabels([])
ax.set_ylabel('[Mg/Fe]')


# Alpha median

ax = fig.add_subplot(212, sharey=ax)
alpha_Fe = (Si_Fe + Ti_Fe + Ca_Fe)/3.
e_alpha_Fe = pow(pow(e_Si_Fe, 2) + pow(e_Ti_Fe, 2) + pow(e_Ca_Fe, 2), 0.5)

ax.scatter(Fe_H, alpha_Fe, color='k')
ax.errorbar(Fe_H, alpha_Fe, xerr=e_Fe_H, yerr=e_alpha_Fe, fmt=None, ecolor='k')

A = np.vstack([Fe_H, np.ones(len(Fe_H))]).T
m, c = np.linalg.lstsq(A, alpha_Fe)[0]

x = np.array([np.min(Fe_H), np.max(Fe_H)])
ax.plot(x, m * x + c, 'k-')
ax.plot(xlim, [0, 0], 'k:') # Zero line

ax.set_ylabel('[$\\alpha$/Fe]')
ax.set_xlabel('[Fe/H]')

ax.set_xlim(xlim)
ax.set_ylim(ylim)



plt.savefig('aquarius-mg-alpha-fe.eps')
plt.savefig('aquarius-mg-alpha-fe.pdf')
