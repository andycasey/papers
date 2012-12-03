import matplotlib.pyplot as plt
import numpy as np

# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Mg_Fe = []
e_Mg_Fe = []

Al_Fe = []
e_Al_Fe = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Mg')
	_Mg_Fe, _e_Mg_Fe = data[idx, 1:]

	idx = list(data[:, 0]).index('Al')
	_Al_Fe, _e_Al_Fe = data[idx, 1:]

	Mg_Fe.append(_Mg_Fe)
	Al_Fe.append(_Al_Fe)
	e_Mg_Fe.append(_e_Mg_Fe)
	e_Al_Fe.append(_e_Al_Fe)


# Clean up Mg_Fe, Al_Fe, e_Mg_Fe, e_Al_Fe
Mg_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Mg_Fe])
Al_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Al_Fe])

e_Mg_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Mg_Fe])
e_Al_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Al_Fe])

Mg_Fe, Al_Fe, e_Mg_Fe, e_Al_Fe = map(np.array, [Mg_Fe, Al_Fe, e_Mg_Fe, e_Al_Fe])


A = np.vstack([Mg_Fe, np.ones(len(Mg_Fe))]).T
m, c = np.linalg.lstsq(A, Al_Fe)[0]

fig = plt.figure()

ax = fig.add_subplot(111)
ax.scatter(Mg_Fe, Al_Fe, color='k')
ax.errorbar(Mg_Fe, Al_Fe, xerr=e_Mg_Fe, yerr=e_Al_Fe, fmt=None, ecolor='k')


x = np.array([np.min(Mg_Fe), np.max(Mg_Fe)])
ax.plot(x, m * x + c, 'k-')

ax.set_xlabel('[Mg/Fe]')
ax.set_ylabel('[Al/Fe]')

ax.set_xlim(-0.5, 1.0)
ax.set_ylim(-0.5, 1.5)

plt.savefig('aquarius-mg-al.eps')
plt.savefig('aquarius-mg-al.pdf')
