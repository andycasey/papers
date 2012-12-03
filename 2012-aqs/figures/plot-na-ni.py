import matplotlib.pyplot as plt
import numpy as np

# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Na_Fe = []
e_Na_Fe = []

Ni_Fe = []
e_Ni_Fe = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Na')
	_Na_Fe, _e_Na_Fe = data[idx, 1:]

	idx = list(data[:, 0]).index('Ni')
	_Ni_Fe, _e_Ni_Fe = data[idx, 1:]

	Na_Fe.append(_Na_Fe)
	Ni_Fe.append(_Ni_Fe)
	e_Na_Fe.append(_e_Na_Fe)
	e_Ni_Fe.append(_e_Ni_Fe)


# Clean up Na_Fe, Ni_Fe, e_Na_Fe, e_Ni_Fe
Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Na_Fe])
Ni_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Ni_Fe])

e_Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Na_Fe])
e_Ni_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Ni_Fe])

Na_Fe, Ni_Fe, e_Na_Fe, e_Ni_Fe = map(np.array, [Na_Fe, Ni_Fe, e_Na_Fe, e_Ni_Fe])


A = np.vstack([Na_Fe, np.ones(len(Na_Fe))]).T
m, c = np.linalg.lstsq(A, Ni_Fe)[0]

fig = plt.figure()

ax = fig.add_subplot(111)
ax.scatter(Na_Fe, Ni_Fe, color='k')
ax.errorbar(Na_Fe, Ni_Fe, xerr=e_Na_Fe, yerr=e_Ni_Fe, fmt=None, ecolor='k')


x = np.array([np.min(Na_Fe), np.max(Na_Fe)])
ax.plot(x, m * x + c, 'k-')

ax.set_xlabel('[Na/Fe]')
ax.set_ylabel('[Ni/Fe]')

ax.set_xlim(-0.1, 0.5)
ax.set_ylim(-0.7, 0.5)

plt.savefig('aquarius-na-ni.eps')
plt.savefig('aquarius-na-ni.pdf')
