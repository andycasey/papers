import matplotlib.pyplot as plt

#NGC;ID;Bmag;Vmag;Icmag;Kmag;RV11;RV13;Teff;[FeI/H];e_[FeI/H];[O/Fe];e_[O/Fe];l_[Na/Fe];[Na/Fe];e_[Na/Fe];RAJ2000;DEJ2000
# ; ;mag;mag;mag;mag;km/s;km/s;K;[Sun];[Sun];[Sun];[Sun]; ;[Sun];[Sun];"h:m:s";"d:m:s"


caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_l_Na_Fe, caretta_Na_Fe, caretta_e_Na_Fe = np.loadtxt('caretta-et-al-2009-O-Na.txt', skiprows=67, delimiter=';', \
														usecols=(9, 10, 11, 12, 13, 14, 15, ), unpack=True, dtype=str)

# Replace invalid items with np.nans

def clean(array):

	new_array = []
	for item in array:
		try:
			item = float(item)

		except ValueError:
			new_array.append(np.nan)

		else:
			new_array.append(item)

	return np.array(new_array)


caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_Na_Fe, caretta_e_Na_Fe = map(clean, [caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_Na_Fe, caretta_e_Na_Fe])



# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Na_Fe = []
e_Na_Fe = []

O_Fe = []
e_O_Fe = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Na')
	_Na_Fe, _e_Na_Fe = data[idx, 1:]

	idx = list(data[:, 0]).index('O')
	_O_Fe, _e_O_Fe = data[idx, 1:]

	Na_Fe.append(_Na_Fe)
	O_Fe.append(_O_Fe)
	e_Na_Fe.append(_e_Na_Fe)
	e_O_Fe.append(_e_O_Fe)


# Clean up Na_Fe, O_Fe, e_Na_Fe, e_O_Fe
Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Na_Fe])
O_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in O_Fe])

e_Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Na_Fe])
e_O_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_O_Fe])

Na_Fe, O_Fe, e_Na_Fe, e_O_Fe = map(np.array, [Na_Fe, O_Fe, e_Na_Fe, e_O_Fe])


fig = plt.figure()
ax = fig.add_subplot(111)

# Oxygen upper limits
upper_limits = np.where(caretta_e_O_Fe == 9.999)[0]
measurements = np.array(list(set(range(len(caretta_O_Fe))).difference(upper_limits)))

ax.errorbar(caretta_O_Fe[upper_limits], caretta_Na_Fe[upper_limits], xerr=0.02, xlolims=True, fmt=None, ecolor='#666666', zorder=-1)
#ax.errorbar(caretta_O_Fe[upper_limits], caretta_Na_Fe[upper_limits], facecolor='none', edgecolor='k')
#ax.errorbar(caretta_O_Fe[measurements], caretta_Na_Fe[measurements], xerr=caretta_e_O_Fe[measurements], yerr=caretta_e_Na_Fe[measurements], fmt=None)

ax.scatter(caretta_O_Fe[measurements], caretta_Na_Fe[measurements], facecolor='none', edgecolor='#666666', zorder=-1)

ax.set_xlabel('[O/Fe]')
ax.set_ylabel('[Na/Fe]')

# Plot aquarius [Na/Fe] and [O/Fe] on caretta plot?

ax.scatter(O_Fe, Na_Fe, marker='o', color='b', s=50, zorder=10)
ax.errorbar(O_Fe, Na_Fe, xerr=e_O_Fe, yerr=e_Na_Fe, fmt=None, ecolor='b', elinewidth=1, zorder=10)

ax.set_xlim(-1.5, 1.0)
ax.set_ylim(-0.6, 1.1)

plt.draw()
plt.savefig('caretta-et-al-2009-o-na.pdf')
plt.savefig('caretta-et-al-2009-o-na.eps')

# Plot [Na/Fe] and [O/Fe] on its own

fig = plt.figure()
ax2 = fig.add_subplot(111)
ax2.scatter(O_Fe, Na_Fe, marker='o', color='k')
ax2.errorbar(O_Fe, Na_Fe, xerr=e_O_Fe, yerr=e_Na_Fe, fmt=None, ecolor='k', edgecolor='k')

ax2.set_xlabel('[O/Fe]')
ax2.set_ylabel('[Na/Fe]')

ax2.set_xlim(-1.5, 1.0)
ax2.set_ylim(-0.6, 1.1)

plt.draw()
plt.savefig('aquarius-o-na.pdf')
plt.savefig('aquarius-o-na.eps')
