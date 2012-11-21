import matplotlib.pyplot as plt

#NGC;ID;Bmag;Vmag;Icmag;Kmag;RV11;RV13;Teff;[FeI/H];e_[FeI/H];[O/Fe];e_[O/Fe];l_[Na/Fe];[Na/Fe];e_[Na/Fe];RAJ2000;DEJ2000
# ; ;mag;mag;mag;mag;km/s;km/s;K;[Sun];[Sun];[Sun];[Sun]; ;[Sun];[Sun];"h:m:s";"d:m:s"


Fe_I, e_Fe_I, O_Fe, e_O_Fe, l_Na_Fe, Na_Fe, e_Na_Fe = np.loadtxt('caretta-et-al-2009-O-Na.txt', skiprows=67, delimiter=';', \
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


Fe_I, e_Fe_I, O_Fe, e_O_Fe, Na_Fe, e_Na_Fe = map(clean, [Fe_I, e_Fe_I, O_Fe, e_O_Fe, Na_Fe, e_Na_Fe])



fig = plt.figure()
ax = fig.add_subplot(111)

# Oxygen upper limits
upper_limits = np.where(e_O_Fe == 9.999)[0]
measurements = np.array(list(set(range(len(O_Fe))).difference(upper_limits)))

ax.errorbar(O_Fe[upper_limits], Na_Fe[upper_limits], xerr=0.02, xlolims=True, fmt=None, ecolor='k')
#ax.errorbar(O_Fe[upper_limits], Na_Fe[upper_limits], facecolor='none', edgecolor='k')
#ax.errorbar(O_Fe[measurements], Na_Fe[measurements], xerr=e_O_Fe[measurements], yerr=e_Na_Fe[measurements], fmt=None)

ax.scatter(O_Fe[measurements], Na_Fe[measurements], facecolor='none', edgecolor='k')

ax.set_xlabel('[O/Fe]')
ax.set_ylabel('[Na/Fe]')

ax.set_xlim(-1.5, 1.0)
ax.set_ylim(-0.6, 1.1)

plt.draw()
plt.savefig('caretta-et-al-2009-o-na.pdf')

