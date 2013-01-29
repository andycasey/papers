

# C2225316, C2306265, J221821, J223504, J223811
Vhels = [-156.4, -221.1, -159.5, -169.7, -235.7]
Verrs = [0.1, 0.1, 0.1, 0.1, 0.1]

Fe_H = [-1.26, -1.17, -1.61, -0.68, -1.45]
Fe_H_errs = [0.01, 0.01, 0.01, 0.02, 0.01]
Fe_H_errs = [0.10, 0.1, 0.1, 0.1, 0.1]


# Wylie de boer data
# C2225316, C2306265, J221821, J223504, J223811

wdb_Vhels = [-155.7, -221.8, -154.1, -166.9, -230.1]
wdb_Verrs = [0.7, 1.7, 1.1, 1.3, 1.9]

wdb_Fe_H = [-1.20, np.NaN, -1.15, -0.98, -1.20]
wdb_Fe_H_errs = [0.14, np.NaN, 0.21, 0.17, 0.20]


# Williams data
# C2225316, C2306265, J221821, J223504, J223811
williams_Vhels = [-155.7, -221.8, -154.1, -166.9, -230.1]
williams_Verrs = [0.7, 1.7, 1.1, 1.3, 1.9]

williams_Fe_H = [-1.29, np.NaN, -1.54, -0.33, -0.78]
williams_Fe_H_errs = [0.2, 0.2, 0.2, 0.2, 0.2]


xlims = (-245, -145)
ylims = (-1.8, -0.2)

print "NOT READING IN METALLICITY OR VELOCITY FROM TABLES: HARD CODED IN FOR FEH-VHEL"
fig = plt.figure()

ax = fig.add_subplot(111)

s = 50
ax.scatter(williams_Vhels, williams_Fe_H, facecolor='g', edgecolor='g', marker='s', s=s, label='Williams et al. (2011)')
ax.scatter(wdb_Vhels, wdb_Fe_H, facecolor='b', edgecolor='b', marker='^', s=s, label='Wylie-de Boer et al. (2012)')
ax.scatter(Vhels, Fe_H, c='k', marker='o', s=s, label='This work')

ax.errorbar(wdb_Vhels, wdb_Fe_H, xerr=wdb_Verrs, yerr=wdb_Fe_H_errs, fmt=None, color='b', ecolor='b')
ax.errorbar(williams_Vhels, williams_Fe_H, xerr=williams_Verrs, yerr=williams_Fe_H_errs, fmt=None, color='g', ecolor='g')


ax.errorbar(Vhels, Fe_H, yerr=Fe_H_errs, xerr=Verrs, fmt=None, color='k', ecolor='k')

#ax.errorbar(xlims[0] + 0.05 * np.abs(xlims[1] - xlims[0]), ylims[1] - 0.1 * np.abs(ylims[1] - ylims[0]), 0.1, 1.0, edgecolor='k', fmt=None)

# Connect the lines
for i in xrange(len(Vhels)):
	data_x = [Vhels[i], wdb_Vhels[i], williams_Vhels[i]]
	data_y = [Fe_H[i], wdb_Fe_H[i], williams_Fe_H[i]]

	# Full triangle
	data_x.append(data_x[0])
	data_y.append(data_y[0])

	ax.plot(data_x, data_y, '--', c='k', zorder=-1)

ax.legend(loc=2)

ax.set_xlim(xlims)
ax.set_ylim(ylims)

ax.set_xlabel('$V_{\\rm hel}$ (km s$^{-1}$)')
ax.set_ylabel('[Fe/H]')

plt.draw()

plt.savefig('aquarius-vhel-feh.pdf')
plt.savefig('aquarius-vhel-feh.eps')