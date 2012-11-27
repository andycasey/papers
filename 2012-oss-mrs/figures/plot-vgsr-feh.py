from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
labelsize, fontsize = 13, 13

# VGSR and [Fe/H] adopted as taken from table 1

data = [
	# Vgsr, Verr, [Fe/H]adopted, probability (1 = high, 2 = medium, 3 = low)
	(73.3, 9.3, -1.78, 3),
	(78.4, 5.2, -1.63, 1),
	(77.0, 4.0, -1.31, 3),
	(74.9, 17.6, -1.40, 1),
	(109.5, 9.0, -1.43, 2),
	(79.2, 3.3, -1.84, 1),
	(93.2, 29.8, -2.82, 1),
	(83.6, 3.5, -1.62, 1),
	(118.9, 11.7, -1.65, 1),
	(124.5, 6.7, -1.48, 3),
	(105.1, 5.1, 1.12, 3),
	(108.2, 9.0, -1.17, 1),
	(109.3, 8.1, -2.37, 2),
	(81.5, 4.6, -2.70, 1),
	(94.7, 5.1, -1.54, 3),
	(109.9, 25.4, -1.06, 3),
	(97.5, 5.9, -0.90, 2),
	(82.7, 5.0, -1.16, 1),
	(66.7, 8.7, -2.10, 2),
	(65.3, 5.4, -1.74, 2)
]

data = np.array(data)
col_vgsr, col_verr, col_feh, col_prob = range(4)

fig = plt.figure()
ax = fig.add_subplot(111)

edgecolors = np.array(['k', '#666666', 'r'])
facecolors = np.array(['k', '#ffffff'])
indices = np.arange(len(data))
indices = np.where(data[:, col_prob] != 3)[0]


idx = np.where(data[:, col_prob] == 1)[0]
ax.errorbar(data[idx, col_vgsr].flatten(), data[idx, col_feh].flatten(), ecolor=edgecolors[0], yerr=0.20, xerr=data[idx, col_verr].flatten(), fmt=None, zorder=-1)

idx = np.where(data[:, col_prob] == 2)
ax.errorbar(data[idx, col_vgsr].flatten(), data[idx, col_feh].flatten(), ecolor=edgecolors[1], yerr=0.20, xerr=data[idx, col_verr].flatten(), fmt=None, zorder=-1)

ax.scatter(data[indices, col_vgsr], data[indices, col_feh], marker='o', facecolor=facecolors[map(int, data[indices, col_prob] - 1)], edgecolor=edgecolors[map(int, data[indices, col_prob] - 1)], zorder=3, s=40)

#idx = np.where(data[:, col_prob] == 3)[0]
#ax.errorbar(data[idx, col_vgsr].flatten(), data[idx, col_feh].flatten(), ecolor=colors[2], yerr=0.20, xerr=data[idx, col_verr].flatten(), fmt=None, zorder=-1)

ax.set_xlabel('$V_{GSR}$ [km s$^{-1}$]', fontsize=labelsize)
ax.set_ylabel('[Fe/H] (adopted)', fontsize=labelsize)

plt.savefig('vgsr-feh.eps')
plt.savefig('vgsr-feh.pdf')


fig = plt.figure()
ax = fig.add_subplot(111)

newberg_bins = np.arange(-1.275, -3.1, -0.125)[::-1]
newberg_vals_in_bins = [1, 1, 1, 0, 1, 6, 6, 7, 5, 2, 3, 2, 1, 1] #[1, 1, 2, 3, 2, 5, 7, 6, 6, 1, 0, 1, 1, 1]

newberg_values = []

for i, val in enumerate(newberg_bins[:-1] + 0.5 * np.diff(newberg_bins)[0]):
	
	for j in xrange(newberg_vals_in_bins[i]):
		newberg_values.append(val)

newberg_bins = np.arange(-0.775, -3.1, -0.125)[::-1]

indices = np.where(data[:, col_prob] == 1)[0]
ax.hist(data[indices, col_feh], bins=newberg_bins, edgecolor='k', facecolor='#666666', lw=1, label='High probability candidates (This work)')

indices = np.where(data[:, col_prob] == 2)[0]
ax.hist(data[indices, col_feh], bins=newberg_bins, edgecolor='k', hatch='//', facecolor='none', lw=1, label='Medium probability candidates (This work)')

#indices = np.where(data[:, col_prob] == 2)[0]
#ax.hist(data[indices, col_feh], bins=newberg_bins, edgecolor='k', facecolor='#333333', lw=1)

ax.hist(newberg_values, bins=newberg_bins, facecolor='none', edgecolor='k', lw=2, label='BHB stars, Newberg et al. (2010)')

ax.set_xlim(-3.2, -0.7)
ax.set_xlabel('[Fe/H]')
ax.set_ylabel('$N$')
ax.set_ylim(0, 10)
ax.legend()

plt.savefig('newberg-feh.eps')
plt.savefig('newberg-feh.pdf')


