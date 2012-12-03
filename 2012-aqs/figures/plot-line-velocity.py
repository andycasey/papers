import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm

rest_wavelengths, line_velocities, is_acceptable = np.loadtxt('HD41667-line-measurements.data', delimiter=';', usecols=(0, 11, -2), unpack=True, dtype=str)

# Get only the acceptable indices

acceptable_indices = np.array([item.strip() for item in is_acceptable]) == 'True'
rest_wavelengths = np.array(rest_wavelengths[acceptable_indices], dtype=np.float)
line_velocities = np.array(line_velocities[acceptable_indices], dtype=np.float)

# Get only finite line velocities
finite_indices = np.where(np.isfinite(line_velocities))
line_velocities = line_velocities[finite_indices]
rest_wavelengths = rest_wavelengths[finite_indices]



fig = plt.figure(figsize=(6,6))
fig.subplots_adjust(hspace=0.25)

# histogram 
ax = fig.add_subplot(211)
n, bins, patches = ax.hist(line_velocities, facecolor='none', normed=True)
ax.set_xlabel('Line velocity, $V_{\lambda}$ [km s$^{-1}$]')
ax.set_ylabel('$P(V_{line})$')


mu, sigma = norm.fit(line_velocities)
print "mu, sigma", mu, sigma
y = mlab.normpdf(bins, mu, sigma)
ax.plot(bins, y, 'k:', linewidth=2)

ax.text(3.0, 0.7, "$N\,=\,%i$\n$\mu\,=\,%1.2f$\n$\sigma\,=\,%1.2f$" % (len(line_velocities), mu, sigma, ), color='k', verticalalignment='top', horizontalalignment='right')

# scatter
ax2 = fig.add_subplot(212)
ax2.scatter(rest_wavelengths, line_velocities, facecolor='none')

xlim = ax2.get_xlim()
ylim = int(np.max(np.abs(ax2.get_ylim()))) - 0.5

ax2.plot(xlim, [0, 0], 'k-', zorder=-1)
ax2.plot(xlim, [np.median(line_velocities)] * 2, 'k--', zorder=-1)

ax2.set_xlim(xlim)
ax2.set_ylim(-ylim, ylim)
ax.set_xlim(-ylim, ylim)

ax2.set_xlabel('Rest wavelength, $\lambda$ [$\AA{}$]')
ax2.set_ylabel('Line velocity, $V_{\lambda}$ [km s$^{-1}$]')

plt.savefig('line-velocity.pdf')
plt.savefig('line-velocity.eps')

