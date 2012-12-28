import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

wavelength, smh_ew, norris_ew = np.loadtxt('SMH-Norris-comparison.data', usecols=(0, 1, 2, ), unpack=True)

fig = plt.figure(figsize=(6,7))
fig.subplots_adjust(hspace=0.0, wspace=0.0)
gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2])

ax1 = fig.add_subplot(gs[0])
#ax1 = plt.subplot2grid((3, 1), (0, 0))

ax1.scatter(smh_ew, smh_ew - norris_ew, facecolor='none', edgecolor='k', marker='+')
ax1.plot([0, 200], [0, 0], 'k-', zorder=-1)

A = np.vstack([smh_ew, np.ones(len(norris_ew))]).T
m, c = np.linalg.lstsq(A, smh_ew - norris_ew)[0]
x = np.array([np.min(smh_ew), np.max(smh_ew)])
ax1.plot(x, m * x + c, 'k:')

ylim = np.max(np.abs(np.array(ax1.get_ylim())))
ax1.set_ylim(-15, 15)

ax1.xaxis.set_visible(False)
ax1.set_ylabel('$\Delta{}W_\lambda$ [m$\AA{}$]')


ax2 = fig.add_subplot(gs[1], sharex=ax1)
#ax2 = plt.subplot2grid((3, 1), (1, 0), rowspan=2)
ax2.scatter(smh_ew, norris_ew, facecolor='none', edgecolor='k', marker='+')

A = np.vstack([norris_ew, np.ones(len(norris_ew))]).T
m, c = np.linalg.lstsq(A, smh_ew)[0]
x = np.array([0, 200])
ax2.plot(x, x, 'k-', zorder=-1)
x = np.array([np.min(smh_ew), np.max(smh_ew)])
ax2.plot(x, m * x + c, 'k:')

# Plot an error cone
error = 10 # percent
bounds = np.array([0, 160])
#ax2.plot(bounds, bounds * (1 + error/100.), '-', c='#aaaaaa', zorder=-5)
#ax2.plot(bounds, bounds * (1 - error/100.), '-', c='#aaaaaa', zorder=-5)


ax1.set_xlim(bounds)
ax2.set_xlim(bounds)
ax2.set_ylim(bounds)

ax2.set_xlabel('$W_\lambda$ (This work, automatic) [m$\AA{}$]')
ax2.set_ylabel('$W_\lambda$ (Norris et al. 1996) [m$\AA{}$]')

ax2.get_yticklabels()[-1].set_visible(False)
ax1.get_yticklabels()[0].set_visible(False)
ax1.get_yticklabels()[-1].set_visible(False)

ax1.text(5, 10, '$\langle{}\Delta{}W_\lambda\\rangle{}\,=\,-0.64\,\pm\,2.78\,$m${\AA}$', color='k', verticalalignment='center')

ax2.text(5, 150, "$a_0\,=\,%1.2f$\n$a_1\,=\,%1.2f$" % (c, m, ), verticalalignment='top')

ax1.set_title('%i lines in HD 140283' % (len(smh_ew), ))

plt.savefig('smh-norris.pdf')
plt.savefig('smh-norris.eps')

