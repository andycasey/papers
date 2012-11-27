 # -*- coding: utf-8 -*-

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)



fontsize = 11
labelsize = 11
binx = 20.

oss_vgsrs = [65, 125]
oss_member_count = 10
oss_vgsr_centroid = oss_vgsrs[0] + 0.5 * (oss_vgsrs[1] - oss_vgsrs[0])




matplotlib.rc('font', **{'sans-serif' : 'Arial',
                         'family' : 'sans-serif'})

obs_vgsr = np.loadtxt('velocities.data', delimiter=',', usecols=(-1,))
besancon_vgsr = np.loadtxt('besancon_sel.txt', usecols=(-2,))


bins = np.arange(-250, 250, binx)

fig = plt.figure()
obs_ax = fig.add_subplot(111)

hist_obs, bins = np.histogram(obs_vgsr, bins=bins)



# Normalise the besancon results
hist_besancon, bins = np.histogram(besancon_vgsr, bins=bins, normed=True)
hist_besancon *= len(besancon_vgsr) * binx/10. *4.

#obs_ax.hist([obs_vgsr, besancon_vgsr], bins=bins, histtype='bar', color=['grey', 'none'], label=['Observed', u'Besançon'])
#obs_ax.bar([bins[:-1], bins[:-1]], [hist_obs, hist_besancon], [np.diff(bins), np.diff(bins)], color=['r', 'g'])



#obs_ax.bar(bins[:-1] + np.diff(bins)[0]/2, hist_besancon, np.diff(bins)/2, color='white', label=u'Besançon (Robin et al. 2003)')
#obs_ax.scatter(bins[:-1] + np.diff(bins)[0]/2, hist_besancon, facecolor='k', zorder=3)
obs_bar = obs_ax.bar(bins[:-1], hist_obs, np.diff(bins), color='#cccccc', label='Observed')
obs_plot = obs_ax.plot(bins[:-1] + np.diff(bins)[0]/2, hist_besancon, '--ok', zorder=2, label=u'Besançon model (Robin et al. 2003)')


obs_ax.set_xlim(-250, 250)
obs_ax.set_ylim(0, obs_ax.get_ylim()[1])
obs_ax.set_xlabel('$V_{GSR}$ [km s$^{-1}$]', fontsize=labelsize)
obs_ax.set_ylabel('$N$', fontsize=labelsize)

patch = Rectangle((oss_vgsrs[0], 0), oss_vgsrs[1] - oss_vgsrs[0], obs_ax.get_ylim()[1], color='#85c078', zorder=-1)
obs_ax.add_patch(patch)

idx = np.searchsorted(bins, oss_vgsr_centroid) - 1
obs_ax.annotate('Orphan Stream', (oss_vgsr_centroid, hist_obs[idx] + 5), (oss_vgsr_centroid, obs_ax.get_ylim()[1]/4), horizontalalignment='center', arrowprops={'arrowstyle': '->'}, fontsize=13)

handles, labels = obs_ax.get_legend_handles_labels()
plt.legend(handles[::-1], labels[::-1])

plt.draw()
plt.show()
                                
plt.savefig('vgsr-histogram.pdf')
plt.savefig('vgsr-histogram.eps')
