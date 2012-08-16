from scipy.interpolate import splrep, splev

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

fontsize = 13
labelsize = 13

xlims = (0.4, 1.2)
ylims = (0, 1.5)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.set_xlim(xlims)
ax.set_ylim(ylims)

# Plot our data
g, r, ew = np.loadtxt('ew-mg-all.data', usecols=(19, 20, 5, ), unpack=True)

# Distinguish upper limits and normal measurements
idx = np.where(ew < 0)
g_upperlims, r_upperlims, ew_upperlims = (item.copy() for item in (g[idx], r[idx], ew[idx]))

idx = np.where(ew > 0)
g, r, ew = g[idx], r[idx], ew[idx]

# Plot everything

# Now plot the upper limits
ax.errorbar(g_upperlims - r_upperlims, -ew_upperlims, fmt='ko', yerr=0.05, lolims=True, zorder=0)
ax.scatter(g_upperlims - r_upperlims, -ew_upperlims, edgecolor='none', marker='o', facecolor='w', zorder=1)

ax.scatter(g - r, ew, marker='o', facecolor='none')


# Now plot the dividing spline

match_points = np.array([
    [0.40, 0.35],
    [0.50, 0.45],
    [0.70, 0.72],
    [0.80, 0.80],
    [0.85, 0.83],
    [0.90, 0.84],
    [1.00, 0.81],
    [1.10, 0.75],
    [1.20, 0.70]
])
# convert to mA

tck = splrep(match_points[:, 0], match_points[:, 1], k=3)

x = np.arange(ax.get_xlim()[0], ax.get_xlim()[1], 0.01)
y = splev(x, tck)

giants = 0
for gi, ri, ewi in zip(g, r, ew):
    idx = np.searchsorted(x, gi - ri) - 1
    if idx >= len(y):
        continue
    
    if y[idx] >= ewi:
        giants += 1
        
print "Total number", len(ew)
print "Giants:", giants, giants/float(len(ew)) * 100.
print "Dwarfs:", len(ew) - giants, (len(ew) - giants)/float(len(ew)) * 100.

#ax.scatter(match_points[:, 0], match_points[:, 1], facecolor='r', s=100)
ax.plot(x, y, 'g', alpha=0.80, lw=5)
#ax.fill_between(x, 0, y, alpha=0.5, facecolor='g')

ax.set_xlim((xlims[0], xlims[1] - 0.05))
ax.set_ylim(ylims)

ax.set_ylabel('$EW_{\lambda8807}$ [A]', fontsize=labelsize)
ax.set_xlabel('$g - r$', fontsize=labelsize)
plt.savefig('ew-mg.pdf')