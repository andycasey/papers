from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


good_candidate_oids = [570, 169, 143, 337]# 571, 76, 264]
# removed "OSS-9" or oid 604

# Using isochrone consistency check
good_candidate_oids.extend([917, 768, 76, 685, 932])

# OSS 1
#good_candidate_oids.append(492)

# OSS 5
#good_candidate_oids.append(672)

#upperlimit_candidate_oids = [798, 570, 932, 604, 424, 672, 685]
upperlimit_candidate_oids = []
dodgy_candidate_oid = False

distance = 22.5 #kpc (my measurement)
my_distance_modulus = 5 * np.log10(distance * 1000) - 5
newberg_distance_modulus = 5 * np.log10(21.4 * 1000) - 5 # Newberg et al 2010 estimate distance at (l, b) = (250, 50) of 21.4 +/- 1.0 kpc

labelsize = 13
fontsize = 13

xlims = (0.1, 1.1)
ylims = (15.5, 20.5)


my_cmap = matplotlib.colors.LinearSegmentedColormap.from_list('mycm', ['#000080','#44BB44', '#AA0000'])


fig = plt.figure()
ax = fig.add_subplot(111)

# Load in all colour data
oid, g0, g, r0 = np.loadtxt('ew-mg-all.data', usecols=(0, 19, 9, 20, ), unpack=True)
ax.scatter(g0 - r0, g0, c='#cccccc', facecolor='none', alpha=0.2)


if dodgy_candidate_oid:
    idx = np.where(oid == dodgy_candidate_oid)[0]
    if len(idx) == 0:
        
        print "Couldn't find the dodgy candidate ID info!"
        #raise a
    
    elif len(idx) > 1:
        print "OKAY we have too many,.."
        raise a
    
    print idx, oid[idx], g0[idx] - r0[idx], g0[idx]
    ax.scatter(g0[idx] - r0[idx], g0[idx], facecolor='white', edgecolor='k', s=50)
    ax.scatter(g0[idx] - r0[idx], g0[idx], facecolor='none', edgecolor='k', s=130, zorder=-1)

# Load in isochrone information

g, r = np.loadtxt('giardi-10gyr-1.5.data', unpack=True)
ax.plot(g - r, g + newberg_distance_modulus, 'k:', zorder=-10, label='[Fe/H] $=-1.5$ at $21.4$ kpc')

g, r = np.loadtxt('giardi+10gyr-1.63.data', unpack=True)
ax.plot(g - r, g + my_distance_modulus, 'k-', lw=2, zorder=-1, label='[Fe/H] $=-1.71$ at $24$ kpc')

g, r = np.loadtxt('giardi-10gyr-2.0.data', unpack=True)
ax.plot(g - r, g + newberg_distance_modulus, 'k--', zorder=-10, label='[Fe/H] $=-2.0$ at $21.4$ kpc')



# Load in our candidates
#oid, g0, r0, ew, feh_cat, feh_iso, feh_adpoted, feh_iso_upperlimit = np.loadtxt('all-candidates.data', usecols=(0, 19, 20, 5, 25, 26, 27, 28)
cols = (0, 1, 2, 19, 20, 4, 3, 5, 25, 26, 27, 28)

oid, ra, dec, g0, r0, vgsr, verr, ew, feh_cat, feh_iso, feh_adopted, feh_iso_upperlimit = np.loadtxt('all-candidates.data', usecols=cols, unpack=True)

table_data = np.loadtxt('all-candidates.data', usecols=cols)

# Plot upper measurements
idx = np.where(ew < 0)
scat_upperlims = ax.scatter(g0[idx] - r0[idx], g0[idx], c=feh_cat[idx], marker='v', edgecolor='k', vmin=-2.8, vmax=-0.8, s=50, cmap=my_cmap)


# Plot normal measurements
idx = np.where(ew > 0)
scat = ax.scatter(g0[idx] - r0[idx], g0[idx], c=feh_cat[idx], edgecolor='k', vmin=-2.8, vmax=-0.8, s=50, cmap=my_cmap)


for g__, gmr__, feh__ in zip(g0[idx], g0[idx] - r0[idx], feh_cat[idx]):
    print g__, gmr__, feh__



gc_gr0 = []
gc_g0 = []
print "OID RA DEC G0 R0 VGSR VERR EW FEH_CAT FEH_ISO FEH_ADOPTED UPPERLIM"
for candidate_oid in good_candidate_oids:
    idx = np.where(oid == candidate_oid)[0]
    
    if len(idx) == 0:
        print "Couldn't find a candidate,..."
        raise a
    
    elif len(idx) > 1:
        print "We have too many good candidates"
        raise a
    
    gc_gr0.append(g0[idx] - r0[idx])
    gc_g0.append(g0[idx])
    
    print "%i   %3.5f   %3.5f   %2.2f   %2.2f   %3.1f   %3.1f   %3.4f   %2.2f   %2.2f   %2.2f  %i" % tuple(table_data[idx][0])
    
ax.scatter(gc_gr0, gc_g0, marker='o', facecolor='none', edgecolor='k', s=130, zorder=-1)

gc_gr0 = []
gc_g0 = []
for candidate_oid in upperlimit_candidate_oids:
    idx = np.where(oid == candidate_oid)[0]

    if len(idx) == 0:
        print "Couldn't find an upper limit candidate"
        raise a
    
    elif len(idx) > 1:
        print "We have too many upper limit matches!"
        raise a
    
    gc_gr0.append(g0[idx] - r0[idx])
    gc_g0.append(g0[idx])

    print "%i   %3.5f   %3.5f   %2.2f   %2.2f   %3.1f   %3.1f   %3.4f   %2.2f   %2.2f   %i" % tuple(table_data[idx][0])

ax.scatter(gc_gr0, gc_g0, marker='v', facecolor='none', edgecolor='k', s=130, zorder=-1)

cbar = plt.colorbar(scat)

# Set labels
cbar.set_label('[Fe/H] (adopted)', fontsize=labelsize)
ax.set_xlabel('$g - r$', fontsize=labelsize)
ax.set_ylabel('$g$', fontsize=labelsize)


# Set limits
ax.set_xlim(xlims)
ax.set_ylim(ylims)

# Invert y axis
ax.set_ylim(ax.get_ylim()[::-1])

    
plt.savefig('cmd.pdf')
plt.close(fig)

# Plot the [Fe/H]_CaT vs [Fe/H]_iso

fig = plt.figure()
ax = fig.add_subplot(111)

# plot the normal ones first
idx = np.where(feh_iso_upperlimit == 0)[0]
ax.scatter(feh_cat[idx], feh_iso[idx], marker='o', facecolor='k')

# now plot the lower limits
idx = np.where(feh_iso_upperlimit == 1)[0]
print "idx ", idx
ax.errorbar(feh_cat[idx], feh_iso[idx], fmt='ok', lolims=True, yerr=0.1)

intersect = -2.25
minimum_isochrone_metallicity = -2.28
pm_window = 0.30
ax.plot([-3.0, -0.5], [-3.0, -0.5], 'k-')
ax.plot([-3.0, -0.5], [-3.0 - pm_window, -0.5 - pm_window], 'k:', zorder=-1)
ax.plot([-3.0, intersect - pm_window, -0.5], [intersect, intersect, -0.5 + pm_window], 'k:', zorder=-1)

ax.fill_between([-3.0, intersect - pm_window, -0.50], [-3.0 - pm_window, intersect - pm_window*2, -0.50 - pm_window], [intersect, intersect, -0.5 + pm_window], facecolor='g', alpha=0.5, zorder=-2)

ax.text(-0.6, -2.5, "Minimum isochrone metallicity", color='#111111', horizontalalignment="right", fontsize=labelsize)
ax.plot([-3.0, -0.5], [minimum_isochrone_metallicity, minimum_isochrone_metallicity], '--', color='#111111', zorder=-1)


# Add a representative error 
ax.errorbar([-2.70], [-0.80], color='k', xerr=0.2, yerr=0.2)

ax.set_xlabel(r'[Fe/H] (Ca II lines)', fontsize=labelsize)
ax.set_ylabel(r'[Fe/H] (Isochrone fitting)', fontsize=labelsize)

ax.set_xlim([-3.0, -0.5])
ax.set_ylim([-3.0, -0.5])
plt.draw()

plt.savefig('feh.pdf')
plt.close(fig)

