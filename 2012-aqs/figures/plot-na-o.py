import matplotlib.pyplot as plt

#NGC;ID;Bmag;Vmag;Icmag;Kmag;RV11;RV13;Teff;[FeI/H];e_[FeI/H];[O/Fe];e_[O/Fe];l_[Na/Fe];[Na/Fe];e_[Na/Fe];RAJ2000;DEJ2000
# ; ;mag;mag;mag;mag;km/s;km/s;K;[Sun];[Sun];[Sun];[Sun]; ;[Sun];[Sun];"h:m:s";"d:m:s"


caretta_cluster, caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_l_Na_Fe, caretta_Na_Fe, caretta_e_Na_Fe = np.loadtxt('caretta-et-al-2009-O-Na.txt', skiprows=67, delimiter=';', \
														usecols=(0, 9, 10, 11, 12, 13, 14, 15, ), unpack=True, dtype=str)

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


caretta_cluster, caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_Na_Fe, caretta_e_Na_Fe = map(clean, [caretta_cluster, caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_Na_Fe, caretta_e_Na_Fe])


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



# Subplots


# Mg-Al in sub-plots


# Delete 6397 because it has less stars than even our stream sample!
unique_clusters = np.unique(caretta_cluster)
idx = list(unique_clusters).index(6397.)
unique_clusters = np.delete(unique_clusters, idx)

yplots = 5
xplots = int(np.ceil(float(len(unique_clusters) + 1) / yplots))

xlims = (-0.9, 0.9)
ylims = (-0.6, 1.1)

fig = plt.figure()
fig.subplots_adjust(hspace=0, wspace=0,right=0.97,top=0.97,bottom=0.08,left=0.09)

for i, cluster in enumerate(unique_clusters):

	ax = fig.add_subplot(xplots, yplots, i + 1)

	index = np.where(caretta_cluster == cluster)[0]

	cluster_na_fe = np.array(caretta_Na_Fe)[index]
	cluster_o_fe = np.array(caretta_O_Fe)[index]

	ax.scatter(cluster_o_fe, cluster_na_fe, facecolor='k')

	#A = np.vstack([cluster_na_fe, np.ones(len(cluster_na_fe))]).T
	#m, c = np.linalg.lstsq(A, cluster_o_fe)[0]

	# Polyfit
	#p = np.poly1d(np.polyfit(cluster_na_fe, cluster_o_fe, 1))

	x_range = np.array([np.min(cluster_o_fe), np.max(cluster_o_fe)])
	y_range = np.array([np.min(cluster_na_fe), np.max(cluster_na_fe)])

	#ax.plot(x_range, m * x_range + c, '-', color='#666666')
	#ax.plot(x_range, p(x_range), '-', color='b')

	ax.set_xlim(xlims)
	ax.set_xticks([-0.5, 0.0, 0.5])

	ax.set_ylim(ylims)
	ax.set_yticks([-0.5, 0.0, 0.5, 1.0])


	ax.set_xlabel('[O/Fe]')
	ax.set_ylabel('[Na/Fe]')

	if i % (yplots): 
		ax.set_ylabel('')
		ax.set_yticklabels([''] * len(ax.get_yticklabels()))
	
	if len(unique_clusters) - i - 1 > xplots:

		ax.set_xlabel('')
		ax.set_xticklabels([''] * len(ax.get_xticklabels()))

		ax.xaxis.set_visible(False)
	
	ax.text(xlims[0] + 0.05 * (xlims[1] - xlims[0]), ylims[0] + 0.05 * (ylims[1] - ylims[0]), 'NGC %i' % (cluster, ), verticalalignment='bottom')


# Add 
ax = fig.add_subplot(xplots, yplots, i + 2)
ax.scatter(O_Fe, Na_Fe, facecolor='k')

A = np.vstack([O_Fe, np.ones(len(O_Fe))]).T
m, c = np.linalg.lstsq(A, Na_Fe)[0]

print 'Aquarius', m

x_range = np.array([np.min(O_Fe), np.max(O_Fe)])
y_range = np.array([np.min(Na_Fe), np.max(Na_Fe)])

ax.plot(x_range, m * x_range + c, '-', color='#666666')

ax.set_xlabel('[O/Fe]')
ax.set_ylabel('[Na/Fe]')

ax.set_xlim(xlims)
ax.set_xticks([-0.5, 0.0, 0.5])

ax.set_ylim(ylims)
ax.set_yticks([-0.5, 0.0, 0.5, 1.0])




if i + 1 % yplots:
	ax.set_ylabel('')
	ax.set_yticklabels([''] * len(ax.get_yticklabels()))
	
ax.text(xlims[0] + 0.05 * (xlims[1] - xlims[0]), ylims[0] + 0.05 * (ylims[1] - ylims[0]), 'Aquarius', verticalalignment='bottom')


plt.savefig('aquarius-o-na-cluster.pdf')



# Get data from Reddy et al 2005
reddy_data = np.loadtxt('Reddy_thick_disk_2005.data', delimiter='|', usecols=(1, 18, 21, 20), dtype=str)

Reddy_ThickDisk_Na_Fe = []
Reddy_ThickDisk_O_Fe = []
Reddy_ThinDisk_Na_Fe = []
Reddy_ThinDisk_O_Fe = []
Reddy_ThickThinDisk_Na_Fe = []
Reddy_ThickThinDisk_O_Fe = []
Reddy_Halo_Na_Fe = []
Reddy_Halo_O_Fe = []
Reddy_ThickDisk_Halo_Na_Fe = []
Reddy_ThickDisk_Halo_O_Fe = []

for classification, fe_h, na_h, o_h in reddy_data:
	
	classification, fe_h, na_h, o_h = [var.strip() for var in (classification, fe_h, na_h, o_h)]


	_Na_Fe, _O_Fe = np.nan, np.nan

	if len(fe_h) == 0:
		fe_h, na_h, o_h = [np.nan] * 3

	else:

		if len(na_h) == 0:
			na_h = np.nan

		else:
			_Na_Fe = float(na_h) #- float(fe_h)


		if len(o_h) == 0:
			o_h = np.nan

		else:
			_O_Fe = float(o_h) #- float(fe_h)


	if classification == 'Thick disc':
		Reddy_ThickDisk_Na_Fe.append(_Na_Fe)
		Reddy_ThickDisk_O_Fe.append(_O_Fe)

	elif classification == 'Thin disc':
		Reddy_ThinDisk_Na_Fe.append(_Na_Fe)
		Reddy_ThinDisk_O_Fe.append(_O_Fe)

	elif classification == 'Halo':
		Reddy_Halo_Na_Fe.append(_Na_Fe)
		Reddy_Halo_O_Fe.append(_O_Fe)

	elif classification == 'Thin/thick disc':
		Reddy_ThickThinDisk_Na_Fe.append(_Na_Fe)
		Reddy_ThickThinDisk_O_Fe.append(_O_Fe)

	elif classification == 'Thick disc/halo':
		Reddy_ThickDisk_Halo_Na_Fe.append(_Na_Fe)
		Reddy_ThickDisk_Halo_O_Fe.append(_O_Fe)



Reddy_ThickDisk_O_Fe, Reddy_ThickDisk_Na_Fe, Reddy_ThinDisk_O_Fe, Reddy_ThinDisk_Na_Fe, Reddy_Halo_O_Fe, Reddy_Halo_Na_Fe, Reddy_ThickThinDisk_O_Fe, Reddy_ThickThinDisk_Na_Fe, Reddy_ThickDisk_Halo_O_Fe, Reddy_ThickDisk_Halo_Na_Fe = [np.array(item) for item in (Reddy_ThickDisk_O_Fe, Reddy_ThickDisk_Na_Fe, Reddy_ThinDisk_O_Fe, Reddy_ThinDisk_Na_Fe, Reddy_Halo_O_Fe, Reddy_Halo_Na_Fe, Reddy_ThickThinDisk_O_Fe, Reddy_ThickThinDisk_Na_Fe, Reddy_ThickDisk_Halo_O_Fe, Reddy_ThickDisk_Halo_Na_Fe)]

fig = plt.figure()
fig.subplots_adjust(hspace=0, wspace=0,right=0.97,top=0.97,bottom=0.08,left=0.10)
ax = fig.add_subplot(111)

ax.scatter(Reddy_ThickDisk_O_Fe, Reddy_ThickDisk_Na_Fe, marker='+', facecolor='none', linewidth=1, edgecolor='k', label='Thick disc (Reddy+ 2005)')
ax.scatter(Reddy_ThinDisk_O_Fe, Reddy_ThinDisk_Na_Fe, marker='x', facecolor='none', linewidth=1, edgecolor='k', label='Thin disc (Reddy+ 2005)')
ax.scatter(Reddy_Halo_O_Fe, Reddy_Halo_Na_Fe, marker='o', facecolor='none', edgecolor='k', label='Halo (Reddy+ 2005)')
ax.scatter(Reddy_ThickThinDisk_O_Fe, Reddy_ThickThinDisk_Na_Fe, marker='+', facecolor='none', edgecolor='k')
ax.scatter(Reddy_ThickDisk_Halo_O_Fe, Reddy_ThickDisk_Halo_Na_Fe, marker='o', facecolor='none', edgecolor='k')

#caretta_cluster, caretta_Fe_I, caretta_e_Fe_I, caretta_O_Fe, caretta_e_O_Fe, caretta_Na_Fe, caretta_e_Na_Fe

ax.scatter(caretta_O_Fe, caretta_Na_Fe, marker='+', edgecolor='#cccccc', facecolor='none', zorder=-1, label='GC stars (Carretta+ 2009)')

ax.scatter(O_Fe, Na_Fe, marker='o', facecolor='k', s=40)
ax.errorbar(O_Fe, Na_Fe, xerr=e_O_Fe, yerr=e_Na_Fe, fmt=None, ecolor='k')

ax.legend(loc=2, prop={'size': 11})

ax.set_xlabel('[O/Fe]')
ax.set_ylabel('[Na/Fe]')

ax.set_xlim(xlims)
ax.set_ylim(ylims)

plt.savefig('aquarius-o-na-halo.pdf')
plt.savefig('aquarius-o-na-halo.eps')