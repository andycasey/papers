import matplotlib.pyplot as plt
import numpy as np


# Caretta et al. 2009 GIRAFFE DATA

Carretta_Cluster, Carretta_Fe_H, Carretta_O_Fe, Carretta_Na_Fe, Carretta_Na_limit = np.loadtxt('Carretta_2009-O-Na.csv', delimiter='|',
	usecols=(0, 1, 3, 6, 5), dtype=str, unpack=True)
#Carretta_O_Fe, Carretta_Na_Fe, Carretta_Cluster, Carretta_Fe_H = np.loadtxt('Carretta_2009a.data', delimiter='|', usecols=(-3, -2, 2, -7), dtype=str, unpack=True)

fontsize=8.

delete_rows = []
for i in xrange(len(Carretta_Cluster)):

	O_Fe = Carretta_O_Fe[i].strip()
	Na_Fe = Carretta_Na_Fe[i].strip()

	if len(O_Fe) * len(Na_Fe) == 0:
		delete_rows.append(i)


Carretta_Cluster = map(float, np.delete(Carretta_Cluster, delete_rows))
Carretta_Fe_H = map(float, np.delete(Carretta_Fe_H, delete_rows))
Carretta_O_Fe = map(float, np.delete(Carretta_O_Fe, delete_rows))
Carretta_Na_limit = map(lambda x: x.count('<'), np.delete(Carretta_Na_limit, delete_rows)) # 1 = True, 0 = False
Carretta_Na_Fe = map(float, np.delete(Carretta_Na_Fe, delete_rows))

unique_clusters = np.unique(Carretta_Cluster)

data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

O_Fe = []
e_O_Fe = []

Na_Fe = []
e_Na_Fe = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('O')
	_O_Fe, _e_O_Fe = data[idx, 1:]

	idx = list(data[:, 0]).index('Na')
	_Na_Fe, _e_Na_Fe = data[idx, 1:]

	O_Fe.append(_O_Fe)
	Na_Fe.append(_Na_Fe)
	e_O_Fe.append(_e_O_Fe)
	e_Na_Fe.append(_e_Na_Fe)


# Clean up O_Fe, Na_Fe, e_O_Fe, e_Na_Fe
O_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in O_Fe])
Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Na_Fe])

e_O_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_O_Fe])
e_Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Na_Fe])

O_Fe, Na_Fe, e_O_Fe, e_Na_Fe = map(np.array, [O_Fe, Na_Fe, e_O_Fe, e_Na_Fe])


std_metallicities = []
mean_metallicities = []
for i, cluster in enumerate(unique_clusters):

	index = np.where(Carretta_Cluster == cluster)[0]

	mean_metallicity = np.mean(np.array(Carretta_Fe_H)[index])
	std_metallicity = np.std(np.array(Carretta_Fe_H)[index])

	print "%s [Fe/H] = %1.2f +/- %1.2f" % (cluster, mean_metallicity, std_metallicity, )

	mean_metallicities.append(mean_metallicity)
	std_metallicities.append(std_metallicity)


idx = np.argsort(mean_metallicities)[::-1]
unique_clusters = unique_clusters[idx]
mean_metallicities = np.array(mean_metallicities)[idx]
std_metallicities = np.array(std_metallicities)[idx]

# del 6397
idx = np.where(unique_clusters == 6397.)[0]
unique_clusters = np.delete(unique_clusters, idx)
mean_metallicities = np.delete(mean_metallicities, idx)
std_metallicities = np.delete(std_metallicities, idx)


fig = plt.figure()
fig.subplots_adjust(hspace=0, wspace=0,right=0.97,top=0.97,bottom=0.12,left=0.12)


yplots = 5
xplots = int(np.ceil(float(len(unique_clusters) + 1) / yplots))

xlims = (-0.9, 0.9)
ylims = (-0.65, 1.15)



for i, (cluster, mean_metallicity, std_metallicity) in enumerate(zip(unique_clusters, mean_metallicities, std_metallicities)):

	ax = fig.add_subplot(xplots, yplots, i + 1)

	index = np.where(Carretta_Cluster == cluster)[0]

	Cluster_O_Fe = np.array(Carretta_O_Fe)[index]
	Cluster_Na_Fe = np.array(Carretta_Na_Fe)[index]

	Cluster_Na_Limits = np.array(Carretta_Na_limit)[index] > 0


	if np.any(Cluster_Na_Limits):
		n = np.sum(Cluster_Na_Limits)
		ax.errorbar(Cluster_O_Fe[Cluster_Na_Limits], Cluster_Na_Fe[Cluster_Na_Limits],fmt=None, ecolor='#666666', solid_capstyle='round', xerr=[[0.15] * n, [0] * n], xlolims=True, zorder=-32)
	
		ax.scatter(Cluster_O_Fe[Cluster_Na_Limits], Cluster_Na_Fe[Cluster_Na_Limits], s=30, marker='o', edgecolor='#666666', facecolor='w', zorder=10)
	
	#ax.scatter(Cluster_O_Fe[Cluster_Na_Limits], Cluster_Na_Fe[Cluster_Na_Limits], marker='o', c='b')


	ax.scatter(Cluster_O_Fe[~Cluster_Na_Limits], Cluster_Na_Fe[~Cluster_Na_Limits], s=30, marker='o', edgecolor='k', facecolor='w', zorder=11)

	
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
	
	ax.text(xlims[0] + 0.05 * (xlims[1] - xlims[0]), ylims[0] + 0.03 * (ylims[1] - ylims[0]), 'NGC %i\n[Fe/H] = %1.2f$\pm$%1.2f' % (cluster, mean_metallicity, std_metallicity), verticalalignment='bottom', size=fontsize)


# Add 
ax = fig.add_subplot(xplots, yplots, i + 2)
aquarius_colors = ['c', 'g', 'm', 'r', 'b']
ax.errorbar(O_Fe, Na_Fe, xerr=e_O_Fe, yerr=e_Na_Fe, fmt=None, s=30, elinewidth=1.2, ecolor='k', zorder=-32)
ax.scatter(O_Fe, Na_Fe, edgecolor='k', facecolor=aquarius_colors, s=30, linewidths=1.2, zorder=10)


ax.set_xlabel('[O/Fe]')
ax.set_ylabel('[Na/Fe]')

ax.set_xlim(xlims)
ax.set_xticks([-0.5, 0.0, 0.5])

ax.set_ylim(ylims)
ax.set_yticks([-0.5, 0.0, 0.5, 1.0])

if i + 1 % yplots:
	ax.set_ylabel('')
	ax.set_yticklabels([''] * len(ax.get_yticklabels()))
	
ax.text(xlims[0] + 0.05 * (xlims[1] - xlims[0]), ylims[0] + 0.03 * (ylims[1] - ylims[0]), 'Aquarius Group\n[Fe/H] = %1.2f$\pm$%1.2f' % (-1.20, 0.33, ), verticalalignment='bottom', size=fontsize)


plt.savefig('aquarius-o-na-cluster.pdf')




# FUck it - do the O-Na normal one here too:

# Get data from Reddy et al 2005
reddy_data = np.loadtxt('Reddy_thick_disk_2005.data', delimiter='|', usecols=(1, 18, 20, 21), dtype=str)

Reddy_ThickDisk_O_Fe = []
Reddy_ThickDisk_Na_Fe = []
Reddy_ThinDisk_O_Fe = []
Reddy_ThinDisk_Na_Fe = []
Reddy_ThickThinDisk_O_Fe = []
Reddy_ThickThinDisk_Na_Fe = []
Reddy_Halo_O_Fe = []
Reddy_Halo_Na_Fe = []
Reddy_ThickDisk_Halo_O_Fe = []
Reddy_ThickDisk_Halo_Na_Fe = []

for classification, fe_h, O_h, Na_h in reddy_data:
	
	classification, fe_h, O_h, Na_h = [var.strip() for var in (classification, fe_h, O_h, Na_h)]


	Na_fe, O_fe = np.nan, np.nan

	if len(fe_h) == 0:
		fe_h, O_h, Na_h = [np.nan] * 3

	else:

		if len(O_h) == 0:
			O_h = np.nan

		else:
			O_fe = float(O_h) #- float(fe_h)


		if len(Na_h) == 0:
			Na_h = np.nan

		else:
			Na_fe = float(Na_h) #- float(fe_h)


	if classification == 'Thick disc':
		Reddy_ThickDisk_Na_Fe.append(Na_fe)
		Reddy_ThickDisk_O_Fe.append(O_fe)

	elif classification == 'Thin disc':
		Reddy_ThinDisk_Na_Fe.append(Na_fe)
		Reddy_ThinDisk_O_Fe.append(O_fe)

	elif classification == 'Halo':
		Reddy_Halo_Na_Fe.append(Na_fe)
		Reddy_Halo_O_Fe.append(O_fe)

	elif classification == 'Thin/thick disc':
		Reddy_ThickThinDisk_Na_Fe.append(Na_fe)
		Reddy_ThickThinDisk_O_Fe.append(O_fe)

	elif classification == 'Thick disc/Halo':
		Reddy_ThickDisk_Halo_Na_Fe.append(Na_fe)
		Reddy_ThickDisk_Halo_O_Fe.append(O_fe)


fig = plt.figure()
fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.10, right=0.95, bottom=0.10, top=0.95)

ax = fig.add_subplot(111)




cluster_idx = np.where(np.array(Carretta_Cluster) == 6171)[0]
cluster_idx = np.arange(len(Carretta_Na_limit))

idx = np.where(np.array(Carretta_Na_limit)[cluster_idx] == 0)[0]


ax.scatter(np.array(Carretta_O_Fe)[cluster_idx][idx], np.array(Carretta_Na_Fe)[cluster_idx][idx], marker='o', s=30, edgecolor='#cccccc', facecolor='none', zorder=10)

idx = np.where(np.array(Carretta_Na_limit)[cluster_idx] != 0)[0]

ax.errorbar(np.array(Carretta_O_Fe)[cluster_idx][idx], np.array(Carretta_Na_Fe)[cluster_idx][idx],fmt=None, ecolor='#cccccc', solid_capstyle='round', xerr=[[0.05] * len(idx), [0] * len(idx)], xlolims=True, zorder=-32)
ax.scatter(np.array(Carretta_O_Fe)[cluster_idx][idx], np.array(Carretta_Na_Fe)[cluster_idx][idx], s=30, marker='o', edgecolor='#cccccc', facecolor='w', zorder=5)




ax.scatter(Reddy_Halo_O_Fe, Reddy_Halo_Na_Fe, marker='x', s=30, edgecolor='k', label='Halo (Reddy et al. 2003)', zorder=50)
ax.scatter(Reddy_ThickDisk_Halo_O_Fe, Reddy_ThickDisk_Halo_Na_Fe, marker='x', s=30, edgecolor='k', zorder=50)

ax.scatter(Reddy_ThickDisk_O_Fe, Reddy_ThickDisk_Na_Fe, marker='+', s=30, edgecolor='k', label='Thick disc (Reddy et al. 2003)', zorder=50)
ax.scatter(Reddy_ThickThinDisk_O_Fe, Reddy_ThickThinDisk_Na_Fe, marker='+', s=30, edgecolor='k', zorder=50)

ax.scatter(Reddy_ThinDisk_O_Fe, Reddy_ThinDisk_Na_Fe, marker='+', s=30, edgecolor='#666666', label='Thin disc (Reddy et al. 2003)', zorder=50)



aquarius_colors = ['c', 'g', 'm', 'r', 'b']
ax.errorbar(O_Fe, Na_Fe, xerr=e_O_Fe, yerr=e_Na_Fe, fmt=None, elinewidth=1.2, ecolor='k', zorder=500)
ax.scatter(O_Fe, Na_Fe, edgecolor='k', facecolor=aquarius_colors, s=500, marker='*', linewidths=1.2, zorder=550)


ax.set_xlabel('[O/Fe]')
ax.set_ylabel('[Na/Fe]')

ax.set_xlim(xlims)
ax.set_xticks([-0.5, 0.0, 0.5])

ax.set_ylim(ylims)
ax.set_yticks([-0.5, 0.0, 0.5, 1.0])



plt.savefig('aquarius-o-na-halo.pdf')

