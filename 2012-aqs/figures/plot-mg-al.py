import matplotlib.pyplot as plt
import numpy as np

# Get data from Reddy et al 2003
reddy_2003 = np.loadtxt('Reddy_disk_2003.data', delimiter='|', usecols=(18, 19), dtype=str)

Reddy_2003_Mg_Fe = []
Reddy_2003_Al_Fe = []

for mg, al in reddy_2003:
	mg = mg.strip()
	al = al.strip()

	if len(mg) == 0:
		mg = np.nan

	else:
		mg = float(mg)

	if len(al) == 0:
		al = np.nan

	else:
		al = float(al)

	Reddy_2003_Al_Fe.append(al)
	Reddy_2003_Mg_Fe.append(mg)


Fulbright_2000_Al_Fe = []
Fulbright_2000_Mg_Fe = []

with open('Fulbright_2000.data', 'r') as fp:
	data = fp.readlines()

	for line in data:
		if line.startswith('#'): continue

		mg = line[26:31].strip()
		al = line[32:37].strip()

		if len(mg) == 0:
			mg = np.nan

		else:
			mg = float(mg)

		if len(al) == 0:
			al = np.nan

		else:
			al = float(al)

		Fulbright_2000_Mg_Fe.append(mg)
		Fulbright_2000_Al_Fe.append(al)


# Get data from Reddy et al 2005
reddy_data = np.loadtxt('Reddy_thick_disk_2005.data', delimiter='|', usecols=(1, 18, 22, 23), dtype=str)

Reddy_ThickDisk_Mg_Fe = []
Reddy_ThickDisk_Al_Fe = []
Reddy_ThinDisk_Mg_Fe = []
Reddy_ThinDisk_Al_Fe = []
Reddy_ThickThinDisk_Mg_Fe = []
Reddy_ThickThinDisk_Al_Fe = []
Reddy_Halo_Mg_Fe = []
Reddy_Halo_Al_Fe = []
Reddy_ThickDisk_Halo_Mg_Fe = []
Reddy_ThickDisk_Halo_Al_Fe = []

for classification, fe_h, mg_h, al_h in reddy_data:
	
	classification, fe_h, mg_h, al_h = [var.strip() for var in (classification, fe_h, mg_h, al_h)]


	al_fe, mg_fe = np.nan, np.nan

	if len(fe_h) == 0:
		fe_h, mg_h, al_h = [np.nan] * 3

	else:

		if len(mg_h) == 0:
			mg_h = np.nan

		else:
			mg_fe = float(mg_h) #- float(fe_h)


		if len(al_h) == 0:
			al_h = np.nan

		else:
			al_fe = float(al_h) #- float(fe_h)


	if classification == 'Thick disc':
		Reddy_ThickDisk_Al_Fe.append(al_fe)
		Reddy_ThickDisk_Mg_Fe.append(mg_fe)

	elif classification == 'Thin disc':
		Reddy_ThinDisk_Al_Fe.append(al_fe)
		Reddy_ThinDisk_Mg_Fe.append(mg_fe)

	elif classification == 'Halo':
		Reddy_Halo_Al_Fe.append(al_fe)
		Reddy_Halo_Mg_Fe.append(mg_fe)

	elif classification == 'Thin/thick disc':
		Reddy_ThickThinDisk_Al_Fe.append(al_fe)
		Reddy_ThickThinDisk_Mg_Fe.append(mg_fe)

	elif classification == 'Thick disc/halo':
		Reddy_ThickDisk_Halo_Al_Fe.append(al_fe)
		Reddy_ThickDisk_Halo_Mg_Fe.append(mg_fe)



Reddy_ThickDisk_Mg_Fe, Reddy_ThickDisk_Al_Fe, Reddy_ThinDisk_Mg_Fe, Reddy_ThinDisk_Al_Fe, Reddy_Halo_Mg_Fe, Reddy_Halo_Al_Fe, Reddy_ThickThinDisk_Mg_Fe, Reddy_ThickThinDisk_Al_Fe, Reddy_ThickDisk_Halo_Mg_Fe, Reddy_ThickDisk_Halo_Al_Fe = [np.array(item) for item in (Reddy_ThickDisk_Mg_Fe, Reddy_ThickDisk_Al_Fe, Reddy_ThinDisk_Mg_Fe, Reddy_ThinDisk_Al_Fe, Reddy_Halo_Mg_Fe, Reddy_Halo_Al_Fe, Reddy_ThickThinDisk_Mg_Fe, Reddy_ThickThinDisk_Al_Fe, Reddy_ThickDisk_Halo_Mg_Fe, Reddy_ThickDisk_Halo_Al_Fe)]


# Caretta et al. 2009

Carretta_Mg_Fe, Carretta_Al_Fe, Carretta_Cluster, Carretta_Fe_H = np.loadtxt('Carretta_2009a.data', delimiter='|', usecols=(-3, -2, 2, -7), dtype=str, unpack=True)

delete_rows = []
for i in xrange(len(Carretta_Cluster)):

	Mg_Fe = Carretta_Mg_Fe[i].strip()
	Al_Fe = Carretta_Al_Fe[i].strip()

	if len(Mg_Fe) * len(Al_Fe) == 0: delete_rows.append(i)

Carretta_Cluster = map(float, np.delete(Carretta_Cluster, delete_rows))
Carretta_Fe_H = map(float, np.delete(Carretta_Fe_H, delete_rows))
Carretta_Mg_Fe = map(float, np.delete(Carretta_Mg_Fe, delete_rows))
Carretta_Al_Fe = map(float, np.delete(Carretta_Al_Fe, delete_rows))



data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Mg_Fe = []
e_Mg_Fe = []

Al_Fe = []
e_Al_Fe = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Mg')
	_Mg_Fe, _e_Mg_Fe = data[idx, 1:]

	idx = list(data[:, 0]).index('Al')
	_Al_Fe, _e_Al_Fe = data[idx, 1:]

	Mg_Fe.append(_Mg_Fe)
	Al_Fe.append(_Al_Fe)
	e_Mg_Fe.append(_e_Mg_Fe)
	e_Al_Fe.append(_e_Al_Fe)


# Clean up Mg_Fe, Al_Fe, e_Mg_Fe, e_Al_Fe
Mg_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Mg_Fe])
Al_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Al_Fe])

e_Mg_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Mg_Fe])
e_Al_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Al_Fe])

Mg_Fe, Al_Fe, e_Mg_Fe, e_Al_Fe = map(np.array, [Mg_Fe, Al_Fe, e_Mg_Fe, e_Al_Fe])


A = np.vstack([Mg_Fe, np.ones(len(Mg_Fe))]).T
m, c = np.linalg.lstsq(A, Al_Fe)[0]

fig = plt.figure()
fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.10, right=0.95, bottom=0.10, top=0.95)


ax = fig.add_subplot(111)

# Draw reddy
ax.scatter(Reddy_ThickDisk_Mg_Fe, Reddy_ThickDisk_Al_Fe, marker='+', edgecolor='k', facecolor='none', label='Thick disc (Reddy+ 2003)')
ax.scatter(Reddy_ThinDisk_Mg_Fe, Reddy_ThinDisk_Al_Fe, marker='+', facecolor='none', edgecolor='#666666', label='Thin disc (Reddy+ 2003)')
ax.scatter(Reddy_Halo_Mg_Fe, Reddy_Halo_Al_Fe, marker='+', facecolor='none', edgecolor='#cccccc', label='Halo (Reddy+ 2003)')

#ax.scatter(Reddy_ThickThinDisk_Mg_Fe, Reddy_ThickThinDisk_Al_Fe, marker='o', facecolor='none', edgecolor='m')
#ax.scatter(Reddy_ThickDisk_Halo_Mg_Fe, Reddy_ThickDisk_Halo_Al_Fe, marker='o', facecolor='none', edgecolor='orange')

#ax.scatter(Reddy_2003_Mg_Fe, Reddy_2003_Al_Fe, marker='o', facecolor='none', edgecolor='#666666')
ax.scatter(Fulbright_2000_Mg_Fe, Fulbright_2000_Al_Fe, marker='o', facecolor='none', edgecolor='#666666', label='Halo & Disc (Fulbright 2000)')

unique_clusters = list(np.unique(Carretta_Cluster))
#ax.scatter(Carretta_Mg_Fe, Carretta_Al_Fe, marker='o', c=[unique_clusters.index(cluster) for cluster in Carretta_Cluster], cmap=matplotlib.cm.jet, vmin=0, vmax=len(unique_clusters))

ax.scatter(Mg_Fe, Al_Fe, marker='o', color='k', s=30)
ax.errorbar(Mg_Fe, Al_Fe, xerr=e_Mg_Fe, yerr=e_Al_Fe, fmt=None, ecolor='k')

# C2225316
ax.text(0.50 + 0.03, 0.68 + 0.05, 'C2225316-14437', color='k', fontsize=10)

x = np.array([np.min(Mg_Fe), np.max(Mg_Fe)])
#ax.plot(x, m * x + c, 'k-')

ax.set_xlabel('[Mg/Fe]')
ax.set_ylabel('[Al/Fe]')

#xlims = (-0.4, 0.9)
#ylims = (-0.4, 1.9)


ax.set_xlim(-0.4, 0.9)
ax.set_ylim(-0.4, 1.9)

ax.set_xticks([0.0, 0.5])
ax.set_yticks([0.0, 0.5, 1.0, 1.5])

ax.legend(loc=2)

plt.savefig('aquarius-mg-al.eps')
plt.savefig('aquarius-mg-al.pdf')

fig = plt.figure(figsize=(4, 10))
fig.subplots_adjust(hspace=0.0, wspace=0.0, left=0.15, right=0.90, bottom=0.05, top=0.95)


ax = fig.add_subplot(311)
ax.set_ylabel('[Al/Fe]')
ax.xaxis.set_visible(False)
ax.scatter(Reddy_ThickDisk_Mg_Fe, Reddy_ThickDisk_Al_Fe, marker='+', facecolor='none', linewidth=1, edgecolor='k', label='Thick disc')
ax.scatter(Reddy_ThinDisk_Mg_Fe, Reddy_ThinDisk_Al_Fe, marker='+', facecolor='none', linewidth=1, edgecolor='#cccccc', label='Thin disc')
ax.scatter(Reddy_Halo_Mg_Fe, Reddy_Halo_Al_Fe, marker='s', facecolor='none', edgecolor='k', label='Halo')
ax.scatter(Reddy_ThickThinDisk_Mg_Fe, Reddy_ThickThinDisk_Al_Fe, marker='+', facecolor='none', edgecolor='#666666')
ax.scatter(Reddy_ThickDisk_Halo_Mg_Fe, Reddy_ThickDisk_Halo_Al_Fe, marker='x', facecolor='none', edgecolor='k')



ax.legend(frameon=False)
l = ax.get_legend()

ax.scatter(Mg_Fe, Al_Fe, marker='*', color='b', s=80)
#ax.errorbar(Mg_Fe, Al_Fe, xerr=e_Mg_Fe, yerr=e_Al_Fe, fmt=None, ecolor='k')


ax = fig.add_subplot(312, sharex=ax, sharey=ax)
ax.set_ylabel('[Al/Fe]')
ax.xaxis.set_visible(False)
ax.scatter(Fulbright_2000_Mg_Fe, Fulbright_2000_Al_Fe, marker='x', edgecolor='k', label='Thick disc \& Halo')
ax.scatter(Mg_Fe, Al_Fe, marker='*', color='b', s=80)


ax = fig.add_subplot(313, sharex=ax, sharey=ax)
ax.scatter(Carretta_Mg_Fe, Carretta_Al_Fe, c=[unique_clusters.index(cluster) for cluster in Carretta_Cluster], cmap=matplotlib.cm.jet, vmin=0, vmax=len(unique_clusters), edgecolor='k')

ax.scatter(Mg_Fe, Al_Fe, marker='*', color='b', s=80)

ax.set_xlabel('[Mg/Fe]')
ax.set_ylabel('[Al/Fe]')

ax.set_xlim(-0.4, 0.7)
ax.set_ylim(-0.4, 1.4)
ax.set_yticks([-0.5, 0.0, 0.5, 1.0])



plt.savefig('aquarius-mg-al-2.pdf')


# Mg-Al in sub-plots
Carretta_Cluster = map(float, np.delete(Carretta_Cluster, delete_rows))
Carretta_Fe_H = map(float, np.delete(Carretta_Fe_H, delete_rows))
Carretta_Mg_Fe = map(float, np.delete(Carretta_Mg_Fe, delete_rows))
Carretta_Al_Fe = map(float, np.delete(Carretta_Al_Fe, delete_rows))

fig = plt.figure()
fig.subplots_adjust(hspace=0, wspace=0,right=0.97,top=0.97,bottom=0.08,left=0.08)
unique_clusters = np.unique(Carretta_Cluster)

yplots = 5
xplots = int(np.ceil(float(len(unique_clusters) + 1) / yplots))

xlims = (-0.4, 0.9)
ylims = (-0.4, 1.9)

for i, cluster in enumerate(unique_clusters):

	ax = fig.add_subplot(xplots, yplots, i + 1)

	index = np.where(Carretta_Cluster == cluster)[0]

	Cluster_Mg_Fe = np.array(Carretta_Mg_Fe)[index]
	Cluster_Al_Fe = np.array(Carretta_Al_Fe)[index]

	ax.scatter(Cluster_Mg_Fe, Cluster_Al_Fe, facecolor='k')

	A = np.vstack([Cluster_Mg_Fe, np.ones(len(Cluster_Mg_Fe))]).T
	m, c = np.linalg.lstsq(A, Cluster_Al_Fe)[0]

	# Polyfit
	p = np.poly1d(np.polyfit(Cluster_Mg_Fe, Cluster_Al_Fe, 1))

	x_range = np.array([np.min(Cluster_Mg_Fe), np.max(Cluster_Mg_Fe)])
	y_range = np.array([np.min(Cluster_Al_Fe), np.max(Cluster_Al_Fe)])

	ax.plot(x_range, m * x_range + c, '-', color='#666666')
	ax.plot(x_range, p(x_range), '-', color='b')


	ax.set_xlim(xlims)
	ax.set_xticks([0.0, 0.5])

	ax.set_ylim(ylims)
	ax.set_yticks([0.0, 0.5, 1.0, 1.5])

	ax.set_xlabel('[Mg/Fe]')
	ax.set_ylabel('[Al/Fe]')

	if i % (yplots): 
		ax.set_ylabel('')
		ax.set_yticklabels([''] * len(ax.get_yticklabels()))
	
	if len(unique_clusters) - i - 1 > xplots:

		ax.set_xlabel('')
		ax.set_xticklabels([''] * len(ax.get_xticklabels()))

		ax.xaxis.set_visible(False)
	
	ax.text(xlims[0] + 0.05 * (xlims[1] - xlims[0]), ylims[1] - 0.05 * (ylims[1] - ylims[0]), 'NGC %i' % (cluster, ), verticalalignment='top')


# Add 
ax = fig.add_subplot(xplots, yplots, i + 2)
ax.scatter(Mg_Fe, Al_Fe, facecolor='k')

A = np.vstack([Mg_Fe, np.ones(len(Mg_Fe))]).T
m, c = np.linalg.lstsq(A, Al_Fe)[0]

print 'Aquarius', m

x_range = np.array([np.min(Mg_Fe), np.max(Mg_Fe)])
y_range = np.array([np.min(Al_Fe), np.max(Al_Fe)])

ax.plot(x_range, m * x_range + c, '-', color='#666666')

ax.set_xlabel('[Mg/Fe]')
ax.set_ylabel('[Al/Fe]')

ax.set_xlim(xlims)
ax.set_xticks([0.0, 0.5])

ax.set_ylim(ylims)
ax.set_yticks([0.0, 0.5, 1.0, 1.5])

if i + 1 % yplots:
	ax.set_ylabel('')
	ax.set_yticklabels([''] * len(ax.get_yticklabels()))
	
ax.text(xlims[0] + 0.05 * (xlims[1] - xlims[0]), ylims[1] - 0.05 * (ylims[1] - ylims[0]), 'Aquarius', verticalalignment='top')


plt.savefig('aquarius-mg-al-cluster.pdf')
