import matplotlib.pyplot as plt
import numpy as np

# Get Na and O from files:
data_files = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

Observed_Na_Fe = []
e_Observed_Na_Fe = []

Observed_Ni_Fe = []
e_Observed_Ni_Fe = []

observed_data = {}
for filename in data_files:

	data = np.loadtxt(filename, dtype=str, usecols=(0, -4, -2,))
	observed_data[filename] = data

	idx = list(data[:,0]).index('Na')
	_Observed_Na_Fe, _e_Observed_Na_Fe = data[idx, 1:]

	idx = list(data[:, 0]).index('Ni')
	_Observed_Ni_Fe, _e_Observed_Ni_Fe = data[idx, 1:]

	Observed_Na_Fe.append(_Observed_Na_Fe)
	Observed_Ni_Fe.append(_Observed_Ni_Fe)
	e_Observed_Na_Fe.append(_e_Observed_Na_Fe)
	e_Observed_Ni_Fe.append(_e_Observed_Ni_Fe)


# Clean up Observed_Na_Fe, Observed_Ni_Fe, e_Observed_Na_Fe, e_Observed_Ni_Fe
Observed_Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Observed_Na_Fe])
Observed_Ni_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in Observed_Ni_Fe])

e_Observed_Na_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Observed_Na_Fe])
e_Observed_Ni_Fe = map(float, [value.replace('$', '').replace('\nodata', 'nan') for value in e_Observed_Ni_Fe])

Observed_Na_Fe, Observed_Ni_Fe, e_Observed_Na_Fe, e_Observed_Ni_Fe = map(np.array, [Observed_Na_Fe, Observed_Ni_Fe, e_Observed_Na_Fe, e_Observed_Ni_Fe])


A = np.vstack([Observed_Na_Fe, np.ones(len(Observed_Na_Fe))]).T
m, c = np.linalg.lstsq(A, Observed_Ni_Fe)[0]


# Nissen & Schuster 2010

NS_Na_Fe, NS_Ni_Fe, kind = np.loadtxt('nissen_schuster_2010.data', delimiter='|', skiprows=9, dtype=str, usecols=(6, 12, -1), unpack=True)

NS_Na_Fe = np.array(map(float, NS_Na_Fe))
NS_Ni_Fe = np.array(map(float, NS_Ni_Fe))

# Thick disk stars
thick_disk = np.where(kind == 'TD')[0]

# High alpha stars
high_alpha = list(np.where(kind == 'high-alpha')[0])
high_alpha.extend(np.where(kind == '(high-alpha)')[0])

high_alpha = np.array(high_alpha)

low_alpha = list(np.where(kind == 'low-alpha')[0])
low_alpha.extend(np.where(kind == '(low-alpha)')[0])

low_alpha = np.array(low_alpha)


# Fulbright

Fulbright_Na_Fe = [] #21-25
Fulbright_Ni_Fe = [] #70-74

with open('Fulbright_2000.data', 'r') as fp:
	Fulbright_data = fp.readlines()

	for line in Fulbright_data:
		if line.startswith('#'): continue

		Na_Fe = line[20:25].strip()
		Ni_Fe = line[69:74].strip()

		if len(Na_Fe) * len(Ni_Fe) > 0:
			Fulbright_Na_Fe.append(float(Na_Fe))
			Fulbright_Ni_Fe.append(float(Ni_Fe))

Fulbright_Na_Fe = np.array(Fulbright_Na_Fe)
Fulbright_Ni_Fe = np.array(Fulbright_Ni_Fe)


# Reddy et al. 2003

# [Na/Fe] = 17
# [Ni/Fe] = 30

Reddy_2003_Na_Fe, Reddy_2003_Ni_Fe = np.loadtxt('Reddy_disk_2003.data', delimiter='|', usecols=(17, 30, ), unpack=True, dtype=str)

delete_rows = []
for i in xrange(len(Reddy_2003_Na_Fe)):

	Na_Fe = Reddy_2003_Na_Fe[i].strip()
	Ni_Fe = Reddy_2003_Ni_Fe[i].strip()

	if len(Na_Fe) * len(Ni_Fe) == 0: delete_rows.append(i)

Reddy_2003_Na_Fe = map(float, np.delete(Reddy_2003_Na_Fe, delete_rows))
Reddy_2003_Ni_Fe = map(float, np.delete(Reddy_2003_Ni_Fe, delete_rows))

# [Na/H] = 21
# [Ni/H] = 32


# Get data from Reddy et al 2005
reddy_data = np.loadtxt('Reddy_thick_disk_2005.data', delimiter='|', usecols=(1, 18, 21, 32), dtype=str)

Reddy_ThickDisk_Na_Fe = []
Reddy_ThickDisk_Ni_Fe = []
Reddy_ThinDisk_Na_Fe = []
Reddy_ThinDisk_Ni_Fe = []
Reddy_ThickThinDisk_Na_Fe = []
Reddy_ThickThinDisk_Ni_Fe = []
Reddy_Halo_Na_Fe = []
Reddy_Halo_Ni_Fe = []
Reddy_ThickDisk_Halo_Na_Fe = []
Reddy_ThickDisk_Halo_Ni_Fe = []

for classification, fe_h, mg_h, al_h in reddy_data:
	
	classification, fe_h, mg_h, al_h = [var.strip() for var in (classification, fe_h, mg_h, al_h)]


	al_fe, mg_fe = np.nan, np.nan

	if len(fe_h) == 0:
		fe_h, mg_h, al_h = [np.nan] * 3

	else:

		if len(mg_h) == 0:
			mg_h = np.nan

		else:
			mg_fe = float(mg_h)# - float(fe_h)


		if len(al_h) == 0:
			al_h = np.nan

		else:
			al_fe = float(al_h)# - float(fe_h)


	if classification == 'Thick disc':
		Reddy_ThickDisk_Ni_Fe.append(al_fe)
		Reddy_ThickDisk_Na_Fe.append(mg_fe)

	elif classification == 'Thin disc':
		Reddy_ThinDisk_Ni_Fe.append(al_fe)
		Reddy_ThinDisk_Na_Fe.append(mg_fe)

	elif classification == 'Halo':
		Reddy_Halo_Ni_Fe.append(al_fe)
		Reddy_Halo_Na_Fe.append(mg_fe)

	elif classification == 'Thin/thick disc':
		Reddy_ThickThinDisk_Ni_Fe.append(al_fe)
		Reddy_ThickThinDisk_Na_Fe.append(mg_fe)

	elif classification == 'Thick disc/halo':
		Reddy_ThickDisk_Halo_Ni_Fe.append(al_fe)
		Reddy_ThickDisk_Halo_Na_Fe.append(mg_fe)



Reddy_ThickDisk_Na_Fe, Reddy_ThickDisk_Ni_Fe, Reddy_ThinDisk_Na_Fe, Reddy_ThinDisk_Ni_Fe, Reddy_Halo_Na_Fe, Reddy_Halo_Ni_Fe, Reddy_ThickThinDisk_Na_Fe, Reddy_ThickThinDisk_Ni_Fe, Reddy_ThickDisk_Halo_Na_Fe, Reddy_ThickDisk_Halo_Ni_Fe = [np.array(item) for item in (Reddy_ThickDisk_Na_Fe, Reddy_ThickDisk_Ni_Fe, Reddy_ThinDisk_Na_Fe, Reddy_ThinDisk_Ni_Fe, Reddy_Halo_Na_Fe, Reddy_Halo_Ni_Fe, Reddy_ThickThinDisk_Na_Fe, Reddy_ThickThinDisk_Ni_Fe, Reddy_ThickDisk_Halo_Na_Fe, Reddy_ThickDisk_Halo_Ni_Fe)]


# Reddy et al 2003
# Reddy et al 2005
# Nissen & Schuster 2010
# Fulbright 2000

text_x = 0.45
text_y = -0.20

fig = plt.figure()

ax = fig.add_subplot(411)
ax.xaxis.set_visible(False)

ax.scatter(Reddy_2003_Na_Fe, Reddy_2003_Ni_Fe, marker='^', color='k', label='Thick disc')
ax.legend(loc=2, prop={'size': 9})
ax.text(text_x, text_y, 'Reddy et al. (2003)', color='k', verticalalignment='bottom', horizontalalignment='right')


ax = fig.add_subplot(412, sharex=ax, sharey=ax)
ax.xaxis.set_visible(False)
# square, diamond, triangle
ax.scatter(Reddy_ThickDisk_Na_Fe, Reddy_ThickDisk_Ni_Fe, marker='+', edgecolor='k', label='Thick disc')
#ax.scatter(Reddy_ThickThinDisk_Na_Fe, Reddy_ThickThinDisk_Ni_Fe, marker='^', color='#666666', label='Thick/thin disc')
ax.scatter(Reddy_ThinDisk_Na_Fe, Reddy_ThinDisk_Ni_Fe, marker='x', edgecolor='k', label='Thin disc')
#ax.scatter(Reddy_ThickDisk_Halo_Na_Fe, Reddy_ThickDisk_Halo_Ni_Fe, marker='d', color='k', label='Thick disc/halo')
ax.scatter(Reddy_Halo_Na_Fe, Reddy_Halo_Ni_Fe, marker='o', edgecolor='k', facecolor='none', label='Halo')
ax.legend(loc=2, prop={'size': 9})
ax.text(text_x, text_y, 'Reddy et al. (2003)', color='k', verticalalignment='bottom', horizontalalignment='right')


ax = fig.add_subplot(413, sharex=ax, sharey=ax)
ax.xaxis.set_visible(False)
ax.scatter(Fulbright_Na_Fe, Fulbright_Ni_Fe, marker='o', edgecolor='k', facecolor='none', label='Halo & disc')
ax.legend(loc=2, prop={'size': 9})
ax.text(text_x, text_y, 'Fulbright (2000)', color='k', verticalalignment='bottom', horizontalalignment='right')


ax = fig.add_subplot(414, sharex=ax, sharey=ax)
ax.scatter(NS_Na_Fe[thick_disk], NS_Ni_Fe[thick_disk], marker='o', edgecolor='k', facecolor='#666666', label='Thick disc')
ax.scatter(NS_Na_Fe[low_alpha], NS_Ni_Fe[low_alpha], marker='o', edgecolor='k', facecolor='#cccccc', label='Halo low-$\\alpha$')
ax.scatter(NS_Na_Fe[high_alpha], NS_Ni_Fe[high_alpha], marker='o', edgecolor='k', facecolor='w', label='Halo high-$\\alpha$')

ax.legend(loc=2, prop={'size': 9})
ax.text(text_x, text_y, 'Nissen & Schuster (2010)', color='k', verticalalignment='bottom', horizontalalignment='right')

ax.scatter(Observed_Na_Fe, Observed_Ni_Fe, marker='s', facecolor='b')

ax.set_xlim(-0.5, 0.5)
ax.set_ylim(-0.25, 0.25)

plt.savefig('na-ni.pdf')

