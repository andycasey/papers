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
			mg_fe = float(mg_h) - float(fe_h)


		if len(al_h) == 0:
			al_h = np.nan

		else:
			al_fe = float(al_h) - float(fe_h)


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

ax = fig.add_subplot(111)

# Draw reddy
ax.scatter(Reddy_ThickDisk_Mg_Fe, Reddy_ThickDisk_Al_Fe, marker='o', edgecolor='r', facecolor='none')
ax.scatter(Reddy_ThinDisk_Mg_Fe, Reddy_ThinDisk_Al_Fe, marker='o', facecolor='none', edgecolor='g')
ax.scatter(Reddy_Halo_Mg_Fe, Reddy_Halo_Al_Fe, marker='o', facecolor='none', edgecolor='b')

ax.scatter(Reddy_ThickThinDisk_Mg_Fe, Reddy_ThickThinDisk_Al_Fe, marker='o', facecolor='none', edgecolor='m')
ax.scatter(Reddy_ThickDisk_Halo_Mg_Fe, Reddy_ThickDisk_Halo_Al_Fe, marker='o', facecolor='none', edgecolor='orange')

#ax.scatter(Reddy_2003_Mg_Fe, Reddy_2003_Al_Fe, marker='o', facecolor='none', edgecolor='#666666')
ax.scatter(Fulbright_2000_Mg_Fe, Fulbright_2000_Al_Fe, marker='o', facecolor='none', edgecolor='k')

ax.scatter(Mg_Fe, Al_Fe, color='k')
ax.errorbar(Mg_Fe, Al_Fe, xerr=e_Mg_Fe, yerr=e_Al_Fe, fmt=None, ecolor='k')


x = np.array([np.min(Mg_Fe), np.max(Mg_Fe)])
ax.plot(x, m * x + c, 'k-')

ax.set_xlabel('[Mg/Fe]')
ax.set_ylabel('[Al/Fe]')

#ax.set_xlim(-0.5, 1.0)
#ax.set_ylim(-0.5, 1.5)

plt.savefig('aquarius-mg-al.eps')
plt.savefig('aquarius-mg-al.pdf')
