
# OMEGA CEN DATA
# Francois 1988
fr_feh, fr_e_feh, fr_bah, fr_e_bah = np.loadtxt('omegaCen_francois_1988.data', usecols=(1,2,3,4,), delimiter=',', unpack=True)
fr_bafe = fr_bah - fr_feh
fr_e_bafe = pow(pow(fr_e_feh, 2) + pow(fr_e_bah, 2), 0.5)

# Smith 2000
sm_logfe, sm_e_logfe, sm_logba, sm_e_logba = np.loadtxt('omegaCen_smith_2000.data', usecols=(1,2,3,4,), delimiter=',', unpack=True)
sm_feh = sm_logfe - 7.5
sm_bah = sm_logba - 2.18
sm_bafe = sm_bah - sm_feh
sm_e_feh = sm_e_logfe # Approx
sm_e_bafe = pow(pow(sm_e_feh, 2) + pow(sm_e_logba, 2), 0.5)

# Norris & Da Costa

# HALO DATA

# Reddy et al. 2003 + 2005
rd03_feh, rd03_bafe = np.loadtxt('Reddy_disk_2003.data', usecols=(4, 35,), delimiter='|', dtype=str, comments='#', unpack=True)

idx = np.where(np.array(map(len, [item.strip() for item in rd03_bafe])) > 0)[0]
rd03_feh = np.array(map(float, rd03_feh[idx]))
rd03_bafe = np.array(map(float, rd03_bafe[idx]))


rd05_feh, rd05_bah = np.loadtxt('Reddy_thick_disk_2005.data', usecols=(18, 36), delimiter='|', dtype=str, comments='#', unpack=True)
idx = np.where(np.array(map(len, [item.strip() for item in rd05_bah])) > 0)
rd05_feh = np.array(map(float, rd05_feh[idx]))
rd05_bah = np.array(map(float, rd05_bah[idx]))

rd05_bafe = rd05_bah - rd05_feh

# Fulbright 2000
with open('Fulbright_2000.data', 'r') as fp:
    content = fp.readlines()
    fb_feh = []
    fb_bafe = []

    for line in content:
        if line.startswith('#'): continue
        try:
            feh, bah = map(float, [line[8:13], line[87:92]])
        except:
            continue
        else:
            bafe = bah-feh

            fb_feh.append(feh)
            fb_bafe.append(bah)


# Marino et al 2011
with open('omegaCen_marino_2011.data', 'r') as fp:
    content = fp.readlines()
    mr11_feh = []
    mr11_feh2 = []
    mr11_bafe = []
    mr11_lafe = []

    for line in content:
        if line.startswith('#'): continue
        try:
            feh, bafe = map(float, [line[52:57], line[70:75]])
        except:
            continue
        else:
            
            mr11_feh.append(feh)
            mr11_bafe.append(bafe)


        try:
            feh, lafe = map(float, [line[52:57], line[76:81]])
        except:
            continue
        else:
            
            mr11_feh2.append(feh)
            mr11_lafe.append(lafe)




fig = plt.figure()
fig.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95)
ax = fig.add_subplot(111)

halo_color = 'k'
omegacen_color = '#bbbbbb'


aquarius_feh = np.array([-1.23, -1.13, -1.58, -0.62, -1.43])
aquarius_e_feh = np.array([0.10, 0.10, 0.10, 0.10, 0.10])

aquarius_ba_fe = np.array([0.62, -0.06, 0.00, 0.10, 0.03])
aquarius_e_ba_fe = np.array([0.15, 0.15, 0.15, 0.15, 0.15])


# Now with updated errors and abundances
aquarius_feh = np.array([-1.22, -1.13, -1.58, -0.63, -1.43])
aquarius_e_feh = np.array([0.13, 0.15, 0.16, 0.14, 0.12])

aquarius_ba_fe = np.array([0.62, -0.10, 0.00, 0.10, 0.03])
aquarius_e_ba_fe = np.array([0.16, 0.17, 0.20, 0.23, 0.13])

aquarius_colors = ['c', 'g', 'm', 'r', 'b']

# Other stars:

# Plot omega-Cen data

# Francois 1988
ax.errorbar(fr_feh, fr_bafe, xerr=fr_e_feh, yerr=fr_e_bafe, fmt=None, ecolor=omegacen_color, edgecolor=omegacen_color, zorder=-1)
ax.scatter(fr_feh, fr_bafe, marker='o', facecolor=omegacen_color, edgecolor=omegacen_color, zorder=1)

# Smith 2000
ax.errorbar(sm_feh, sm_bafe, xerr=sm_e_feh, yerr=sm_e_bafe, fmt=None, ecolor=omegacen_color, edgecolor=omegacen_color, zorder=-1)
ax.scatter(sm_feh, sm_bafe, marker='o', facecolor=omegacen_color, edgecolor=omegacen_color, zorder=1)

# Marino
ax.errorbar(mr11_feh, mr11_bafe, xerr=0.10, yerr=0.15, fmt=None, ecolor=omegacen_color, edgecolor=omegacen_color, zorder=-1)
ax.scatter(mr11_feh, mr11_bafe, marker='s', facecolor=omegacen_color, edgecolor=omegacen_color, zorder=1)


# Plot Halo data
ax.scatter(rd03_feh, rd03_bafe, marker='+', facecolor=halo_color, edgecolor=halo_color, s=20, zorder=2)
ax.scatter(rd05_feh, rd05_bah, marker='x', facecolor=halo_color, edgecolor=halo_color,s=20, zorder=2)
ax.scatter(fb_feh, fb_bafe, marker='o', facecolor=halo_color, edgecolor='none', s=20, zorder=2)

#ax.errorbar([c2225316[0]], [c2225316[1]], xerr=e_c2225316[0], yerr=e_c2225316[1], fmt=None, elinewidth=2, ecolor='k')
#ax.scatter([c2225316[0]], [c2225316[1]], marker='*', s=500, facecolor='c', edgecolor='k', zorder=10)
ax.errorbar(aquarius_feh, aquarius_ba_fe, xerr=aquarius_e_feh, yerr=aquarius_e_ba_fe, fmt=None, elinewidth=2, ecolor='k')
ax.scatter(aquarius_feh, aquarius_ba_fe, marker='*', s=500, facecolor=aquarius_colors, edgecolor='k', zorder=10)


# Lines?

ax.set_xlabel('[Fe/H]')
ax.set_ylabel('[Ba/Fe]')


ax.set_xlim(-2.5, 0.)
ax.set_ylim(-1., 1.5)

plt.savefig('omegaCen-c222531.pdf')
plt.savefig('omegaCen-c222531.eps')
plt.savefig('omegaCen-c222531.png')


fig = plt.figure()
fig.subplots_adjust(left=0.10, bottom=0.10, right=0.95, top=0.95)
ax = fig.add_subplot(111)



# Marino
ax.errorbar(mr11_feh2, mr11_lafe, xerr=0.10, yerr=0.15, fmt=None, ecolor=omegacen_color, edgecolor=omegacen_color, zorder=-1)
ax.scatter(mr11_feh2, mr11_lafe, marker='s', facecolor=omegacen_color, edgecolor=omegacen_color, zorder=1)

plt.savefig('omegaCen-lafe.pdf')


