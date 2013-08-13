import matplotlib.pyplot as plt
import numpy as np

# Get Na and O from files:
filenames = ['c2225316-14437_abundances.data', 'c2306265-085103_abundances.data', 'j221821-183424_abundances.data', 'j223504-152834_abundances.data', 'j223811-104126_abundances.data']

aquarius_colors = ['c', 'g', 'm', 'r', 'b']

# Make sure 'Fe I' is one of the species
species_to_plot = ['Sc II', 'V I', 'Cr I', 'Mn I', 'Fe I']
species_to_plot = ['Na I', 'Al I', 'K I', 'Fe II', 'Co I', 'Cu I', 'Zn I', 'Sr I', 'Y II', 'Zr I|Zr II',
'Ba II', 'La II', 'Ce II', 'Nd II', 'Eu II', 'Fe I']

species_to_plot = ['Zn I', 'Sr I', 'Y II', 'Zr I|Zr II',
'Ba II', 'La II', 'Ce II', 'Nd II', 'Eu II', 'Fe I']


abundances = np.zeros((len(filenames), len(species_to_plot)))
uncertainties = np.zeros((len(filenames), len(species_to_plot)))

# Fill with nan's as default
abundances[:, :] = np.nan
uncertainties[:, :] = np.nan

Fe_I_index = species_to_plot.index('Fe I')

for i, filename in enumerate(filenames):

    data = np.loadtxt(filename, dtype=str, usecols=(0, 1, -6, -4, -2, ))

    for j, species in enumerate(species_to_plot):

        full_species = [' '.join(row).replace('\\textsc{', '').rstrip('}') for row in data[:, :2]]

        # Is this a mixed species?
        mixed_species = species.split('|')
        for mixed_specie in mixed_species:
            idx = full_species.index(mixed_specie)

            if mixed_specie == 'Fe I':
                data_slice = [data[idx, 2], data[idx, 4]]
            else:
                data_slice = data[idx, 3:5]

            abundance, uncertainty = map(float, [value.replace('$', '').replace('\\', '').replace('nodata', 'nan') for value in data_slice])

            # Only finite values
            if np.isfinite(abundance):
                abundances[i, j] = abundance
                uncertainties[i, j] = uncertainty


fig = plt.figure(figsize=(8, 12))
fig.subplots_adjust(left=0.12, bottom=0.05, right=0.99, top=0.99, hspace=0.0)

xlim = np.array([-1.7, -0.5])

# limits and ticks
ylims = {
    'Na I': ( [-0.05, 0.4],     [0.0, 0.2, 0.4] ),
    'Al I': ( [-0.05, 0.80],    [0.0, 0.4, 0.8] ),
    'Zn I': ( [-0.15, 0.45],    [0.0, 0.2, 0.4]),
    'Sr I': ([-1.00, 0.10],     [-0.8, -0.4, 0.0]),
    'K I': ([-0.05, 0.95], [0.0, 0.4, 0.8]),
    'Y II': ([-0.4, 0.95],    [-0.4, 0.0, 0.4, 0.8]),
    'Zr I|Zr II': ([-0.30, 0.9], [-0.3, 0.0, 0.3, 0.6]),
    'Ba II': ([-0.30, 0.90], [-0.3, 0.0, 0.3, 0.6]),
    'La II': ([-0.30, 0.90], [-0.3, 0.0, 0.3, 0.6]),
    'Nd II': ([-0.50, 0.9], [-0.3, 0.0, 0.3, 0.6]),
    'Eu II': ([-0.1, 0.6], [0.0, 0.3, 0.6])
}


for i, species in enumerate(species_to_plot):

    if species == 'Fe I': continue

    label_species = species.split('|')
    if len(label_species) > 1:
        label_species = label_species[0] + ',II'

    else:
        label_species = label_species[0]

    axes = fig.add_subplot(len(species_to_plot) - 1, 1, i + 1)

    x_data, y_data = abundances[:, Fe_I_index], abundances[:, i]
    x_error, y_error = uncertainties[:, Fe_I_index], uncertainties[:, i]


    axes.errorbar(x_data, y_data, xerr=x_error, yerr=y_error, fmt=None, ecolor='k', elinewidth=1.2, zorder=-32)
    axes.scatter(x_data, y_data, facecolor=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)

    #ax.errorbar(Fe_H, Mg_Fe, xerr=e_Fe_H, yerr=e_Mg_Fe, fmt=None, ecolor='k', elinewidth=1.2, zorder=-32)
    #ax.scatter(Fe_H, Mg_Fe, facecolor=aquarius_colors, marker='*', edgecolor='k', linewidths=1.2, s=500, zorder=10)


    idx = np.isfinite(y_data)
    finite_x_data, finite_y_data = x_data[idx], y_data[idx]

    A = np.vstack([finite_x_data, np.ones(len(finite_x_data))]).T
    m, c = np.linalg.lstsq(A, finite_y_data)[0]

    # Exclude C2225316
    c222_idx = np.where(finite_x_data != -1.23)

    finite_x_data2 = finite_x_data[c222_idx]
    finite_y_data2 = finite_y_data[c222_idx]


    A2 = np.vstack([finite_x_data2, np.ones(len(finite_x_data2))]).T
    m2, c2 = np.linalg.lstsq(A2, finite_y_data2)[0]

    print "%s = %1.2f x [Fe I/H] + %1.2f" % (label_species, m, c)

    x = np.array([np.min(finite_x_data), np.max(finite_x_data)])
    axes.plot(x, m * x + c, '-', c="#666666", lw=1.2)   # Trend line
    axes.plot(x, m2 * x + c2, '-', c='c', lw=1.2), # Trend line without C2225316

    axes.plot(xlim, [0, 0], 'k:')        # Zero line

    # Y-axis label and limits
    axes.set_ylabel('[%s/Fe]' % (label_species, ), fontsize=16)
    ylim_max = np.max(np.abs(axes.get_ylim())) * 1.05
    axes.set_ylim(-ylim_max, ylim_max)

    ytick_max = np.round(ylim_max, decimals=1)
    ytick_max -= 0.2
    if species == 'Fe II':
        ytick_max = 0.05    

        axes.set_ylim(-0.10, 0.10)

    axes.set_yticks([-ytick_max, 0, ytick_max])
    axes.set_yticklabels(['$-$%1.1f' % (ytick_max, ), 0, ' %1.1f' % ytick_max], fontsize=16)


    axes.set_xlim(xlim)
    if ylims.has_key(species):
        axes.set_ylim(ylims[species][0])
        axes.set_yticks(ylims[species][1])

        labels = []
        for item in ylims[species][1]:
            if item == 0.0:
                labels.append(' $\,$0.0')
            elif item < 0.0:
                labels.append('$-$%1.1f' % abs(item))
            else:
                labels.append(' $\,$%1.1f' % item)

        axes.set_yticklabels(labels)

    if len(species_to_plot) - 1 != i + 1:
        axes.set_xticklabels([''])

    else:
        axes.set_xticklabels(['$-$%1.1f' % (np.abs(num), ) for num in axes.get_xticks()], fontsize=16)
        axes.set_xlabel('[Fe I/H]', fontsize=16)


plt.savefig('plot-x-on-fe.pdf')
plt.savefig('plot-x-on-fe.eps')
plt.savefig('plot-x-on-fe.png')
