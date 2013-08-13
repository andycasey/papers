
from load_data import *

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 17}

matplotlib.rc('font', **font)

program_stars_scatter_args = {
    "facecolor": "w",
    "marker": "s",
    "edgecolor": "k",
    "lw": 2,
    "zorder": 10,
    "s": 50
}

program_stars_not_OSS_scatter_args = {
    "marker": "x",
    "zorder": 10,
    "edgecolor": "k",
    "s": 50 
}

standard_stars_scatter_args = {
    "marker": "o",
    "s": 50,
    "lw": 2,
    "facecolor": "#666666",
    "zorder": 10,
    "edgecolor": "k"
}


def plot_ngc2419(ax):


    color = '#8b8b8b'
    cohen_scatter_args = {
        "marker": "^",
        "edgecolor": color,
        "facecolor": "w",
        "s": 50,
        "zorder": 10,
    }

    mucciarelli_scatter_args = {
        "marker": "v",
        "edgecolor": color,
        "facecolor": "w",
        "s": 50,
        "zorder": 10,
    }

    mucciarelli_headers = 'ID Teff logg vt [Fe/H] [Fe/H]_err [Mg/Fe] [Mg/Fe]_err [K/Fe] [K/Fe]_err [Ti/Fe] [Ti/Fe]_err [Ca/Fe] [Ca/Fe]_err'

    mg_fe, mg_fe_err, k_fe, k_fe_err = range(4)
    mucciarelli_columns = ['[Mg/Fe]', '[Mg/Fe]_err', '[K/Fe]', '[K/Fe]_err']

    mucciarelli_column_indices = [mucciarelli_headers.split().index(item) for item in mucciarelli_columns]
    mucciarelli = np.loadtxt('ngc2419-mucciarelli.txt', usecols=mucciarelli_column_indices)

    ax.errorbar(mucciarelli[:, mg_fe], mucciarelli[:, k_fe], xerr=mucciarelli[:, mg_fe_err], yerr=mucciarelli[:, k_fe_err], fmt=None, lw=1, ecolor=color, zorder=-1)
    ax.scatter(mucciarelli[:, mg_fe], mucciarelli[:, k_fe], **mucciarelli_scatter_args)

    # Cohen data
    cohen = np.loadtxt('ngc2419-cohen.txt', usecols=(1,2,3,4, ))

    ax.errorbar(cohen[:, mg_fe], cohen[:, k_fe], xerr=cohen[:, mg_fe_err], yerr=cohen[:, k_fe_err], fmt=None, lw=1, ecolor=color, zorder=-1)
    ax.scatter(cohen[:, mg_fe], cohen[:, k_fe], **cohen_scatter_args)

    ax.set_ylabel('[K/Fe]')
    ax.set_xlabel('[Mg/Fe]')



def plot_A_B(element_a, element_b, ax, standard_stars, program_stars, program_stars_not_OSS, ylim=None, limscale=0.10):

    print [star[abundance]['%s/Fe' % element_a] for star in program_stars],
    print [star[abundance]['%s/Fe' % element_b] for star in program_stars],
        

    ax.scatter(
        [star[abundance]['%s/Fe' % element_a] for star in program_stars],
        [star[abundance]['%s/Fe' % element_b] for star in program_stars],
        **program_stars_scatter_args)

    ax.scatter(
        [star[abundance]['%s/Fe' % element_a] for star in program_stars_not_OSS],
        [star[abundance]['%s/Fe' % element_b] for star in program_stars_not_OSS],
        **program_stars_not_OSS_scatter_args)
    
    ax.scatter(
        [star[abundance]['%s/Fe' % element_a] for star in standard_stars],
        [star[abundance]['%s/Fe' % element_b] for star in standard_stars],
        **standard_stars_scatter_args
        )


    if ylim is None:
        ylim = ax.get_ylim()

    # Program stars

    # A limit = lower limit
    # B limit = upper limit
    # A limit + B limit = 

    ax.errorbar(
        [star[abundance]['%s/Fe' % element_a] for star in program_stars],
        [star[abundance]['%s/Fe' % element_b] for star in program_stars],
        xerr=[star[uncertainty]['%s/Fe' % element_a] for star in program_stars],
        yerr=[star[uncertainty]['%s/Fe' % element_b] for star in program_stars],
        fmt=None,
        zorder=-1,
        lw=2,
        ecolor='k'
        )



    # Standard stars

    ax.errorbar(
        [star[abundance]['%s/Fe' % element_a] for star in standard_stars],
        [star[abundance]['%s/Fe' % element_b] for star in standard_stars],
        xerr=[star[uncertainty]['%s/Fe' % element_a] for star in standard_stars],
        yerr=[star[uncertainty]['%s/Fe' % element_b] for star in standard_stars],
        fmt=None,
        zorder=-1,
        lw=2,
        ecolor='k',
        )



def plot_X_Y(element_a, element_b, ax, standard_stars, program_stars, program_stars_not_OSS, ylim=None, limscale=0.05):

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        **program_stars_scatter_args
        )

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars_not_OSS],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars_not_OSS],
        **program_stars_not_OSS_scatter_args
        )

    '''
    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        facecolor='b',
        zorder=10,
        s=50
        )

    ax.scatter(
        hx,
        hy,
        facecolor='none',
        edgecolor='k',
        zorder=-1,
        s=130,
        )
    '''

    if ylim is None:
        ylim = ax.get_ylim()

    # Program stars

    # A limit = lower limit
    # B limit = upper limit
    # A limit + B limit = 
    lolimits = []
    uplimits = []
    yup = []
    ydown = []

    for A_limit, B_limit, yerr in zip(
        [star[limit]['%s/H' % element_a] for star in program_stars],
        [star[limit]['%s/H' % element_b] for star in program_stars],
        [(np.array(star[uncertainty]['%s/H' % element_a]) + np.array(star[uncertainty]['%s/H' % element_b]))/2. for star in program_stars]
        ):

        if A_limit and not B_limit:
            # lower limit
            lolimits.append(True)
            uplimits.append(False)

            yup.append(0)
            ydown.append(limscale * np.diff(ylim))

        elif B_limit and not A_limit:
            # Upper limit
            uplimits.append(True)
            lolimits.append(False)

            yup.append(limscale * np.diff(ylim))
            ydown.append(0)


        elif not A_limit and not B_limit:
            uplimits.append(False)
            lolimits.append(False)

            yup.append(yerr)
            ydown.append(yerr)

        else:
            yup.append(0)
            ydown.append(limscale * np.diff(ylim))

            uplimits.append(False)
            lolimits.append(True)

    ax.errorbar(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        xerr=[star[uncertainty]['Fe I/H'] for star in program_stars],
        yerr=(ydown, yup),
        fmt=None,
        zorder=-1,
        lw=2,
        ecolor='k',
        lolims=lolimits,
        uplims=uplimits,
        )



    # Standard stars


    # A limit = lower limit
    # B limit = upper limit
    # A limit + B limit = ??
    lolimits = []
    uplimits = []
    yup = []
    ydown = []
    for A_limit, B_limit, yerr in zip(
        [star[limit]['%s/H' % element_a] for star in standard_stars],
        [star[limit]['%s/H' % element_b] for star in standard_stars],
        [(star[uncertainty]['%s/H' % element_a] + star[uncertainty]['%s/H' % element_b])/2 for star in standard_stars]):
        if A_limit and not B_limit:
            # lower limit
            lolimits.append(True)
            uplimits.append(False)

            yup.append(0)
            ydown.append(limscale * np.diff(ylim))


        elif B_limit and not A_limit:
            # Upper limit
            uplimits.append(True)
            lolimits.append(False)

            yup.append(limscale * np.diff(ylim))
            ydown.append(0)


        elif not A_limit and not B_limit:
            uplimits.append(False)
            lolimits.append(False)

            yup.append(yerr)
            ydown.append(yerr)

        else:
            yup.append(0)
            ydown.append(limscale * np.diff(ylim))

            uplimits.append(False)
            lolimits.append(True)

    ax.errorbar(
        [star[abundance]['Fe I/H'] for star in standard_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in standard_stars],
        xerr=[star[uncertainty]['Fe I/H'] for star in standard_stars],
        yerr=(ydown, yup),
        fmt=None,
        zorder=-1,
        ecolor='k',
        lw=2,
        lolims=lolimits,
        uplims=uplimits,
        )

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in standard_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in standard_stars],
        **standard_stars_scatter_args
        )

    

venn_data = get_venn_data('Venn_2004_Erratum.data')


# [ALPHA/FE]

# [Alpha/Fe]::[Mg/Fe]
xlim = (-4, 0.5)
ylim = (-0.5, 1.0)

xticks = []
yticks = [-0.4, 0., 0.4, 0.8]

#fig = plt.figure()
fig = plt.figure(figsize=(8.625, 12.462))

fig.subplots_adjust(hspace=0.05, right=0.95, top=0.95)

ax = fig.add_subplot(411)
plot_venn_data(venn_data, ax, '[Fe/H]', '[Mg/Fe]')
ax.plot(xlim, [0, 0], 'k:', zorder=-1)
ax.plot([0, 0], ylim, 'k:', zorder=-1)
ax.xaxis.set_visible(False)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)

ax.errorbar(
    [star[abundance]['Fe I/H'] for star in program_stars],
    [star[abundance]['Mg I/Fe'] for star in program_stars],
    xerr=[star[uncertainty]['Fe I/H']   for star in program_stars],
    yerr=[star[uncertainty]['Mg I/Fe']  for star in program_stars],
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Mg I/Fe'] for star in program_stars]
    )

ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['Mg I/Fe'] for star in program_stars],
        **program_stars_scatter_args
        )

ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars_not_OSS],
        [star[abundance]['Mg I/Fe'] for star in program_stars_not_OSS],
        **program_stars_not_OSS_scatter_args
        )


ax.errorbar(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    [star[abundance]['Mg I/Fe'] for star in standard_stars],
    xerr=[star[uncertainty]['Fe I/H']   for star in standard_stars],
    yerr=[star[uncertainty]['Mg I/Fe']  for star in standard_stars],
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Mg I/Fe'] for star in standard_stars]
    )

ax.scatter(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    [star[abundance]['Mg I/Fe'] for star in standard_stars],
    **standard_stars_scatter_args
    )




# [Alpha/Fe]::[Ca/Fe]
ax = fig.add_subplot(412)
plot_venn_data(venn_data, ax, '[Fe/H]', '[Ca/Fe]')
ax.plot(xlim, [0, 0], 'k:', zorder=-1)
ax.plot([0, 0], ylim, 'k:', zorder=-1)
ax.xaxis.set_visible(False)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)

ax.errorbar(
    [star[abundance]['Fe I/H'] for star in program_stars],
    [star[abundance]['Ca I/Fe'] for star in program_stars],
    xerr=[star[uncertainty]['Fe I/H']   for star in program_stars],
    yerr=[star[uncertainty]['Ca I/Fe']  for star in program_stars],
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Ca I/Fe'] for star in program_stars]
    )


ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['Ca I/Fe'] for star in program_stars],
        **program_stars_scatter_args
        )

ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars_not_OSS],
        [star[abundance]['Ca I/Fe'] for star in program_stars_not_OSS],
        **program_stars_not_OSS_scatter_args
        )


ax.errorbar(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    [star[abundance]['Ca I/Fe'] for star in standard_stars],
    xerr=[star[uncertainty]['Fe I/H']   for star in standard_stars],
    yerr=[star[uncertainty]['Ca I/Fe']  for star in standard_stars],
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Ca I/Fe'] for star in standard_stars]
    )

ax.scatter(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    [star[abundance]['Ca I/Fe'] for star in standard_stars],
    **standard_stars_scatter_args
    )



# [Alpha/Fe]::[Ti/Fe]
ax = fig.add_subplot(413)
plot_venn_data(venn_data, ax, '[Fe/H]', '[Ti/Fe]')
ax.plot(xlim, [0, 0], 'k:', zorder=-1)
ax.plot([0, 0], ylim, 'k:', zorder=-1)
ax.xaxis.set_visible(False)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)

ax.errorbar(
    [star[abundance]['Fe I/H'] for star in program_stars],
    [star[abundance]['Ti II/Fe'] for star in program_stars],
    xerr=[star[uncertainty]['Fe I/H']   for star in program_stars],
    yerr=[star[uncertainty]['Ti II/Fe']  for star in program_stars],
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Ti II/Fe'] for star in program_stars]
    )


ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['Ti II/Fe'] for star in program_stars],
        **program_stars_scatter_args
        )

ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars_not_OSS],
        [star[abundance]['Ti II/Fe'] for star in program_stars_not_OSS],
        **program_stars_not_OSS_scatter_args
        )


ax.errorbar(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    [star[abundance]['Ti II/Fe'] for star in standard_stars],
    xerr=[star[uncertainty]['Fe I/H']   for star in standard_stars],
    yerr=[star[uncertainty]['Ti II/Fe']  for star in standard_stars],
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Ti II/Fe'] for star in standard_stars]
    )

ax.scatter(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    [star[abundance]['Ti II/Fe'] for star in standard_stars],
    **standard_stars_scatter_args
    )



# [Alpha/Fe]::[alpha/Fe]
ax = fig.add_subplot(414)
y_data = (venn_data['[Ca/Fe]'] + venn_data['[Mg/Fe]'] + venn_data['[Ti/Fe]'])/3.
plot_venn_data(venn_data, ax, '[Fe/H]', y_data=y_data)
ax.plot(xlim, [0, 0], 'k:', zorder=-1)
ax.plot([0, 0], ylim, 'k:', zorder=-1)
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)




program_y_data = np.array([star[abundance]['Mg I/Fe'] for star in program_stars]) \
    + np.array([star[abundance]['Ca I/Fe'] for star in program_stars]) \
    + np.array([star[abundance]['Ti II/Fe'] for star in program_stars])
program_y_data /= 3.

program_y_err = np.array([star[uncertainty]['Mg I/Fe'] for star in program_stars]) \
    + np.array([star[uncertainty]['Ca I/Fe'] for star in program_stars]) \
    + np.array([star[uncertainty]['Ti II/Fe'] for star in program_stars])
program_y_err /= 3.




standard_y_data = np.array([star[abundance]['Mg I/Fe'] for star in standard_stars]) \
    + np.array([star[abundance]['Ca I/Fe'] for star in standard_stars]) \
    + np.array([star[abundance]['Ti II/Fe'] for star in standard_stars])
standard_y_data /= 3.

standard_y_err = np.array([star[uncertainty]['Mg I/Fe'] for star in standard_stars]) \
    + np.array([star[uncertainty]['Ca I/Fe'] for star in standard_stars]) \
    + np.array([star[uncertainty]['Ti II/Fe'] for star in standard_stars])
standard_y_err /= 3.


ax.errorbar(
    [star[abundance]['Fe I/H'] for star in program_stars],
    program_y_data,
    xerr=[star[uncertainty]['Fe I/H']   for star in program_stars],
    yerr=program_y_err,
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2
    )

ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        program_y_data,
        **program_stars_scatter_args
        )

program_y_data = np.array([star[abundance]['Mg I/Fe'] for star in program_stars_not_OSS]) \
    + np.array([star[abundance]['Ca I/Fe'] for star in program_stars_not_OSS]) \
    + np.array([star[abundance]['Ti II/Fe'] for star in program_stars_not_OSS])
program_y_data /= 3.

ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars_not_OSS],
        program_y_data,
        **program_stars_not_OSS_scatter_args
        )



ax.errorbar(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    standard_y_data,
    xerr=[star[uncertainty]['Fe I/H']   for star in standard_stars],
    yerr=standard_y_err,
    fmt=None,
    zorder=-1,
    ecolor='k',
    lw=2,
    lolims=[star[limit]['Ti II/Fe'] for star in standard_stars]
    )

ax.scatter(
    [star[abundance]['Fe I/H']  for star in standard_stars],
    standard_y_data,
    **standard_stars_scatter_args
    )


ax.set_ylabel('[(Mg+Ca+Ti)/3Fe]')
plt.savefig('alpha-Fe.pdf')
plt.savefig('alpha-Fe.png')



# [Ba/Y]
xlim = (-4, 0.5)
ylim = (-1, 1.5)
fig = plt.figure()
fig.subplots_adjust(hspace=0.05, right=0.95, top=0.95)

ax = fig.add_subplot(111)
plot_venn_data(venn_data, ax, '[Fe/H]', '[Ba/Fe]-[Y/Fe]')
ax.plot(xlim, [0, 0], 'k:', zorder=-10)
ax.plot([0, 0], ylim, 'k:', zorder=-10)

plot_X_Y('Ba II', 'Y II', ax, standard_stars, program_stars, program_stars_not_OSS, ylim=ylim)

ax.set_xlabel('[Fe/H]')
ax.set_ylabel('[Ba/Y]')
ax.set_xlim(xlim)
ax.set_ylim(ylim)

plt.savefig('Ba-Y.pdf')
plt.savefig('Ba-Y.png')


# [Na/Fe] vs [Ni/Fe]
xlim = (-1, 0.5)
ylim = (-0.5, 0.5)
fig = plt.figure()
fig.subplots_adjust(hspace=0.05, right=0.95, top=0.95)

ax = fig.add_subplot(111)
plot_venn_data(venn_data, ax, '[Na/Fe]', '[Ni/Fe]')
ax.plot(xlim, [0, 0], 'k:', zorder=-10)
ax.plot([0, 0], ylim, 'k:', zorder=-10)

plot_A_B('Na I', 'Ni I', ax, standard_stars, program_stars, program_stars_not_OSS, ylim=ylim)

ax.set_xlabel('[Na/Fe]')
ax.set_ylabel('[Ni/Fe]')
ax.set_xlim(xlim)
ax.set_ylim(ylim)
plt.savefig('Na-Ni.pdf')
plt.savefig('Na-Ni.png')


# [Mg/K]



fig = plt.figure()
fig.subplots_adjust(hspace=0.05, right=0.95, top=0.95)

ax = fig.add_subplot(111)
plot_ngc2419(ax)
plot_A_B('Mg I', 'K I', ax, standard_stars, program_stars, program_stars_not_OSS, ylim=ylim)
xlim = ax.get_xlim()
ylim = ax.get_ylim()

ax.plot(xlim, [0, 0], 'k:', zorder=-10)
ax.plot([0, 0], ylim, 'k:', zorder=-10)
ax.set_xlim(xlim)
ax.set_ylim(ylim)

plt.savefig('Mg-K.pdf')
plt.savefig('Mg-K.png')



