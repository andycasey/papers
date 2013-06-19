import os
import numpy as np

import matplotlib.pyplot as plt

def load_abundances(filename, minimum_uncertainty=0.05):
    """Loads the [X/H] and [X/Fe] abundance ratios and uncertainties
    for a single star from the filename provided."""

    if not os.path.exists(filename):
        raise IOError("file does not exist (%s)" % (filename, ))

    if not np.isfinite(minimum_uncertainty) or minimum_uncertainty < 0:
        raise ValueError("minimum uncertainty must be a positive finite value")

    abundances = {}
    uncertainties = {}
    limits = {}

    with open(filename, "r") as fp:
        content = fp.read().split("\n")

        for line in content:
            if line.startswith("#") or line.startswith("%") or len(line) == 0: continue

            line = line.replace("&", "").replace("_nodata", "nan").rstrip("\\").rstrip("&")
            if line.split()[0] == "C":
                element, transition = ("C", "(CH)", )
                n_lines, log_eps, sigma_eps, x_on_h, x_on_fe, uncertainty = line.split()[2:]
                limit = False

            else:
                line_split = line.split()
                if len(line_split) > 6:
                    element, transition, n_lines, log_eps, sigma_eps, x_on_h, x_on_fe, uncertainty = line_split
                    limit = False

                else:
                    element, transition, n_lines, log_eps, x_on_h, x_on_fe = line_split
                    uncertainty, sigma_eps = ('0.15', 'nan', )
                    limit = True

                    if not "<" in log_eps:

                        print filename, element
                        assert "<" in log_eps


            # Clean up transition
            transition = transition.strip('_textsc{}')
            species_repr = '%s %s' % (element, transition, )

            # Fix types
            n_lines = int(n_lines)
            log_eps, sigma_eps, x_on_h, x_on_fe, uncertainty = map(np.float, [item.replace("$", "").replace("<", "") for item in [log_eps, sigma_eps, x_on_h, x_on_fe, uncertainty]])

            # x_on_h
            abundances['%s/H' % (species_repr, )] = x_on_h
            abundances['%s/Fe' % (species_repr, )] = x_on_fe

            limits['%s/H' % (species_repr, )] = limit
            limits['%s/Fe' % (species_repr, )] = limit
            
            uncertainties['%s/H' % (species_repr, )] = np.max([uncertainty, minimum_uncertainty])
            uncertainties['%s/Fe' % (species_repr, )] = np.max([uncertainty, minimum_uncertainty])

    return (abundances, uncertainties, limits)


hd76932, hd122563, hd44007, hd41667, hd136316, hd141531, hd142948, oss3, oss6, oss8, oss14, oss18 = \
    [load_abundances(filename) for filename in ["%s.abundances" % (star, ) for star in "hd76932, hd122563, hd44007, hd41667, hd136316, hd141531, hd142948, oss-1, oss-2, oss-3, oss-4, oss-5".split(", ")]]

# Useful
abundance, uncertainty, limit = (0, 1, 2)

program_stars = (oss3, oss6, oss8, oss14, oss18, )
standard_stars = (hd44007, hd41667, hd76932, hd122563, hd136316, hd141531, hd142948, )


# [alpha elements/Fe] vs [Fe/H]
fig = plt.figure()

ylim = (-0.10, 0.75)
xlim = (-3.1, -0.6)
yticks = [0.0, 0.2, 0.4, 0.6]
elements = "Mg I, Ca I, Ti II".split(", ")

for i, element in enumerate(elements):
    if i > 0:
        ax = fig.add_subplot(len(elements) + 1, 1, i + 1, sharex=ax)

    else:
        ax = fig.add_subplot(len(elements) + 1, 1, i + 1)
        fig.subplots_adjust(hspace=0)


    ax.errorbar(
        [star[abundance]['Fe I/H']  for star in program_stars],
        [star[abundance]['%s/Fe' % element] for star in program_stars],
        xerr=[star[uncertainty]['Fe I/H']   for star in program_stars],
        yerr=[star[uncertainty]['%s/Fe' % element]  for star in program_stars],
        fmt=None,
        zorder=-1,
        ecolor='k',
        lolims=[star[limit]['%s/Fe'% element] for star in program_stars]
        )

    ax.scatter(
        [star[abundance]['Fe I/H']  for star in program_stars],
        [star[abundance]['%s/Fe' % element] for star in program_stars],
        facecolor='b',
        zorder=10,
        )

    ax.errorbar(
        [star[abundance]['Fe I/H']  for star in standard_stars],
        [star[abundance]['%s/Fe' % element] for star in standard_stars],
        xerr=[star[uncertainty]['Fe I/H']   for star in standard_stars],
        yerr=[star[uncertainty]['%s/Fe' % element]  for star in standard_stars],
        fmt=None,
        zorder=-1,
        ecolor='k',
        lolims=[star[limit]['%s/Fe'% element] for star in standard_stars]
        )

    ax.scatter(
        [star[abundance]['Fe I/H']  for star in standard_stars],
        [star[abundance]['%s/Fe' % element] for star in standard_stars],
        facecolor='r',
        zorder=10,
        )

    ax.axhline(0.0, xmin=-10, xmax=10, c='k', linestyle=':')

    ax.xaxis.set_visible(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yticks(yticks)
    ax.set_ylabel('[%s/Fe]' % element)


# [Alpha/Fe]
ax = fig.add_subplot(len(elements) + 1, 1, len(elements) + 1)

#fig = plt.figure()
#ax = fig.add_subplot(111)

prog_alpha = []
prog_fe = []
prog_alpha_err = []
prog_fe_err = []

for star in program_stars:
    prog_alpha.append((star[abundance]['Mg I/Fe'] + star[abundance]['Ca I/Fe'] + star[abundance]['Ti II/Fe'])/3.)
    prog_fe.append(star[abundance]['Fe I/H'])

    prog_alpha_err.append(pow(pow(star[uncertainty]['Mg I/Fe'], 2) + pow(star[uncertainty]['Ca I/Fe'], 2) + pow(star[uncertainty]['Ti II/Fe'], 2), 0.5))
    prog_fe_err.append(star[uncertainty]['Fe I/H'])

std_alpha = []
std_fe = []
std_alpha_err = []
std_fe_err = []

for star in standard_stars:
    std_alpha.append((star[abundance]['Mg I/Fe'] + star[abundance]['Ca I/Fe'] + star[abundance]['Ti II/Fe'])/3.)
    std_fe.append(star[abundance]['Fe I/H'])

    std_alpha_err.append(pow(pow(star[uncertainty]['Mg I/Fe'], 2) + pow(star[uncertainty]['Ca I/Fe'], 2) + pow(star[uncertainty]['Ti II/Fe'], 2), 0.5))
    std_fe_err.append(star[uncertainty]['Fe I/H'])

ax.errorbar(std_fe, std_alpha, xerr=std_fe_err, yerr=std_alpha_err, ecolor='k', fmt=None, zorder=-1)
ax.errorbar(prog_fe, prog_alpha, xerr=prog_fe_err, yerr=prog_alpha_err, ecolor='k', fmt=None, zorder=-1)

ax.scatter(prog_fe, prog_alpha, facecolor='b', label='program',zorder=10,)
ax.scatter(std_fe, std_alpha, facecolor='r', label='std',zorder=10,)

print "Standard [Fe/H]", std_fe
print "Standard alpha", std_alpha

print "Program [Fe/H]", prog_fe
print "Program alpha", prog_alpha
ax.set_xlim(xlim)
ax.set_ylim(ylim)
ax.set_yticks(yticks)
ax.axhline(0.0, xmin=-10, xmax=10, c='k', linestyle=':')

ax.set_xlabel('[Fe/H]')
ax.set_ylabel('[$\\alpha$/Fe]')

plt.savefig('alpha-fe.pdf')


# Fe-peak vs [Fe/H]
fig = plt.figure()

ylim = (-0.65, 0.65)
yticks = (-0.4, 0.0, 0.4)
xlim = (-3.1, -0.6)
elements = "V I, Cr I, Mn I, Co I, Ni I, Cu I, Zn I".split(", ")

"""
  V _textsc{I} &   2 &    2.34 &    0.16 & $-$1.59 & $-$0.01 &    0.11 \
  Cr _textsc{I} &  17 &    3.90 &    0.13 & $-$1.74 & $-$0.16 &    0.03 \
 Cr _textsc{II} &   3 &    4.27 &    0.04 & $-$1.37 &    0.21 &    0.02 \
  Mn _textsc{I} &   5 &    3.59 &    0.19 & $-$1.84 & $-$0.26 &    0.08 \
  Fe _textsc{I} &  61 &    5.92 &    0.15 & $-$1.58 &    0.00 &    0.02 \
 Fe _textsc{II} &  12 &    5.93 &    0.16 & $-$1.57 &    0.01 &    0.05 \
  Co _textsc{I} &   4 &    3.45 &    0.16 & $-$1.54 &    0.04 &    0.08 \
  Ni _textsc{I} &  21 &    4.59 &    0.17 & $-$1.63 & $-$0.05 &    0.04 \
  Cu _textsc{I} &   1 &    1.72 &    0.00 & $-$2.47 & $-$0.89 &    0.00 \
  Zn _textsc{I} &   1 &    2.94 &    0.00 & $-$1.62 & $-$0.04 &    0.00 \
 Sr _textsc{II} &   2 &    0.52 &    0.36 & $-$2.35 & $-$0.77 &    0.25 \
  Y _textsc{II} &   1 & $<$0.12 &         &$<-$2.09 &$<-$0.51 &         \
% Zr _textsc{II} &   0 & _nodata & _nodata & _nodata & _nodata & _nodata \
 Ba _textsc{II} &   2 &    0.14 &    0.22 & $-$2.04 & $-$0.46 &    0.16 \
 La _textsc{II} &   1 &   -0.07 &    0.00 & $-$1.17 &    0.41 &    0.00 \
 Nd _textsc{II} &   1 &    0.13 &    0.00 & $-$1.29 &    0.29 &    0.00 \
 Eu _textsc{II} &   1 &$<-$0.62 &         &$<-$1.14 & $<$0.44 &         \
"""
for i, element in enumerate(elements):
    if i > 0:
        ax = fig.add_subplot(len(elements), 1, i + 1)

    else:
        ax = fig.add_subplot(len(elements), 1, i + 1)
        fig.subplots_adjust(hspace=0)

    ax.scatter(
        [star[abundance]['Fe I/H']  for star in program_stars],
        [star[abundance]['%s/Fe' % element] for star in program_stars],
        facecolor='b',
        zorder=10,
        )

    ax.scatter(
        [star[abundance]['Fe I/H']  for star in standard_stars],
        [star[abundance]['%s/Fe' % element] for star in standard_stars],
        facecolor='r',
        zorder=10,
        )

    #ylim = ax.get_ylim()
    if element.strip() == "Cu I":
        ylim_thisax = (-1.1, 1.1)
        print "This one is Cu I"

    else:
        ylim_thisax = (ylim[0], ylim[1])

    limits = [star[limit]['%s/Fe'% element] for star in program_stars]
    ax.errorbar(
        [star[abundance]['Fe I/H']  for star in program_stars],
        [star[abundance]['%s/Fe' % element] for star in program_stars],
        xerr=[star[uncertainty]['Fe I/H']   for star in program_stars],
        yerr=[[star[uncertainty]['%s/Fe' % element] if not is_limit else 0.2 * np.diff(ylim_thisax) for star, is_limit in zip(program_stars, limits)], [star[uncertainty]['%s/Fe' % element] if not is_limit else 0 for star, is_limit in zip(program_stars, limits)]],
        fmt=None,
        zorder=-1,
        ecolor='k',
        lolims=limits
        )

    limits = [star[limit]['%s/Fe'% element] for star in standard_stars]
    ax.errorbar(
        [star[abundance]['Fe I/H']  for star in standard_stars],
        [star[abundance]['%s/Fe' % element] for star in standard_stars],
        xerr=[star[uncertainty]['Fe I/H']   for star in standard_stars],
        yerr=[[star[uncertainty]['%s/Fe' % element] if not is_limit else 0.2 * np.diff(ylim_thisax) for star, is_limit in zip(standard_stars, limits)], [star[uncertainty]['%s/Fe' % element] if not is_limit else 0 for star, is_limit in zip(standard_stars, limits)]],
        fmt=None,
        zorder=-1,
        ecolor='k',
        lolims=[star[limit]['%s/Fe'% element] for star in standard_stars]
        )
    
    ax.axhline(0.0, xmin=-10, xmax=10, c='k', linestyle=':')

    ax.xaxis.set_visible(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim_thisax)
    if element.strip() == "Cu I":
        ax.set_yticks([-0.8, -0.4, 0.0, 0.4, 0.8])
    else:
        ax.set_yticks(yticks)
    ax.set_ylabel('[%s/Fe]' % element)

# Final plot should have an x-axis
ax.xaxis.set_visible(True)
ax.set_xlabel('[Fe/H]')

plt.savefig('fe-peak.pdf')

# Ba/Y

def plot_X_Y(element_a, element_b, ax, standard_stars, program_stars):

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        facecolor='b',
        zorder=10
        )

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        facecolor='b',
        zorder=10
        )

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
        [pow(pow(star[uncertainty]['%s/H' % element_a], 2) + pow(star[uncertainty]['%s/H' % element_b], 2), 0.5) for star in program_stars]
        ):
        if A_limit and not B_limit:
            # lower limit
            lolimits.append(True)
            uplimits.append(False)

            yup.append(0)
            ydown.append(0.05 * np.diff(ylim))

        elif B_limit and not A_limit:
            # Upper limit
            uplimits.append(True)
            lolimits.append(False)

            yup.append(0.05 * np.diff(ylim))
            ydown.append(0)


        elif not A_limit and not B_limit:
            uplimits.append(False)
            lolimits.append(False)

            yup.append(yerr)
            ydown.append(yerr)

        else:
            yup.append(0)
            ydown.append(yerr)

            uplimits.append(False)
            lolimits.append(True)

    ax.errorbar(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        xerr=[star[uncertainty]['Fe I/H'] for star in program_stars],
        yerr=(ydown, yup),
        fmt=None,
        zorder=-1,
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
        [pow(pow(star[uncertainty]['%s/H' % element_a], 2) + pow(star[uncertainty]['%s/H' % element_b], 2), 0.5) for star in standard_stars]):
        if A_limit and not B_limit:
            # lower limit
            lolimits.append(True)
            uplimits.append(False)

            yup.append(0)
            ydown.append(0.05 * np.diff(ylim))


        elif B_limit and not A_limit:
            # Upper limit
            uplimits.append(True)
            lolimits.append(False)

            yup.append(0.05 * np.diff(ylim))
            ydown.append(0)


        elif not A_limit and not B_limit:
            uplimits.append(False)
            lolimits.append(False)

            yup.append(yerr)
            ydown.append(yerr)

        else:
            yup.append(0)
            ydown.append(yerr)

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
        lolims=lolimits,
        uplims=uplimits,
        )

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in standard_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in standard_stars],
        facecolor='r',
        zorder=10
        )

    ax.axhline(0, xmin=-10, xmax=10, c='k', linestyle=':')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('[%s/%s]' % (element_a.split()[0], element_b.split()[0], ))



# [Y/Fe]
fig = plt.figure()

fig.subplots_adjust(hspace=0)
ax = fig.add_subplot(511)
plot_X_Y('Y II', 'Fe I', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [Ba/Fe]
ax = fig.add_subplot(512)
plot_X_Y('Ba II', 'Fe I', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [La/Fe]
ax = fig.add_subplot(513)
plot_X_Y('La II', 'Fe I', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [Eu/Fe]
ax = fig.add_subplot(514)
plot_X_Y('Eu II', 'Fe I', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [Sr/Fe]
ax = fig.add_subplot(515)
plot_X_Y('Sr II', 'Fe I', ax, standard_stars, program_stars)
# Sr II not affected much by non-LTE -- of the order +/- 0.05 dex
plt.savefig('heavy-Fe.pdf')


# [Y/Eu]
fig = plt.figure()
fig.subplots_adjust(hspace=0)
ax = fig.add_subplot(411)
plot_X_Y('Y II', 'Eu II', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [Ba/Eu]
ax = fig.add_subplot(412)
plot_X_Y('Ba II', 'Eu II', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [La/Eu]
ax = fig.add_subplot(413)
plot_X_Y('La II', 'Eu II', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)

# [Ba/Y]
ax = fig.add_subplot(414)
plot_X_Y('Ba II', 'Y II', ax, standard_stars, program_stars)
ax.xaxis.set_visible(False)


# S-process
fig = plt.figure()
ax = fig.add_subplot(111)
plot_X_Y('Ba II', 'Sr II', ax, standard_stars, program_stars)
#elements = "Ba II - Sr II, Ba II, Sr II, La II, Eu II".split(", ")



# [Na, Ni]
fig = plt.figure()
ax = fig.add_subplot(111)



ax.scatter(
    [star[abundance]['Na I/Fe']  for star in standard_stars],
    [star[abundance]['Ni I/Fe'] for star in standard_stars],
    facecolor='r',
    zorder=10,
    )


ax.scatter(
    [star[abundance]['Na I/Fe']  for star in program_stars],
    [star[abundance]['Ni I/Fe'] for star in program_stars],
    facecolor='b',
    zorder=10,
    )

ax.set_xlabel('[Na/Fe]')
ax.set_ylabel('[Ni/Fe]')

