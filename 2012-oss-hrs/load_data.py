import os
import numpy as np


def get_venn_data(filename):

    with open(filename, 'r') as fp:
        contents = fp.readlines()


    # Skip first 43 lines
    del contents[:43]

    # Skip the name (first 10) and reference (135+) in each line
    contents = [map(float, line[11:134].split()) for line in contents]
    
    keys = (
            'U', 'V', 'W', 'Tn', 'Tk', 'Ha', 'Gind',
            '[Fe/H]', '[Mg/Fe]', '[Ca/Fe]', '[Ti/Fe]',
            '[a/Fe]', '[Na/Fe]', '[Ni/Fe]', '[Y/Fe]',
            '[Ba/Fe]', '[La/Fe]', '[Eu/Fe]'
            )

    data = {}
    for i, key in enumerate(keys):
        data[key] = np.array([line[i] for line in contents])

    # Assign star to each component
    assigned_component = []
    possible_components = ('Thin', 'Thick', 'Halo', 'dSph')

    probabilities = np.array([data['Tn'], data['Tk'], data['Ha']]).T
    for i, probability in enumerate(probabilities):

        if np.max(probability) > 0:
            max_prob_index = int(np.argmax(probability))

        else:
            max_prob_index = 3

        # Overwrite components
        if max_prob_index != 3 and data["V"][i] < -200:
            # Retro
            assigned_component.append("Retro")

        elif max_prob_index != 3 and (data["U"][i]**2 + data["W"][i]**2)**0.5 > 350:
            # High velocity halo
            assigned_component.append("High")

        else:
            assigned_component.append(possible_components[max_prob_index])

    data['component'] = assigned_component

    return data




def plot_venn_data(data, ax, x_column, y_column=None, y_data=None):

    if y_column is None and y_data is None: raise TypeError

    scatter_args = {
        "Thin": {
            "marker": "o",
            "facecolor": "#e62d24",
            "edgecolor": "#e62d24",
            "alpha": 0.4,
        },
        "Thick": {
            "marker": "o",
            "facecolor": "#6fc06a",
            "edgecolor": "#6fc06a",
            "alpha": 0.4,
        },
        "Halo": {
            "marker": "o",
            "edgecolor": "#4cc6f8",
            "facecolor": "#4cc6f8",
            "alpha": 0.4,
        },
        "Retro": {
            "marker": "o",
            "facecolor": "k",
            "edgecolor": "k",
            "alpha": 0.4,
        },
        "High": {
            "marker": "o",
            "facecolor": "#137dca",
            "edgecolor": "#137dca",
            "alpha": 0.4,
        },
        "dSph": {
            "marker": "s",
            "facecolor": "none",
            "edgecolor": "k",
            "s": 50,
            "alpha": 0.4
        }
    }

    if "-" in x_column:
        x_element_a, x_element_b = [word.strip() for word in x_column.split("-")]
        x_data = data[x_element_a] - data[x_element_b]

    else:
        x_data = data[x_column]
        ax.set_xlabel(x_column)

    if y_column != None:
        if "-" in y_column:
            y_element_a, y_element_b = [word.strip() for word in y_column.split("-")]
            y_data = data[y_element_a] - data[y_element_b]

        else:
            y_data = data[y_column]
            ax.set_ylabel(y_column)


    # Separate by components
    possible_components = list(set(data["component"]))

    for component in possible_components:
        indices = np.where([row_component == component for row_component in data["component"]])[0]

        scatter_arguments = scatter_args[component]

        ax.scatter(x_data[indices], y_data[indices], **scatter_arguments)




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


def plot_X_Y(element_a, element_b, ax, standard_stars, program_stars, ylim=None, limscale=0.10):

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        facecolor='b',
        zorder=10,
        s=50
        )

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in program_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in program_stars],
        facecolor='b',
        zorder=10,
        s=50
        )

    '''
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
        [pow(pow(star[uncertainty]['%s/H' % element_a], 2) + pow(star[uncertainty]['%s/H' % element_b], 2), 0.5) for star in program_stars]
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
        lolims=lolimits,
        uplims=uplimits,
        )

    ax.scatter(
        [star[abundance]['Fe I/H'] for star in standard_stars],
        [star[abundance]['%s/H' % element_a] - star[abundance]['%s/H' % element_b] for star in standard_stars],
        facecolor='r',
        zorder=10,
        s=50
        )

    ax.axhline(0, xmin=-10, xmax=10, c='k', linestyle=':')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('[%s/%s]' % (element_a.split()[0], element_b.split()[0], ))


hd76932, hd122563, hd44007, hd41667, hd136316, hd141531, hd142948, oss3, oss6, oss8, oss14, oss18 = \
    [load_abundances(filename) for filename in ["%s.abundances" % (star, ) for star in "hd76932, hd122563, hd44007, hd41667, hd136316, hd141531, hd142948, oss-1, oss-2, oss-3, oss-4, oss-5".split(", ")]]

# Useful
abundance, uncertainty, limit = (0, 1, 2)

program_stars = (oss3, oss6, oss8, oss14, oss18, )

program_stars_not_OSS = (oss3, oss18, )
program_stars_that_are_OSS = (oss6, oss8, oss14, )
standard_stars = (hd44007, hd41667, hd76932, hd122563, hd136316, hd141531, hd142948, )



