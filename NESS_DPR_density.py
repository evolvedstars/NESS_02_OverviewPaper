from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt

def texify(nums, decimals = 2):
    #Account for any numbers that might not be in scientific notation
    f = '%.' + str(decimals) + 'E'
    from decimal import Decimal
    numbers = np.array([f % Decimal(n) for n in nums])
    #For each value in /numbers/, return a TeX string of the form $mantissa \times 10^{exponent}$
    me = np.array([str(uu).split('E') if 'E' in str(uu) else [str(uu), '0'] for uu in numbers])
    mantissa = [np.round(float(mant), decimals = decimals) for mant in me[:, 0]]
    exponent = [int(exp) for exp in me[:, 1]]
    return np.array(['$' + str(m) + '\\times 10^{' + str(e) + '}$' if e != 0 else str(m) for m, e in zip(mantissa, exponent)])

def setPlotParams(use_tex = True):
    #the latex preamble setup here may result in errors for some.
    plt.figure(figsize = (8, 8))
    params = {'legend.fontsize': 'x-large',
              'axes.labelsize':20,
              'axes.titlesize':20,
              'xtick.labelsize':20,
              'ytick.labelsize':20}
    if use_tex:
        params['text.usetex'] = True
        params['text.latex.preamble'] = r'\usepackage{bm}'
        plt.rcParams.update(params)
    else:
        plt.rcParams.update(params)
    return plt

def NESS_DPR_density(infile = 'NESS_Table1_scaled_newer.vot', use_tex = True):
    """
        Reproduce Figure 9 and Table 4 from Scicluna et al. (2021)
        Required input: VOTable containing data for Table 1 from the paper
    """
    t = Table.read(infile, format = 'votable')
    tiers = ['very low', 'low', 'intermediate', 'high', 'extreme']
    ktiers = []
    for tier in tiers:
        ktiers.append(np.nonzero(t['sample'] == tier)[0])

    """ Figure 9 from the paper """
    plt = setPlotParams(use_tex = use_tex)
    logdpr = [np.log10(t['GRAMS_DPR'][k]) for k in ktiers]
    _ = plt.boxplot(logdpr, meanline = True, showmeans = True, labels = tiers, whis = (2.5, 97.5))
    try:
        plt.ylabel(r'$\log{({\rm DPR}\ [{\rm M}_\odot\ {\rm yr}^{-1}])}$')
    except:
        plt.ylabel('log(DPR [Msun/yr])')
    plt.savefig('NESS_DPR_density.pdf', bbox_inches = 'tight')
    plt.close()
    # total DPR for each bin
    DPR_total = np.array([t['GRAMS_DPR'][k].sum() for k in ktiers])

    """Table 4 from the paper """
    #Volumes for the higher tiers are affected by extinction in the Galactic Plane (|b| < 1.5 deg).
    #   For these tiers, volume_total = volume_inner + volume_outer, where
    #   volume_inner = (Tiers 2 and 3) = volume of sphere with inner radius
    #                = (Tier 4) = volume of cylinder with inner radius and height = scale height
    #                               for Tier 4 as estimated from, e.g., Figure 8 in the paper.
    #   volume_outer = volume of two cones of inner base radius `radius_old` and outer base radius `radius`
    #                   with vertex angles (90 - b) deg.
    radius_old = np.array([0.0, 0.0, 0.400, 0.800, 2.0]) #old values used for original source list
    radius = np.array([0.250, 0.300, 0.600, 1.200, 3.0]) #values used for extension list
    scale_height = 0.1 #scale height in kpc for Tier 4 sources as estimated from, e.g., Figure 8 in the paper.
    volume_inner = 4/3 * np.pi * radius_old**3
    volume_inner[-1] = np.pi * radius_old[-1]**2 * scale_height
    b = 1.5 #Galactic latitude that is excluded due to the Plane
    volume_outer = 2/3 * np.pi * (radius**3 - radius_old**3) * np.tan(np.pi * b / 180)
    volume = volume_inner + volume_outer
    volume[:2] = 4/3 * np.pi * radius[:2]**3 #spherical volume only applicable to Tiers 0 and 1
    #Now compute the densities
    DPR_total_str = texify(DPR_total, decimals = 1)
    DPR_total_sum_str = texify(np.array([DPR_total.sum()]), decimals = 1)
    surface_area = np.pi * radius**2
    DPR_disc_integ = DPR_total / surface_area
    DPR_disc_integ_str = texify(DPR_disc_integ, decimals = 1)
    DPR_disc_integ_sum_str = texify(np.array([DPR_disc_integ.sum()]), decimals = 1)
    DPR_density = DPR_total / volume
    DPR_density_str = texify(DPR_density, decimals = 1)
    DPR_density_sum_str = texify(np.array([DPR_density.sum()]), decimals = 1)
    num_str = [str(len(kk)) for kk in ktiers]
    num_sum_str = str(np.array([len(kk) for kk in ktiers]).sum())

    print(r"    \begin{tabular}{llccc}")
    print(r"    \hline Tier & No. & Total DPR & Disc-averaged DPR & Volume-averaged DPR \\")
    print(r"     & & M$_\odot$\,yr$^{-1}$& M$_\odot$\,yr$^{-1}$\,kpc$^{-2}$ & M$_\odot$\,yr$^{-1}$\,kpc$^{-3}$ \\")
    print(r"    \hline\hline")
    rows = ['    ' + str(i) + '  &' + num_str[i] + '&' + DPR_total_str[i] + '& ' + DPR_disc_integ_str[i] + '& ' + DPR_density_str[i] + \
            r'\\' for i in range(len(DPR_total_str))]
    row_totals = '    Total  &' + num_sum_str + '&' + DPR_total_sum_str[0] + '& ' + DPR_disc_integ_sum_str[0] + '& ' + DPR_density_sum_str[0] + r'\\'
    for row in rows:
        print(row)
    print('    \hline')
    print(row_totals)
    print('    \hline')
    print("    \end{tabular}")
