def NESS_DPR_density(infile = 'NESS_Table1_scaled.vot'):
    t = Table.read(infile, format = 'votable')
    #t = t[t['Ext'] == 1].copy()
    #
    tiers = ['very low', 'low', 'intermediate', 'high', 'extreme']
    ktiers = []
    for tier in tiers:
        ktiers.append(np.nonzero(t['sample'] == tier)[0])
    plt = setPlotParams()
    logdpr = [np.log10(t['GRAMS_DPR'][k]) for k in ktiers]
    _ = plt.boxplot(logdpr, meanline = True, showmeans = True, labels = tiers, whis = (2.5, 97.5))
    plt.ylabel(r'$\log{({\rm DPR}\ [{\rm M}_\odot\ {\rm yr}^{-1}])}$')
    plt.savefig('NESS_DPR_density.pdf', bbox_inches = 'tight')
    plt.close()
    #total DPR for each bin
    DPR_total = np.array([t['GRAMS_DPR'][k].sum() for k in ktiers])
    """Computing volumes for each distance bin"""
    radius_old = np.array([0.0, 0.0, 0.400, 0.800, 2.0]) #old values used for original source list
    radius = np.array([0.250, 0.300, 0.600, 1.200, 3.0]) #values used for extension list
    scale_height = 0.1 #scale height in kpc for Tier 4 sources from Libby's plots.
    #Volumes for the higher tiers are affected by extinction in the Galactic Plane (|b| < 1.5 deg).
    #   For these tiers, volume_total = volume_inner + volume_outer, where
    #   volume_inner = (Tiers 2 and 3) = volume of sphere with inner radius
    #                = (Tier 4) = volume of cylinder with inner radius and height = scale height
    #                               for Tier 4 from Libby's plots.
    #   volume_outer = volume of two cones of inner base radius `radius_old` and outer base radius `radius`
    #                   with vertex angles (90 - b) deg.
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
