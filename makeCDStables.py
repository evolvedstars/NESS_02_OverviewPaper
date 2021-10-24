"""
    Generate machine-readable versions of Tables 1, 3, and 5 from the NESS overview paper
    (Scicluna et al. 2021).
    The script requires the cdspyreadme package, which is installed if missing.

    To execute from the Unix command line:
        python3 makeCDStables.py
    To execute from the Python command line:
        from makeCDStables import *
"""
import numpy as np
from astropy.table import Table, MaskedColumn
from astropy import units
import subprocess, os
try:
    import cdspyreadme
except:
    print("Package cdspyreadme not found, installing...")
    subprocess.run(['pip3', 'install', 'cdspyreadme'])
    import cdspyreadme

def preptables():
    """ Table 3 """
    t3 = Table.read('table3.csv', format = 'csv')
    t3.remove_column('col0')
    cols = np.array(t3.colnames)
    k = np.nonzero(['F' in c for c in cols])[0]
    newcols = np.array(['SCUBA2_' + c for c in cols[k]])
    t3.rename_columns(cols[k].tolist(), newcols.tolist())
    t3['IRASPSC'] = [tt['IRASPSC'].replace('IRAS', '') for tt in t3]
    desc = ['IRAS PSC identifier', \
            'JCMT/SCUBA2 450 um continuum flux', 'uncertainty in JCMT/SCUBA2 450 um continuum flux', \
            'JCMT/SCUBA2 450 um continuum flux quality flag (-1 if flux is 3\sigma upper limit)', \
            'JCMT/SCUBA2 850 um continuum flux', 'uncertainty in JCMT/SCUBA2 850 um continuum flux', \
            'JCMT/SCUBA2 850 um continuum flux quality flag (-1 if flux is 3\sigma upper limit)', \
            'spectral index (see Section 4.4 in paper)', 'lower bound of 68\% confidence interval for spectral index', \
            'upper bound of 68\% confidence interval for spectral index', \
            'quality flag for confidence interval (-99, 0, 1)']
    # u1 = units.def_unit('dimensionless', 1 * units.dimensionless_unscaled)
    u1 = units.dimensionless_unscaled
    u2 = units.Jy
    unit = [u1, u2, u2, u1, u2, u2, u1, u1, u1, u1, u1]
    unit = ['-', 'Jy', 'Jy', '-', 'Jy', 'Jy', '-', '-', '-', '-', '-']
    for i, c in enumerate(t3.columns):
        t3[c].description = desc[i]
        t3[c].unit = unit[i]

    """ Table 5 """
    t5 = Table.read('table5.csv', format = 'csv')
    t5.remove_column('SIMBAD ID')
    t5.rename_columns(['FWHM (arcseconds)', 'Extent (arcseconds)', 'Notes'], \
                      ['HARP_CO32_FWHM', 'HARP_CO32_Extent', 'HARP_CO32_Extent_Notes'])
    t5['IRASPSC'] = [tt['IRASPSC'].replace('IRAS ', '') for tt in t5]
    t5['HARP_CO32_Extent'] = ['-1' if e.strip() == 'N' else '-2' if e.strip() == '?' else \
                              e for e in t5['HARP_CO32_Extent']]
    t5['HARP_CO32_Extent_flag'] = [int(e) if '-' in e else 0 for e in t5['HARP_CO32_Extent']]
    t5['HARP_CO32_Extent'] = MaskedColumn(data = t5['HARP_CO32_Extent'].astype(float), \
                                          mask = t5['HARP_CO32_Extent_flag'] < 0)
    # reorder
    t5 = t5['IRASPSC', 'HARP_CO32_FWHM', 'HARP_CO32_Extent', 'HARP_CO32_Extent_flag', 'HARP_CO32_Extent_Notes'].copy()

    desc = ['IRAS PSC identifier', 'FWHM of Gaussian fit to radial profile of velocity-integrated intensity', \
            'Full extent of Gaussian fit to radial profile of velocity-integrated intensity', \
            '-1: extended emission is marginal, -2: extended emission not detected', \
            'Notes']
    u1 = units.dimensionless_unscaled
    u2 = units.arcsecond
    unit = [u1, u2, u2, u1, u1]
    unit = ['-', 'arcsec', 'arcsec', '-', '-']
    for i, c in enumerate(t5.columns):
        t5[c].description = desc[i]
        t5[c].unit = unit[i]

    return t3, t5

table1 = Table.read('NESS_Table1_scaled_newer.vot', format = 'votable')
table3, table5 = preptables()

tmkr = cdspyreadme.CDSTablesMaker()
tmkr.title = """The Nearby Evolved Stars Survey II: Constructing a volume-limited sample \
    and first results from the James Clerk Maxwell Telescope (Scicluna et al. 2021)"""
tmkr.authors = """P. Scicluna, F. Kemper, I. McDonald, S. Srinivasan, A. Trejo, \
                  S. H. J. Wallstr\\"om, J. G. A. Wouterloot, J. Cami, J. Greaves, \
                  Jinhua He, D. T. Hoai, Hyosun Kim, O. C. Jones, H. Shinnaga, \
                  C. J. R. Clark, T. Dharmawardena, W. Holland, H. Imai, J. Th. van Loon, \
                  K. M. Menten, R. Wesson, H. Chawner, S. Feng, S. Goldman, F.C. Liu, \
                  H. MacIsaac, J. Tang, S. Zeegers, K. Amada, V. Antoniou, A. Bemis, \
                  M. L. Boyer, S. Chapman, X. Chen, S.-H. Cho, L. Cui, F. Dell'Agli, \
                  P. Friberg, S. Fukaya, H. Gomez, Y. Gong, M. Hadjara, C. Haswell, \
                  N. Hirano, S. Hony, H. Izumiura, M. Jeste, X. Jiang, T. Kaminski, \
                  N. Keaveney, J. Kim, K. E. Kraemer, Y.-J. Kuan, E. Lagadec, C.F. Lee, \
                  D. Li, S.-Y. Liu, T. Liu, I. de Looze, F. Lykou, C. Maraston, \
                  J. P. Marshall, M. Matsuura, C. Min, M. Otsuka, M. Oyadomari, \
                  H. Parsons, N. A. Patel, E. Peeters, T. A. Pham, J. Qiu, \
                  S. Randall, G. Rau, M. P. Redman, A. M. S. Richards, S. Serjeant, \
                  C. Shi, G. C. Sloan, M. W. L. Smith, J. A. Toal\\'{a}, S. Uttenthaler, \
                  P. Ventura, B. Wang, I. Yamamura, T. Yang, Y. Yun, F. Zhang, \
                  Y. Zhang, G. Zhao, M. Zhu and A. A. Zijlstra
                """
tmkr.date = 2021
tmkr.abstract = """The Nearby Evolved Stars Survey (NESS) is a volume-complete sample of \
                   $\sim$850 Galactic evolved stars within 3\,kpc at (sub-)mm wavelengths, \
                   observed in the CO $J = $ (2--1) and (3--2) rotational lines, and the \
                   sub-mm continuum, using the James Clark Maxwell Telescope and Atacama \
                   Pathfinder Experiment. NESS consists of five tiers, based on distances \
                   and dust-production rate (DPR). We define a new metric for estimating the \
                   distances to evolved stars and compare its results to \emph{Gaia} EDR3. \
                   Replicating other studies, the most-evolved, highly enshrouded objects \
                   in the Galactic Plane dominate the dust returned by our sources, and we \
                   initially estimate a total DPR of $4.7\times 10^{-5}$ M$_\odot$ yr$^{-1}$ \
                   from our sample. Our sub-mm fluxes are systematically higher and spectral \
                   indices are typically shallower than dust models typically predict. The \
                   450/850 $\mu$m spectral indices are consistent with the blackbody \
                   Rayleigh--Jeans regime, suggesting a large fraction of evolved stars have \
                   unexpectedly large envelopes of cold dust."""
tmkr.putRef("https://ui.adsabs.harvard.edu/abs/2021MNRAS.NNN.NNNNS/abstract", "NESS overview paper")
tmkr.author = "P. Scicluna+"
tmkr.catalogue = ""
tmkr.keywords = "surveys – catalogues – stars: AGB and post-AGB – stars: mass-loss – stars: winds, outflows"

tmkr.addTable(table1, name='table1.mrt', description='machine-readable table')
tmkr.addTable(table3, name='table3.mrt', description='machine-readable table')
tmkr.addTable(table5, name='table5.mrt', description='machine-readable table')

tmkr.toMRT()
tmkr.makeReadMe()
