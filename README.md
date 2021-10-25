# NESS_02_OverviewPaper
A repository containing code used to produce tables and plots in Scicluna et al (2021)

To generate Tables 1, 3, and 5 from the paper, use `python3 makeCDStables.py` from the Unix 
command line or `from makeCDStables import *` at the python prompt.

To generate Figure and Table 4 from the paper, type the following at the python prompt:
`from NESS_DPR_density import *`
`NESS_DPR_density()`
If you encounter errors regarding TeX, try `NESS_DPR_density(use_tex = False)`.
