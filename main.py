'''
This is the main function of this program. It will demonstrate the whole process of TMDPDF analysis.

The main steps are:
1. g.s. fit
2. renormalization
3. extrapolation in the coordinate space and FT
4. matching
5. large Pz and physical pion mass extrapolation
6. combine the gamma t and gamma z to get the final result

'''

# %%
import numpy as np
import gvar as gv
import multiprocessing as mp
import os
from tqdm import tqdm

from funcs import *
from read_raw_module import Read_Raw
from gs_fit_module import Gs_Fit
from prior_setting import two_state_fit


# %%
#############################################################
'''
do the g.s. fit to all sets, save the log and figures, output the chi2_ls, p_value_ls, pdf_re_ls, pdf_im_ls
'''
#############################################################

#! bs fits
if False:
    from gs_fit_module import read_and_fit_bs

    #* parallel processing
    with mp.Pool() as pool:
        import warnings

        # Filter out RuntimeWarning messages
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        loop_paras_ls = [(gamma, mass, mom, ll, b, z) for ll in [6] for gamma in ['t', 'z'] for mass in [220, 310] for mom in [8, 10, 12] for b in range(1, 6) for z in range(13)]

        # pool.map(read_and_fit, loop_paras_ls)

        results = list(tqdm(pool.imap(read_and_fit_bs, loop_paras_ls), total=len(loop_paras_ls))) #* use tqdm to show the progress

#! gvar fit
if True:
    from gs_fit_module_Lfit import read_and_fit_gvar
    read_and_fit_gvar()

#############################################################
#############################################################




# %%
#############################################################
'''
do the renormalization, make z dependence plots
'''
#############################################################
if False:
    #* read the gs fit result
    re_dic = {'L6':[], 'L8':[], 'L10':[]}
    im_dic = {'L6':[], 'L8':[], 'L10':[]}

    b = 1
    mom = 8

    for ll in [6, 8, 10]:
        for z in range(13):
            temp = gv.load('dump/gs_fit_bs/220t_P{}_L{}_b{}_z{}_tmax8_cut1_Q_chi_re_im'.format(mom, ll, b, z))

            re = gv.dataset.avg_data(temp['re'], bstrap=True)
            im = gv.dataset.avg_data(temp['im'], bstrap=True)

            re_dic['L{}'.format(ll)].append(re)
            im_dic['L{}'.format(ll)].append(im)


    y_ls = [ [v.mean for v in re_dic['L{}'.format(ll)]] for ll in [6, 8, 10] ]
    yerr_ls = [ [v.sdev for v in re_dic['L{}'.format(ll)]] for ll in [6, 8, 10] ]
    label_ls = ['L6', 'L8', 'L10']

    errorbar_ls_plot([np.arange(13) for i in range(3)], y_ls, yerr_ls, label_ls, 'z dependence mix L at b={}, mom=8, real'.format(b))


    y_ls = [ [v.mean for v in im_dic['L{}'.format(ll)]] for ll in [6, 8, 10] ]
    yerr_ls = [ [v.sdev for v in im_dic['L{}'.format(ll)]] for ll in [6, 8, 10] ]
    label_ls = ['L6', 'L8', 'L10']

    errorbar_ls_plot([np.arange(13) for i in range(3)], y_ls, yerr_ls, label_ls, 'z dependence mix L at b={}, mom=8, imag'.format(b))

#############################################################
#############################################################






# %%
