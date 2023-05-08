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
if False:
    from gs_fit_module import read_and_fit_gvar
    read_and_fit_gvar()

#############################################################
#############################################################




# %%
#############################################################
'''
do the renormalization, make z dependence plots
'''
#############################################################
if True:
    #* read the gs fit result
    all_after_gs_fit = gv.load('read_from_here/all_after_gs_fit_2pt_tmin2_b45_cut0.pkl')

    re_dic = {'b1':[], 'b2':[], 'b3':[], 'b4':[], 'b5':[]}
    im_dic = {'b1':[], 'b2':[], 'b3':[], 'b4':[], 'b5':[]}

    mass = 220
    gamma = 't'
    mom = 8
    ll = 6

    for b in range(1, 6):
        for z in range(13):
            fit_id = '{}{}_P{}_L{}_b{}_z{}'.format(mass, gamma, mom, ll, b, z)

            temp = all_after_gs_fit[fit_id]

            #! renormalization
            zO_fix = 1.05
            wl_sqrt = np.sqrt( readin_wloop(b, z) )

            re = temp['re'] / (zO_fix * wl_sqrt)
            im = temp['im'] / (zO_fix * wl_sqrt)

            re_dic['b{}'.format(b)].append(re)
            im_dic['b{}'.format(b)].append(im)


    y_ls = [ [v.mean for v in re_dic['b{}'.format(b)]] for b in range(1, 6) ]
    yerr_ls = [ [v.sdev for v in re_dic['b{}'.format(b)]] for b in range(1, 6) ]
    label_ls = ['b1', 'b2', 'b3', 'b4', 'b5']

    errorbar_ls_plot([np.arange(13) for i in range(5)], y_ls, yerr_ls, label_ls, 'z dependence mix b, {}{}_P{}_L{}, real'.format(mass, gamma, mom, ll), ylim=[-0.5, 1.8])


    y_ls = [ [v.mean for v in im_dic['b{}'.format(b)]] for b in range(1, 6) ]
    yerr_ls = [ [v.sdev for v in im_dic['b{}'.format(b)]] for b in range(1, 6) ]
    label_ls = ['b1', 'b2', 'b3', 'b4', 'b5']

    errorbar_ls_plot([np.arange(13) for i in range(5)], y_ls, yerr_ls, label_ls, 'z dependence mix b, {}{}_P{}_L{}, imag'.format(mass, gamma, mom, ll), ylim=[-0.5, 1.8])

#############################################################
#############################################################



# %%
