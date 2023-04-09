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
def read_and_fit(loop_paras):
    #! here is the fitting parameter setting
    ra_tmax = 8
    tau_cut = 1

    gamma, mass, mom, ll, b, z = loop_paras
    fit_id='{}{}_P{}_L{}_b{}_z{}_tmax{}_cut{}'.format(mass, gamma, mom, ll, b, z, ra_tmax, tau_cut)

    #* if the fit result already exists, skip it
    if os.path.exists('dump/gs_fit/{}_Q_chi_re_im'.format(fit_id)):
        return

    read_raw = Read_Raw('data_raw/')

    data_dic = {}
    temp_2pt = read_raw.read_2pt_bs(mass, mom)
    data_dic['2pt_re'] = np.real( temp_2pt )
    data_dic['2pt_im'] = np.imag( temp_2pt )


    for tseq in range(4, 9): # this tseq range comes from the raw data
        temp_ra_re, temp_ra_im = read_raw.read_ratio_bs(gamma, mass, mom, ll, b, z, tseq=tseq)

        data_dic['ra_re_tseq_{}'.format(tseq)] = temp_ra_re[:, 1:tseq]
        data_dic['ra_im_tseq_{}'.format(tseq)] = temp_ra_im[:, 1:tseq]

    gs_fit = Gs_Fit(two_state_fit(), fit_id)
    gs_fit.para_set(pt2_tmin=3, pt2_tmax=9, ra_tmin=4, ra_tmax=ra_tmax, tau_cut=tau_cut) #! here is the fitting parameter setting

    p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls = gs_fit.main(data_dic)

    Q_chi_re_im = {}
    Q_chi_re_im['Q'] = np.array(p_value_ls)
    Q_chi_re_im['chi2'] = np.array(chi2_ls)
    Q_chi_re_im['re'] = np.array(pdf_re_ls)
    Q_chi_re_im['im'] = np.array(pdf_im_ls)

    gv.dump(Q_chi_re_im, 'dump/gs_fit/{}_Q_chi_re_im'.format(fit_id))

    return



#* parallel processing
with mp.Pool() as pool:
    import warnings

    # Filter out RuntimeWarning messages
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    loop_paras_ls = [(gamma, mass, mom, ll, b, z) for ll in [8, 10] for gamma in ['t'] for mass in [220] for mom in [8] for b in [1,3,5] for z in range(13)]

    # pool.map(read_and_fit, loop_paras_ls)

    # results = list(tqdm(pool.imap(read_and_fit, loop_paras_ls), total=len(loop_paras_ls))) #* use tqdm to show the progress
#############################################################
#############################################################




# %%
#############################################################
'''
do the renormalization, make z dependence plots
'''
#############################################################

#* read the gs fit result
re_dic = {'L6':[], 'L8':[], 'L10':[]}
im_dic = {'L6':[], 'L8':[], 'L10':[]}

b = 1
mom = 8

for ll in [6, 8, 10]:
    for z in range(13):
        temp = gv.load('dump/gs_fit/220t_P{}_L{}_b{}_z{}_tmax8_cut1_Q_chi_re_im'.format(mom, ll, b, z))

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
