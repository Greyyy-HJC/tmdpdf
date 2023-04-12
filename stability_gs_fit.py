'''
This code is used to make stability plot for the ground state fits.
Pick some (b, z) points, vary the tmin of the fit, and plot the comparison plot.
'''

# %%
import os
import numpy as np
import lsqfit as lsf

from funcs import *
from read_raw_module import Read_Raw
from gs_fit_module import Gs_Fit
from prior_setting import two_state_fit



def each_fit(tmin, b, z):
    #! bs fits
    gamma = 'z'
    mass = 220
    mom = 8
    ll = 6


    read_raw = Read_Raw('data_raw/')

    data_dic = {}
    temp_2pt = read_raw.read_2pt_bs(mass, mom)
    data_dic['2pt_re'] = np.real( temp_2pt )
    data_dic['2pt_im'] = np.imag( temp_2pt )


    for tseq in range(4, 9):
        temp_ra_re, temp_ra_im = read_raw.read_ratio_bs(gamma, mass, mom, ll, b, z, tseq=tseq)

        data_dic['ra_re_tseq_{}'.format(tseq)] = temp_ra_re[:, 1:tseq]
        data_dic['ra_im_tseq_{}'.format(tseq)] = temp_ra_im[:, 1:tseq]


    gs_fit = Gs_Fit(two_state_fit(), fit_id='{}{}_P{}_L{}_b{}_z{}_tmax9_cut0'.format(mass, gamma, mom, ll, b, z))
    gs_fit.para_set(pt2_tmin=3, pt2_tmax=9, ra_tmin=tmin, ra_tmax=9, tau_cut=0)


    #! bs fit
    # p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls = gs_fit.main_bs(data_dic)

    # p_value = np.mean(p_value_ls)
    # chi2 = np.mean(chi2_ls)
    # re_gv = gv.dataset.avg_data(pdf_re_ls, bstrap=True)
    # im_gv = gv.dataset.avg_data(pdf_im_ls, bstrap=True)


    #! gvar fit
    data_dic_avg = gv.dataset.avg_data(data_dic, bstrap=True)
    fit_res = gs_fit.main_gvar(data_dic_avg)


    p_value = fit_res.Q
    chi2 = fit_res.chi2 / fit_res.dof
    re_gv = fit_res.p['pdf_re']
    im_gv = fit_res.p['pdf_im']


    return p_value, chi2, re_gv, im_gv



tmin_ls = np.arange(4, 7)
p_value_ls = []
chi2_ls = []
re_gv_ls = []
im_gv_ls = []

for tmin in tmin_ls:
    p_value, chi2, re_gv, im_gv = each_fit(tmin, b=1, z=1)
    p_value_ls.append(p_value)
    chi2_ls.append(chi2)
    re_gv_ls.append(re_gv)
    im_gv_ls.append(im_gv)


stability_plot(tmin_ls, re_gv_ls, p_value_ls, chi2_ls, 'b1_z1_real', chose_idx=0, save=True)



# %%
