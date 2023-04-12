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
    p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls = gs_fit.main(data_dic)

    re_gv = gv.dataset.avg_data(pdf_re_ls, bstrap=True)
    im_gv = gv.dataset.avg_data(pdf_im_ls, bstrap=True)

    return np.mean(p_value_ls), np.mean(chi2_ls), re_gv, im_gv


# def each_fit(tmin, b, z):
#     #! single gvar fit
#     gamma = 'z'
#     mass = 220
#     mom = 8
#     ll = 6


#     pt2_tmin = 3
#     pt2_tmax = 9
#     ra_tmin = tmin
#     ra_tmax = 9
#     tau_cut = 0


#     read_raw = Read_Raw('data_raw/')

#     data_dic = {}
#     temp_2pt = read_raw.read_2pt_bs(mass, mom)
#     data_dic['2pt_re'] = np.real( temp_2pt )
#     data_dic['2pt_im'] = np.imag( temp_2pt )


#     for tseq in range(4, 9):
#         temp_ra_re, temp_ra_im = read_raw.read_ratio_bs(gamma, mass, mom, ll, b, z, tseq=tseq)

#         data_dic['ra_re_tseq_{}'.format(tseq)] = temp_ra_re[:, 1:tseq]
#         data_dic['ra_im_tseq_{}'.format(tseq)] = temp_ra_im[:, 1:tseq]

#     data_dic_avg = gv.dataset.avg_data(data_dic, bstrap=True)

#     #* set the x values of fit
#     x = {}
    
#     # 2pt
#     x['2pt_re'] = np.arange(pt2_tmin, pt2_tmax)
#     x['2pt_im'] = np.arange(pt2_tmin, pt2_tmax)

#     # ratio
#     ra_t = []
#     ra_tau = []
#     for tseq in range(ra_tmin, ra_tmax):
#         for tau in range(1+tau_cut, tseq - tau_cut): #* because the tau in the data dic is from 1 to tseq - 1 without tseq, so tau_cut = 0 means tau from 1 to tseq - 1
#             ra_t.append(tseq)
#             ra_tau.append(tau)

#     x['ra_re'] = [ra_t, ra_tau]
#     x['ra_im'] = [ra_t, ra_tau]


#     #* set the y values of fit
#     y = gv.BufferDict()
#     y['2pt_re'] = data_dic_avg['2pt_re'][pt2_tmin:pt2_tmax]
#     y['2pt_im'] = data_dic_avg['2pt_im'][pt2_tmin:pt2_tmax]


#     ra_re = []
#     ra_im = []
#     for tseq in range(ra_tmin, ra_tmax):
#         for tau in range(1+tau_cut, tseq - tau_cut): #* because the tau in the data dic is from 1 to tseq - 1 without tseq, so tau_cut = 0 means tau from 1 to tseq - 1
#             tau_idx = tau - 1
#             ra_re.append(data_dic_avg['ra_re_tseq_{}'.format(tseq)][tau_idx])
#             ra_im.append(data_dic_avg['ra_im_tseq_{}'.format(tseq)][tau_idx])

#     y['ra_re'] = ra_re
#     y['ra_im'] = ra_im


#     gs_fit = Gs_Fit(two_state_fit(), fit_id='{}{}_P{}_L{}_b{}_z{}_tmax9_cut0'.format(mass, gamma, mom, ll, b, z))
#     gs_fit.para_set(pt2_tmin, pt2_tmax, ra_tmin, ra_tmax, tau_cut)

#     fit_res = lsf.nonlinear_fit(data=(x, y), prior=gs_fit.prior, fcn=gs_fit.get_fcn(), maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')


#     return fit_res.Q, fit_res.logGBF, fit_res.p['pdf_re'], fit_res.p['pdf_im']



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


# %%
stability_plot(tmin_ls, re_gv_ls, p_value_ls, chi2_ls, 'b1_z1_real', chose_idx=0, save=False)



# %%
