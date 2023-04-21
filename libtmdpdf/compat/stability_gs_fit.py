'''
This code is used to make stability plot for the ground state fits.
Pick some (b, z) points, vary the tmin of the fit, and plot the comparison plot.
'''

# %%
import os

import lsqfit as lsf
import numpy as np
from funcs import *
from gs_fit_module import Gs_Fit
from prior_setting import two_state_fit
from read_raw_module import Read_Raw

'''
define a stability plot function with varying tmin, 3 subplots: matrix element, p value, logGBF

can be used to plot both real and imaginary parts

The input should be tmin list, gvar list, Q list, logGBF list, chose_idx.

chose_idx is the index in the x list, which indicates the fit that you choose to use.
'''
def tmin_stability_plot(tmin_ls, gv_y_ls, Q_ls, logGBF_ls, title, chose_idx, save=True):
    
    #* Create the subplots and set the height ratios
    fig, axs = plt.subplots(3, 1, sharex=True, figsize=fig_size_Qi_An, gridspec_kw=gridspec_tmin)

    label = r'$b = z = 1 a$'

    #* Plot the data on each subplot
    axs[0].errorbar(tmin_ls, [v.mean for v in gv_y_ls], [v.sdev for v in gv_y_ls], label=label, color='dodgerblue', **errorb)

    #* Plot the chosen fit
    upper = gv_y_ls[chose_idx].mean + gv_y_ls[chose_idx].sdev
    lower = gv_y_ls[chose_idx].mean - gv_y_ls[chose_idx].sdev

    axs[0].fill_between(tmin_ls, np.ones_like(tmin_ls) * upper, np.ones_like(tmin_ls) * lower, color=grey, alpha=0.4)

    axs[1].scatter(tmin_ls, Q_ls, marker='X', facecolors='none', edgecolors='k', s=20)
    axs[1].plot(tmin_ls, 0.1 * np.ones_like(tmin_ls), 'r--', linewidth=1)
    axs[2].scatter(tmin_ls, logGBF_ls, marker='o', facecolors='none', edgecolors='k', s=20)


    # Add labels to the x- and y-axes
    axs[0].set_ylabel(r'$\tilde{h}_{\Gamma}^{0}$', font_Qi_An)
    axs[1].set_ylabel(r'$Q$', font_Qi_An)
    axs[2].set_ylabel(r'$\rm logGBF$')

    # set the ylim of the first subplot to be the mean +- 6 sigma
    mid = (upper + lower) / 2
    sigma = (upper - lower) / 2
    axs[0].set_ylim(mid-6*sigma, mid+6*sigma)
    axs[1].set_ylim(-0.3, 1.1)

    # set the ylim of the last subplot to be the avg +- 3 gap
    avg = np.mean(logGBF_ls)
    gap = np.max(logGBF_ls) - np.min(logGBF_ls)
    axs[2].set_ylim(avg-3*gap, avg+3*gap)

    for i in range(3):
        axs[i].tick_params(direction='out', **ls_p)
        axs[i].grid(linestyle=':')

    plt.subplots_adjust(hspace=0)
    # axs[0].set_title(title, font)
    axs[2].set_xlabel(r'$t_{\rm min} / a$', font_Qi_An)

    axs[2].set_xticks(np.arange(4, 7, 1))
    axs[2].set_xticklabels(np.arange(4, 7, 1, dtype=int))

    axs[0].legend(loc='upper left', fontsize = fontsize_Qi_An)

    if save == True:
        plt.savefig('fig/'+title+'.pdf', transparent=True)
    # Display the plot
    plt.show()


def each_fit(tmin, gamma, mass, mom, ll , b, z):
    #! bs fits

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
    logGBF = fit_res.logGBF
    re_gv = fit_res.p['pdf_re']
    im_gv = fit_res.p['pdf_im']


    return p_value, logGBF, re_gv, im_gv





gamma = 'z'
mass = 220
mom = 8
ll = 6
b = 1
z = 1

tmin_ls = np.arange(4, 7)
p_value_ls = []
logGBF_ls = []
re_gv_ls = []
im_gv_ls = []

for tmin in tmin_ls:
    p_value, logGBF, re_gv, im_gv = each_fit(tmin, gamma, mass, mom, ll, b, z)
    p_value_ls.append(p_value)
    logGBF_ls.append(logGBF)
    re_gv_ls.append(re_gv)
    im_gv_ls.append(im_gv)


tmin_stability_plot(tmin_ls, re_gv_ls, p_value_ls, logGBF_ls, '{}{}_P{}_L{}_b{}_z{}_real'.format(mass, gamma, mom, ll, b, z), chose_idx=0, save=True)

tmin_stability_plot(tmin_ls, im_gv_ls, p_value_ls, logGBF_ls, '{}{}_P{}_L{}_b{}_z{}_imag'.format(mass, gamma, mom, ll, b, z), chose_idx=0, save=True)

# %%
