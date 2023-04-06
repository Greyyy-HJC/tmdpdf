# %%
import numpy as np
import gvar as gv
import lsqfit as lsf
from fit_module import Fit
from prior_setting import *
from funcs import *
from plot import *


ens = 'a09m310'
data_set_tidy = gv.load('dump/data_set_tidy')
# print([key for key in data_set_tidy])


#!# check meff plot
def data_check_meff(data_set_tidy, hadron):
    pt2_ls = data_set_tidy[hadron]
    t_ls = np.arange(len(pt2_ls))
    meff_ls = pt2_to_meff(pt2_ls)

    x_ls = [t_ls[:-1]]
    y_ls = [[v.mean for v in meff_ls]]
    yerr_ls = [[v.sdev for v in meff_ls]]

    errorbar_plot(x_ls, y_ls, yerr_ls, '{}_meff'.format(hadron), ylim=[0, 1.5], save=False)
    plt.show()

    meff = np.mean(meff_ls[6:12])


    # print('meff p_sq = {}'.format(p_sq))
    # print(meff)

    return meff

# data_check_meff(data_set_tidy, hadron='proton')


#!# check ratio plot
# fit_on_data_R(data_set_tidy, mom=2, current='V4', title='fit on data R mom 2', ylim=None)


#!# check disp relation



#!# fit 
def do_fit(mom_ls, pt2_n, pt3_n):
    log_folder = 'fit_log/{}/'.format(ens)

    a09m310_fit = Fit(prior_ho_a09m310, pt2_n, pt3_n, include_2pt=True, include_3pt=False)

    #* 2pt variable t
    pt2_t = {}
    pt2_t['proton'] = np.arange(3, 20)


    #* 3pt variable t and tau
    t_ls = []
    tau_ls = []
    tau_cut = 2
    for t in range(4, 10):
        for tau in range(tau_cut, t+1-tau_cut):
            t_ls.append(t)
            tau_ls.append(tau)
    pt3_A3 = {}
    pt3_A3['proton'] = [t_ls, tau_ls]
  

    t_ls = []
    tau_ls = []
    tau_cut = 2
    for t in range(4, 10):
        for tau in range(tau_cut, t+1-tau_cut):
            t_ls.append(t)
            tau_ls.append(tau)
    pt3_V4 = {}
    pt3_V4['proton'] = [t_ls, tau_ls]

    #* add p0 if necessary
    # temp = gv.load('dump/post')
    # p0 = {}
    # for key in temp:
    #     p0[key] = temp[key].mean

    fit_res, corr = a09m310_fit.fit(data_set_tidy, pt2_t, pt3_A3, pt3_V4, mom_ls, best_p0=None, corr=False)
    post = fit_res.p

    #!# record fit results into a txt log file
    log = open(log_folder+"{}.txt".format('proton'), mode="w", encoding="utf-8")
    print(fit_res.format(100), file=log)
    log.close()

    print('Q = {}'.format(fit_res.Q))

    return fit_res



pt2_n = 3
pt3_n = 3
mom_ls = [0]
fit_res = do_fit(mom_ls, pt2_n, pt3_n)

pt2_ls = data_set_tidy['proton'][:20]
mom_plot = 0

meff_plot(pt2_ls, ti=3, tf=15, fit_res=fit_res, mom_ls=mom_ls, mom_plot=mom_plot, title='proton_meff_fit_on_data')

