# %%
import matplotlib.pyplot as plt
import numpy as np
from funcs import *

grey = "#808080" 
red = "#FF6F6F" 
peach = "#FF9E6F" 
orange = "#FFBC6F" 
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C" 
turquoise = "#4AAB89"
blue = "#508EAD" 
grape = "#635BB1"
violet = "#7C5AB8" 
fuschia = "#C3559F"

color_ls = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green']

t_label = r'$\rm{t (a) }$'
meff_label = r'$m_{eff}$'
gev_fm = 0.1973269631 # 1 = 0.197 GeV . fm


def meff_plot(pt2_ls, ti, tf, fit_res, mom_ls, mom_plot, title):
    meff_ls = pt2_to_meff(pt2_ls)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    ax.errorbar(np.arange(len(meff_ls)), [val.mean for val in meff_ls], [val.sdev for val in meff_ls], color=orange, marker='D', label='meff', **errorb)

    t1_fit = np.linspace(ti, tf-1, 100)
    t2_fit = np.linspace(ti+1, tf, 100)

    x = {}
    x['proton'] = t1_fit

    c1_fit = fit_res.fcn( x, fit_res.p )['proton']

    x = {}
    x['proton'] = t2_fit

    c2_fit = fit_res.fcn( x, fit_res.p )['proton']

    meff_fit = []
    for i in range(100):
        meff = np.log( c1_fit[i] / c2_fit[i] )
        meff_fit.append(meff)

    ax.fill_between( t1_fit, [v.mean + v.sdev for v in meff_fit], [v.mean - v.sdev for v in meff_fit], color=blue, alpha=0.4, label='fit' )

    ax.set_ylim([0, 1])
    # ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(t_label, **fs_p)
    ax.set_ylabel(meff_label, **fs_p)
    ax.set_title(title, font)
    ax.legend(loc='upper right')
    ax.tick_params(direction='in', **ls_p)
    ax.grid(linestyle=':')
    plt.savefig('fig/'+title+'.pdf', transparent=True)
    plt.show()




def fit_on_data_R(data_set_tidy, mom, current, title, ylim=None):
    mo = '_' + str(mom)
    hash_key = 'p_sq_{}_pz_0'.format(mom)
    pt2_0_ls = data_set_tidy['p_sq_0_pz_0']['2pt']
    pt2_mom_ls = data_set_tidy[hash_key]['2pt']


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for tsep in range(3, 10):
        pt3_ls = data_set_tidy[hash_key][current+'_tsep_{}'.format(tsep)][1:-1]
        tau_ls = np.arange(tsep+1)[1:-1]
        R_tsep = pt2_pt3_to_R(tsep, tau_ls, pt2_0_ls, pt2_mom_ls, pt3_ls)

        ax.errorbar(tau_ls - tsep/2, [v.mean for v in R_tsep], [v.sdev for v in R_tsep], color=color_ls[tsep], label='tsep {}'.format(tsep), **errorb)


    ax.tick_params(direction='in', **ls_p)
    ax.grid(linestyle=':')
    ax.set_ylim(ylim)
    plt.title(title, font)
    plt.legend()
    plt.show()

    return
