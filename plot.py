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


def fit_on_data_plot_2pt(x, gv_y, fit_res, key, title, log_folder, ylim=None, save=True):
    #* this is a general plot function to make effective mass plot with fit band on the data

    #* data part
    mx = x[:-1]
    gv_my = pt2_to_meff(gv_y)
    my = [v.mean for v in gv_my]
    myerr = [v.sdev for v in gv_my]

    #* fit part
    fit_x = np.linspace(mx[0], mx[-1], 100)
    input_x = {}
    input_x[key] = fit_x
    fit_y = fit_res.fcn(input_x, fit_res.p)[key]
    fit_mx = []
    fit_my = []
    fit_myerr = []

    for i in range(len(fit_x)-1):
        val = np.log(fit_y[i]) - np.log(fit_y[i+1])
        val = val / (fit_x[i+1] - fit_x[i])
        fit_mx.append(fit_x[i])
        fit_my.append(val.mean)
        fit_myerr.append(val.sdev)


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(mx, my, myerr, label='data', **errorb)
    ax.fill_between(fit_mx, [fit_my[i]+fit_myerr[i] for i in range(len(fit_my))], [fit_my[i]-fit_myerr[i] for i in range(len(fit_my))], alpha=0.4, label='fit')
    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')
    ax.set_ylim(ylim)
    plt.title(title, font)
    plt.legend()
    if save == True:
        plt.savefig(log_folder+title+'.pdf', transparent=True)
    # plt.show()




    

def fit_on_data_plot_ratio(ra_t, ra_tau, ra_re_gv, ra_im_gv, fit_res, title, log_folder): 
    '''
    This function is used to make a plot of the 3pt / 2pt ratio with fit results on the data points
    Plot both real and imag parts
    ra_t and ra_tau are just t_ls and tau_ls used for fits, ra_re_gv and ra_im_gv are the y values for fits
    '''

    tmin = min(ra_t)
    tmax = max(ra_t) + 1

    input_x = {}
    input_x['2pt_re'] = np.arange(3, 10)
    input_x['2pt_im'] = np.arange(3, 10)
    input_x['ra_re'] = [ra_t, ra_tau]
    input_x['ra_im'] = [ra_t, ra_tau]

    fit_re_val = fit_res.fcn(input_x, fit_res.p)['ra_re']
    fit_im_val = fit_res.fcn(input_x, fit_res.p)['ra_im']

    tau_dic = {}
    ra_re_dic = {}
    ra_im_dic = {}
    fit_re_dic = {}
    fit_im_dic = {}

    for idx in range(len(ra_t)):
        key = 'tseq_{}'.format(ra_t[idx])
        if key not in tau_dic:
            tau_dic[key] = []
            ra_re_dic[key] = []
            ra_im_dic[key] = []
            fit_re_dic[key] = []
            fit_im_dic[key] = []

        tau_dic[key].append(ra_tau[idx])
        ra_re_dic[key].append(ra_re_gv[idx])
        ra_im_dic[key].append(ra_im_gv[idx])
        fit_re_dic[key].append(fit_re_val[idx])
        fit_im_dic[key].append(fit_im_val[idx])



    #* plot real part
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for tseq in range(tmin, tmax):
        key = 'tseq_{}'.format(tseq)
        ax.errorbar(np.array( tau_dic[key] ) - tseq/2, [v.mean for v in ra_re_dic[key]], [v.sdev for v in ra_re_dic[key]], label='tseq = {}'.format(tseq), color=color_ls[tseq - tmin], **errorb)

        ax.fill_between(np.array( tau_dic[key] ) - tseq/2, [v.mean+v.sdev for v in fit_re_dic[key]], [v.mean-v.sdev for v in fit_re_dic[key]], color=color_ls[tseq - tmin], alpha=0.4)

    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')

    plt.title(title+'_real', font)
    plt.legend(ncol=3)
    plt.xlabel(r'$\tau - t/2$', font)
    plt.ylabel(r'g.s.', font)
    plt.savefig(log_folder+title+'_real.pdf', transparent=True)
    plt.show()


    #* plot imag part
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for tseq in range(tmin, tmax):
        key = 'tseq_{}'.format(tseq)
        ax.errorbar(np.array( tau_dic[key] ) - tseq/2, [v.mean for v in ra_im_dic[key]], [v.sdev for v in ra_im_dic[key]], label='tseq = {}'.format(tseq), color=color_ls[tseq - tmin], **errorb)

        ax.fill_between(np.array( tau_dic[key] ) - tseq/2, [v.mean+v.sdev for v in fit_im_dic[key]], [v.mean-v.sdev for v in fit_im_dic[key]], color=color_ls[tseq - tmin], alpha=0.4)

    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')

    plt.title(title+'_imag', font)
    plt.legend(ncol=3)
    plt.xlabel(r'$\tau - t/2$', font)
    plt.ylabel(r'g.s.', font)
    plt.savefig(log_folder+title+'_imag.pdf', transparent=True)
    plt.show()
