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




    

def fit_on_data_plot_ratio(x_ls, gv_y_ls, fit_res, re_im, title, log_folder, ylim=None): 
    #todo
    '''
    Still need to work on this!
    Plot both real and imag
    better way to get input
    '''


    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    #* data part
    for x, gv_y, color in zip(x_ls, gv_y_ls, color_ls):
        ax.errorbar(x, [v.mean for v in gv_y], [v.sdev for v in gv_y], color=color, **errorb)

    #* fit part
    input_x = {} # as the input to the fcn, so that to get the band of fit results
    input_x['2pt_re'] = np.arange(3, 10)
    input_x['2pt_im'] = np.arange(3, 10)

    ra_t = []
    ra_tau = []
    for tseq in range(4, 8):
        for tau in range(1, tseq):
            ra_t.append(tseq)
            ra_tau.append(tau)

    input_x['ra_re'] = [ra_t, ra_tau]
    input_x['ra_im'] = [ra_t, ra_tau]

    temp = fit_res.fcn(input_x, fit_res.p)['ra_'+re_im]


    fit_y_ls = []
    for tseq in range(4, 8):
        fit_y_ls.append([])
        for tau in range(1, tseq):
            fit_y_ls[tseq-4].append(temp[0])
            temp.pop(0)

    for x, fit_y, color in zip(x_ls, fit_y_ls, color_ls):
        ax.fill_between(x, [v.mean+v.sdev for v in fit_y], [v.mean-v.sdev for v in fit_y], alpha=0.4, color=color)


    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')

    plt.title(title, font)
    plt.legend(ncol=3)
    plt.ylim(ylim)
    plt.xlabel(r'$\tau - t/2$', font)
    plt.ylabel(r'g.s.', font)
    plt.savefig(log_folder+title+'.pdf', transparent=True)
    plt.show()

