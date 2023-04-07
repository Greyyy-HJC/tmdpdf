# %%
import numpy as np
import gvar as gv
import matplotlib.pyplot as plt
from head import *


def bootstrap(conf_ls, times=500, seed_path=None):
    '''
    make sure conf_ls.shape = (N_conf, ...)
    return conf_ls
    '''

    #*# If using a fixed bs_ls
    if seed_path is not None:
        bs_ls = gv.load(seed_path)
        conf_bs = np.mean(conf_ls[bs_ls], axis=1)

    #*# If generating a random bs_ls
    else:
        idx_ls = np.random.randint(len(conf_ls), size=(times, len(conf_ls)))
        conf_bs = np.mean(conf_ls[idx_ls], axis=2)

    return conf_bs


def bs_ls_to_gvar_ls(bs_ls):
    '''
    This function is used to convert the bootstrap list to gvar list by combining each sample with the sdev of all samples

    The shape of bs_ls should be (N_samp, ...)
    '''

    avg = gv.dataset.avg_data(bs_ls, bstrap=True)
    sdev = gv.sdev(avg)

    return [gv.gvar(v, sdev) for v in bs_ls]


def jackknife(data):
    '''
    make sure data.shape = (N_conf * n_t)
    '''
    nf, nt = data.shape # data shape: (N_conf * n_t)
    cv = np.mean(data, axis=0, keepdims=True) # average all conf, cv shape: (1 * nt)
    jac = (nf * cv - data) / (nf - 1) # drop one data each time then average: (mean[N,:] * N_conf - data[N,:]) / (N_conf-1) 
    return jac # jac shape: (n_conf * n_t)


def jk_conf_avg(conf_ls):
    N_conf = len(conf_ls)
    mean = np.mean(conf_ls, axis=0)
    cov = np.cov(conf_ls, rowvar=False) * (N_conf - 1)

    return gv.gvar(mean, cov)

def jk_dic_avg(dic):
    l_dic = {}
    for key in dic:
        l_dic[key] = len(dic[key][0])

    conf_ls = []
    for n_conf in range(len(dic[key])):
        temp = []
        for key in dic:
            temp.append(list(dic[key][n_conf]))

        conf_ls.append( sum(temp, []) ) ## flat

    gv_ls = list(jk_conf_avg(conf_ls))
    
    gv_dic = {}
    for key in l_dic:
        gv_dic[key] = []
        for i in range(l_dic[key]):
            temp = gv_ls.pop(0)
            gv_dic[key].append(temp)

    return gv_dic

def bs_conf_avg(conf_ls):
    mean = np.mean(conf_ls, axis=0)
    cov = np.cov(conf_ls, rowvar=False)

    return gv.gvar(mean, cov)

def bs_dic_avg(dic):
    l_dic = {}
    for key in dic:
        l_dic[key] = len(dic[key][0])

    conf_ls = []
    for n_conf in range(len(dic[key])):
        temp = []
        for key in dic:
            temp.append(list(dic[key][n_conf]))

        conf_ls.append( sum(temp, []) ) ## flat

    gv_ls = list(bs_conf_avg(conf_ls))
    
    gv_dic = {}
    for key in l_dic:
        gv_dic[key] = []
        for i in range(l_dic[key]):
            temp = gv_ls.pop(0)
            gv_dic[key].append(temp)
        
        gv_dic[key] = np.array(gv_dic[key])

    return gv_dic

def gv_to_samples_corr(gv_ls, N_samp):
    '''
    transform gvar to bs samples with correlation
    shape = (N_samp, len(ls))
    '''
    mean = [v.mean for v in gv_ls]
    cov_m = gv.evalcov(gv_ls)
    rng = np.random.default_rng()

    samp_ls = rng.multivariate_normal(mean, cov_m, size=N_samp)

    return samp_ls

def gv_dic_to_samples_corr(gv_dic, N_samp):
    l_dic = {} # record length
    for key in gv_dic:
        l_dic[key] = len(gv_dic[key])

    flatten_ls = []
    for key in gv_dic:
        flatten_ls.append(list(gv_dic[key]))

    flatten_ls = sum(flatten_ls, []) ## flat

    samp_all = gv_to_samples_corr(flatten_ls, N_samp)
    samp_all = list(np.swapaxes(samp_all, 0, 1)) # shape = len(all), 1000

    samp_dic = {}
    for key in l_dic:
        samp_ls = []
        for i in range(l_dic[key]):
            temp = samp_all.pop(0)
            samp_ls.append(temp)

        samp_ls = np.swapaxes(np.array(samp_ls), 0, 1) # shape = 1000, len(key)
        samp_dic[key] = samp_ls

    return samp_dic

def errorbar_plot(x, y, yerr, title, ylim=None, save=True):
    #* this is a general plot function, so put it here

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(x, y, yerr, **errorb)
    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')
    ax.set_ylim(ylim)
    plt.title(title, font)
    # plt.legend()
    if save == True:
        plt.savefig('fig/'+title+'.pdf', transparent=True)
    # plt.show()

def errorbar_ls_plot(x_ls, y_ls, yerr_ls, title, ylim=None, save=True):
    #* this is a general plot function

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for x, y, yerr in zip(x_ls, y_ls, yerr_ls):
        ax.errorbar(x, y, yerr, **errorb)
    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')
    ax.set_ylim(ylim)
    plt.title(title, font)
    # plt.legend()
    if save == True:
        plt.savefig('fig/'+title+'.pdf', transparent=True)
    # plt.show()


def fill_between_plot(x, y, yerr, title, ylim=None, save=True):
    #* this is a general plot function, so put it here

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.fill_between(x, [y[i]+yerr[i] for i in range(len(y))], [y[i]-yerr[i] for i in range(len(y))], alpha=0.4)
    ax.tick_params(direction='in', top='on', right='on', **ls_p)
    ax.grid(linestyle=':')
    ax.set_ylim(ylim)
    plt.title(title, font)
    # plt.legend()
    if save == True:
        plt.savefig('fig/'+title+'.pdf', transparent=True)
    # plt.show()


def pt2_to_meff(pt2_ls):
    meff_ls = []
    for i in range(len(pt2_ls)-1):
        val = np.log(pt2_ls[i]) - np.log(pt2_ls[i+1])
        meff_ls.append(val)
    return np.array(meff_ls)

