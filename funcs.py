# %%
import numpy as np
import gvar as gv
import matplotlib.pyplot as plt


fig_width = 6.75 # in inches, 2x as wide as APS column
gr        = 1.618034333 # golden ratio
fig_size  = (fig_width, fig_width / gr)
plt_axes = [0.1,0.12,0.85,0.8] # 调整画图时 plot 的边框位置在画面中的比例，比如这里是左侧从 0.12 处开始，plot 占 画面的 0.8，所以右侧到 0.92 处，下侧从 0.15 处开始，上侧到 0.97 处

errorp = {"markersize": 5, "mfc": "none", "linestyle": "none"} # circle
errorb = {"markersize": 5, "mfc": "none", "linestyle": "none", "capsize": 3, "elinewidth": 1} # circle
errorl = {"markersize": 5, "mfc": "none", "capsize": 3, "elinewidth": 1} # circle with line
fs_p = {"fontsize": 13} # font size of text, label, ticks
ls_p = {"labelsize": 13}

font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 13}


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
    return meff_ls

def fit_on_data_plot(x, gv_y, fit_res, key, title, log_folder, ylim=None, save=True):
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

def pt2_pt3_to_R(tsep, tau_ls, pt2_0_ls, pt2_mom_ls, pt3_ls):
    #!# here pt2 list should be completed one started from tsep=0
    #!# tau list matches with pt3 list

    ratio1 = np.array( pt3_ls / pt2_0_ls[tsep] )

    #!# different from the formula 19 in the paper, coz we have zero momentum at sink
    ratio2 = []
    for tau in tau_ls:
        val1 = pt2_0_ls[tau] * pt2_0_ls[tsep] * pt2_mom_ls[tsep - tau]
        val2 = pt2_mom_ls[tau] * pt2_mom_ls[tsep] * pt2_0_ls[tsep - tau]

        ratio2.append( val1 / val2 )

    ratio2 = np.array(ratio2) 

    R_ls = ratio1 * (ratio2**0.5)

    return R_ls