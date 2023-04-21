'''
This code is used to make the dispersion relation plot for the 2pt with different momentum.
'''

import lsqfit as lsf
# %%
import numpy as np
from funcs import *
from read_raw_module import Read_Raw

a = 0.12 # lattice spacing
Ls = 48

def meff_fit(t_ls, meff_ls):
    '''
    constant fit of meff
    input t_ls and meff_ls, return the gvar fit result of meff
    '''
    def fcn(x, p):
        return p['meff'] + x * 0

    priors = gv.BufferDict()
    priors['meff'] = gv.gvar(1, 10)

    fit_res = lsf.nonlinear_fit(data=(t_ls, meff_ls), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    return fit_res.p['meff']


#* read from the raw data
read_raw = Read_Raw('data_raw/')

data_dic = {}
for mass in [220, 310]:
    for mom in range(0, 14, 2):
        temp_2pt = np.real( read_raw.read_2pt_bs(mass=mass, mom=mom) )

        data_dic['m{}_P{}'.format(mass, mom)] = temp_2pt[:,1:] #* take the real part only, cut t=0 point

pt2_avg_dic = gv.dataset.avg_data(data_dic, bstrap=True)

meff_avg_dic = {}
for key in pt2_avg_dic.keys():
    meff_avg_dic[key] = pt2_to_meff(pt2_avg_dic[key])

#* make a effective mass plot
key = 'm{}_P{}'.format(310, 8)

errorbar_ls_plot([np.arange(len(meff_avg_dic[key]))], [gv.mean(meff_avg_dic[key])], [gv.sdev(meff_avg_dic[key])], ['meff'], title='meff', save=False)

# %%
#! meff fit parameters
tmin = 4
tmax = 8

m220_fit_ls = []
m310_fit_ls = []

for mom in range(0, 14, 2):
    m220_fit_ls.append( meff_fit(np.arange(tmin, tmax), meff_avg_dic['m{}_P{}'.format(220, mom)][tmin:tmax]) )

    m310_fit_ls.append( meff_fit(np.arange(tmin, tmax), meff_avg_dic['m{}_P{}'.format(310, mom)][tmin:tmax]) )

errorbar_ls_plot( [np.arange(0, 14, 2)], [gv.mean(m220_fit_ls)], [gv.sdev(m220_fit_ls)], ['m220'], title='m220 dispersion relation', save=False )

errorbar_ls_plot( [np.arange(0, 14, 2)], [gv.mean(m310_fit_ls)], [gv.sdev(m310_fit_ls)], ['m310'], title='m310 dispersion relation', save=False )


# %%
#todo make a nice plot of the dispersion relation
#todo convert lattice unit to GeV
#todo fit dispersion relation
#todo plot the fit results on the plot

def disp_relation_plot(mom_ls, meff_ls, title):
    a = 0.12 # lattice spacing in fm

    p_ls = lat_unit_convert(mom_ls, a, Ls, dimension='P')
    E_ls = lat_unit_convert(meff_ls, a, Ls, dimension='M')
    E2_ls = [v**2 for v in E_ls]

    def fcn(x, p):
        return np.sqrt( p['m']**2 + p['c1'] * x**2 + p['c2'] * x**4 * a**2 / (gev_fm**2) )

    
    priors = gv.BufferDict()
    priors['m'] = gv.gvar(0.1, 10)
    priors['c1'] = gv.gvar(1, 10)
    priors['c2'] = gv.gvar(0, 10)

    fit_res = lsf.nonlinear_fit(data=(p_ls, E_ls), prior=priors, fcn=fcn, maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')

    print(fit_res.format(100))

    fit_x = np.arange(p_ls[0], p_ls[-1], 0.01)
    fit_y = fcn(fit_x, fit_res.p)

    label = title[:4]

    fig = plt.figure(figsize=fig_size_Qi_An)
    ax = plt.axes(plt_axes)
    ax.errorbar(p_ls, [v.mean for v in E_ls], [v.sdev for v in E_ls], color='dodgerblue', label=label, **errorb)
    ax.fill_between(fit_x, [v.mean + v.sdev for v in fit_y], [v.mean - v.sdev for v in fit_y], color='dodgerblue', alpha=0.5)

    ax.tick_params(direction='out', **ls_p)
    ax.grid(linestyle=':')
    ax.set_xlabel(r'$P$ / GeV', font_Qi_An)
    ax.set_ylabel(r'$E$ / GeV', font_Qi_An)
    plt.legend(fontsize = fontsize_Qi_An)
    # plt.title(title, font_Qi_An)
    plt.savefig('fig/'+title+'.pdf', transparent=True)
    plt.show()

    return fit_res



fit_res = disp_relation_plot(np.arange(0, 14, 2), np.array(m220_fit_ls), 'm220 dispersion relation')
fit_res = disp_relation_plot(np.arange(0, 14, 2), np.array(m310_fit_ls), 'm310 dispersion relation')

# %%
