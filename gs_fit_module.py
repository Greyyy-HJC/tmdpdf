'''
This is the module to do the ground state fit with bootstrap data.

The input will be the prior setting dict.
The ouput will be the lists of p value, chi2, and the real / imag bare g.s. matrix elements, the lists will be in the same order as the input bs list.

All the fitting parameters will be set here with function para_set.
Call the main function with data dic as the input, it will return the fit result.

Here we do the joint fit of 2pt and ratio, FH can be added if necessary.
Here we set pt2_n = ra_n = 2.

Data dic should be a dict as:
data_dic['2pt_re'] = shape (N_bs, N_t)
data_dic['2pt_im'] = shape (N_bs, N_t)
data_dic['ra_re_tseq_{}'] = shape (N_bs, N_tau), here tau should start from 1 to tseq(not included)
data_dic['ra_im_tseq_{}'] = shape (N_bs, N_tau)

Both the fit results and the plot of the first fit with bs_idx = 0 will be saved in log folder. The log file named as '220z_P8_L6b1z1'.

Example usage can be found in the end.
'''

# %%
import numpy as np
import gvar as gv
import lsqfit as lsf
import os
from tqdm import tqdm

from prior_setting import two_state_fit
from funcs import *
from plot import *

class Gs_Fit():
    def __init__(self, prior_dic, fit_id):
        self.prior = prior_dic
        self.fit_id = fit_id #* the id of the fit, will be used to save the log, should be a string like '220z_P8_L6b1z1'

    def para_set(self, pt2_tmin, pt2_tmax, ra_tmin, ra_tmax, tau_cut=0):
        self.pt2_tmin = pt2_tmin
        self.pt2_tmax = pt2_tmax
        self.ra_tmin = ra_tmin
        self.ra_tmax = ra_tmax
        self.tau_cut = tau_cut

    def main(self, data_dic):
        N_bs = len(data_dic['2pt_re']) # number of bootstrap samples

        #* set the x values of fit
        x = {}
        
        # 2pt
        x['2pt_re'] = np.arange(self.pt2_tmin, self.pt2_tmax)
        x['2pt_im'] = np.arange(self.pt2_tmin, self.pt2_tmax)

        # ratio
        ra_t = []
        ra_tau = []
        for tseq in range(self.ra_tmin, self.ra_tmax):
            for tau in range(1+self.tau_cut, tseq - self.tau_cut): #* because the tau in the data dic is from 1 to tseq - 1 without tseq, so tau_cut = 0 means tau from 1 to tseq - 1
                ra_t.append(tseq)
                ra_tau.append(tau)

        x['ra_re'] = [ra_t, ra_tau]
        x['ra_im'] = [ra_t, ra_tau]


        #* transform the bs list to gvar list
        data_dic_gv = {}
        for key in data_dic:
            data_dic_gv[key] = np.array( bs_ls_to_gvar_ls(data_dic[key]) )

        #! bootstrap loop fit
        p_value_ls = []
        chi2_ls = []
        pdf_re_ls = []
        pdf_im_ls = []

        print('\n>>> Start bs g.s. fits of {}: \n'.format(self.fit_id))
        for bs_idx in tqdm( range(N_bs) ):
            #* set the y values of fit
            y = gv.BufferDict()

            y['2pt_re'] = data_dic_gv['2pt_re'][bs_idx, self.pt2_tmin:self.pt2_tmax]

            ra_re = []
            ra_im = []
            for tseq in range(self.ra_tmin, self.ra_tmax):
                for tau in range(1+self.tau_cut, tseq - self.tau_cut): #* because the tau in the data dic is from 1 to tseq - 1 without tseq, so tau_cut = 0 means tau from 1 to tseq - 1
                    tau_idx = tau - 1 #! this because the tau=0 and tau=tseq in the data dic has already been cut

                    ra_re.append(data_dic_gv['ra_re_tseq_{}'.format(tseq)][bs_idx, tau_idx])
                    ra_im.append(data_dic_gv['ra_im_tseq_{}'.format(tseq)][bs_idx, tau_idx])

            y['ra_re'] = ra_re
            y['ra_im'] = ra_im

            fit_res = lsf.nonlinear_fit(data=(x, y), prior=self.prior, fcn=self.get_fcn(), maxit=10000, svdcut=1e-100, fitter='scipy_least_squares')


            #todo bad fits warning
            if fit_res.chi2 / fit_res.dof > 2:
                print('>>> Warning: bad fit for bs_idx = {} in fit {}'.format(bs_idx, self.fit_id))
                print('>>> chi2/dof = {}'.format(fit_res.chi2 / fit_res.dof))


            ### * ### for log
            if bs_idx == 0:
                #* res file path
                log_folder = 'log/gs_fit/{}/'.format(self.fit_id)
                if not os.path.exists(log_folder):
                    os.mkdir(log_folder)

                #* save the first fit result
                log = open(log_folder+"bs_idx=0.txt", mode="w", encoding="utf-8")
                print(fit_res.format(100), file=log)
                log.close()

                #* add plot
                ra_re_gv = np.array(ra_re)
                ra_im_gv = np.array(ra_im)
                title = self.fit_id + '_bs_idx=0'

                fit_on_data_plot_ratio(ra_t, ra_tau, ra_re_gv, ra_im_gv, fit_res, title, log_folder)


            p_value_ls.append(fit_res.Q)
            chi2_ls.append(fit_res.chi2 / fit_res.dof)
            pdf_re_ls.append(fit_res.p['pdf_re'].mean) #* only save the mean value of the fit result
            pdf_im_ls.append(fit_res.p['pdf_im'].mean)

            #* save the p value and chi2 for all bs fit results
            log = open(log_folder+"bs_collection.txt", mode="w", encoding="utf-8")
            print('\n>>> p values: \n', file=log)
            print(p_value_ls, file=log)
            print('\n>>> chi2: \n', file=log)
            print(chi2_ls, file=log)
            print('\n>>> real: \n', file=log)
            print(pdf_re_ls, file=log)
            print('\n>>> imag: \n', file=log)
            print(pdf_im_ls, file=log)
            log.close()

            ### * ###



        return p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls


    def pt2_re_fcn(self, pt2_t, p):
        #! checked to be consistent with the paper
        de = p['dE1']

        val = p['re_c0'] * np.exp(-p['E0']*pt2_t) * (1 + p['re_c1'] * np.exp( -de * pt2_t ))

        return val

    def pt2_im_fcn(self, pt2_t, p):
        de = p['dE1']

        val = p['im_c0'] * np.exp(-p['E0']*pt2_t) * (1 + p['im_c1'] * np.exp( -de * pt2_t ))

        return val

    def ra_re_fcn(self, ra_t, ra_tau, p):
        de = p['dE1']

        numerator = p['pdf_re'] + p['re_c2'] * ( np.exp( -de * (ra_t - ra_tau) ) + np.exp( -de * ra_tau ) ) + p['re_c3'] * np.exp(-de * ra_t)
        val = numerator / ( 1 + p['re_c1'] * np.exp( -de * ra_t ) )

        return val

    def ra_im_fcn(self, ra_t, ra_tau, p):
        de = p['dE1']

        numerator = p['pdf_im'] + p['im_c2'] * ( np.exp( -de * (ra_t - ra_tau) ) + np.exp( -de * ra_tau ) ) + p['im_c3'] * np.exp(-de * ra_t)
        val = numerator / (1 + p['im_c1'] * np.exp( -de * ra_t ))

        return val

    def get_fcn(self):
        #* x = [2pt_re, 2pt_im, ra_re, ra_im], ra_re = [ra_t, ra_tau]
        #* ra_t like [3, 3, 4, 4, 4, 5, 5, 5, 5, ...]
        #* ra_tau like [1, 2, 1, 2, 3, 1, 2, 3, 4, ...]
        def fcn(x, p):
            val = {}
            val['2pt_re'] = self.pt2_re_fcn(x['2pt_re'], p)
            val['2pt_im'] = self.pt2_im_fcn(x['2pt_im'], p)

            val['ra_re'] = []
            val['ra_im'] = []
            for idx in range(len(x['ra_re'][0])):
                tsep = x['ra_re'][0][idx]
                tau = x['ra_re'][1][idx]
                val['ra_re'].append( self.ra_re_fcn(tsep, tau, p) )

                tsep = x['ra_im'][0][idx]
                tau = x['ra_im'][1][idx]
                val['ra_im'].append( self.ra_im_fcn(tsep, tau, p) )

            return val
        return fcn

if __name__ == '__main__':
    from read_raw_module import Read_Raw

    read_raw = Read_Raw('data_raw/')

    data_dic = {}
    temp_2pt = read_raw.read_2pt_bs(mass=220, mom=8)
    data_dic['2pt_re'] = np.real( temp_2pt )
    data_dic['2pt_im'] = np.imag( temp_2pt )


    for tseq in range(4, 8):
        temp_ra_re, temp_ra_im = read_raw.read_ratio_bs(gamma='z', mass=220, mom=8, ll=6, b=1, z=1, tseq=tseq)

        data_dic['ra_re_tseq_{}'.format(tseq)] = temp_ra_re[:, 1:tseq]
        data_dic['ra_im_tseq_{}'.format(tseq)] = temp_ra_im[:, 1:tseq]



    gs_fit = Gs_Fit(two_state_fit(), fit_id='220z_P8_L6b1z1')
    gs_fit.para_set(pt2_tmin=3, pt2_tmax=10, ra_tmin=4, ra_tmax=8, tau_cut=0)
    p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls = gs_fit.main(data_dic)
    


# %%
'''
Old fit results (totally 624000 fits)
    bs_800, tseq = 4, 5, 6, 7, tau_cut = 0

chi2 > 2:
    total: 12617
    b2_z6(include): 355

chi2 > 1.6:
    total: 49042
    b2_z6(include): 1194

    
New fit results (totally 624000 fits)
    bs_800
    within b2_z6: tseq = 4, 5, 6, 7, 8 tau_cut = 0

'''