import numpy as np 
import gvar as gv  
import lsqfit as lsf

class Fit():
    def __init__(self, prior, pt2_n, pt3_n, include_2pt=True, include_3pt=True):
        self.pt2_n = pt2_n
        self.pt3_n = pt3_n
        
        self.include_2pt = include_2pt
        self.include_3pt = include_3pt
        
        self.prior = prior(self.pt2_n, self.pt3_n)

    def pt2_fit_func(self, pt2_t, p, mom):
        pt2_n = self.pt2_n
        mo = '_'+str(mom) # to indicate the value of momentum

        E_list = {}
        for i in range(pt2_n): # initialize       
            E_list['E'+str(i)] = p['E0'+mo]

        for i in range(1, pt2_n): # define Ei      
            for j in range(1, i+1):
                if i == pt2_n-1 and j == i:
                    E_list['E'+str(i)] += p['dEmax_pt2'+mo]
                else:
                    E_list['E'+str(i)] += p['dE'+str(j)+mo]

        val = p['z0'+mo]*p['z0'+mo]*np.exp(-E_list['E0']*pt2_t)
        
        if pt2_n == 1:
            return val

        for i in range(1, pt2_n):
            val += p['z'+str(i)+mo]*p['z'+str(i)+mo]*np.exp(-E_list['E'+str(i)]*pt2_t)

        return val

    def pt3_fit_func(self, pt3_t_A3, pt3_tau_A3, pt3_t_V4, pt3_tau_V4, p, mom):
        #!# the momentum at the sink is fixed to zero in all three- point functions.
        pt3_n = self.pt3_n
        mo = '_'+str(mom) # to indicate the value of momentum

        E_list = {}
        E_list_0 = {}
        for i in range(pt3_n): # initialize       
            E_list['E'+str(i)] = p['E0'+mo]
            E_list_0['E'+str(i)] = p['E0_0'] #* for zero momentum

        for i in range(1, pt3_n): # define Ei    
            for j in range(1, i+1):
                if i == pt3_n-1 and j == i:
                    if self.pt3_n == self.pt2_n: #* if pt2_n = pt3_n, use same trash can
                        E_list['E'+str(i)] += p['dEmax_pt2'+mo]
                        E_list_0['E'+str(i)] += p['dEmax_pt2'+'_0']
                    else:
                        E_list['E'+str(i)] += p['dEmax_pt3'+mo]
                        E_list_0['E'+str(i)] += p['dEmax_pt3'+'_0']
                else:
                    E_list['E'+str(i)] += p['dE'+str(j)+mo]
                    E_list_0['E'+str(i)] += p['dE'+str(j)+'_0']


        val = {}
        val['pt3_A3'] = 0
        val['pt3_V4'] = 0

        for i in range(self.pt3_n):    
            for j in range(self.pt3_n):
                if mom == 0: 
                    mi = str(np.minimum(j, i))
                    ma = str(np.maximum(j, i))
                else: #!# non-zero mom's A_ij is not symmetric
                    mi = str(i)
                    ma = str(j)

                val['pt3_A3'] += p['A3_'+mi+ma+mo]*p['z'+str(j)+mo]*p['z'+str(i)+'_0'] * np.exp( - E_list['E'+str(i)]*pt3_tau_A3 - E_list_0['E'+str(j)]*(pt3_t_A3-pt3_tau_A3) ) # exp(-Ei*tau - Ej*(t-tau))

                val['pt3_V4'] += p['V4_'+mi+ma+mo]*p['z'+str(j)+mo]*p['z'+str(i)+'_0'] * np.exp( - E_list['E'+str(i)]*pt3_tau_V4 - E_list_0['E'+str(j)]*(pt3_t_V4-pt3_tau_V4) ) # exp(-Ei*tau - Ej*(t-tau))

        return val

    def fit_func(self, mom_ls):
        def fcn(x, p):
            val = {}
            
            if self.include_2pt == True:
                #* mom = 0
                pt2_t = x['proton']
                val['proton'] = self.pt2_fit_func(pt2_t, p, 0)

                for mom in mom_ls:
                    pass

            if self.include_3pt == True:
                #* mom = 0
                pt3_t_A3 = x['pt3_A3_0'][0]
                pt3_tau_A3 = x['pt3_A3_0'][1]
                pt3_t_V4 = x['pt3_V4_0'][0]
                pt3_tau_V4 = x['pt3_V4_0'][1]


                temp = self.pt3_fit_func(pt3_t_A3, pt3_tau_A3, pt3_t_V4, pt3_tau_V4, p, 0)

                val['pt3_A3_0'] = temp['pt3_A3']
                val['pt3_V4_0'] = temp['pt3_V4']

                for mom in mom_ls:
                    mo = '_'+str(mom)

                    pt3_t_A3 = x['pt3_A3'+mo][0]
                    pt3_tau_A3 = x['pt3_A3'+mo][1]
                    pt3_t_V4 = x['pt3_V4'+mo][0]
                    pt3_tau_V4 = x['pt3_V4'+mo][1]


                    temp = self.pt3_fit_func(pt3_t_A3, pt3_tau_A3, pt3_t_V4, pt3_tau_V4, p, mom)

                    val['pt3_A3'+mo] = temp['pt3_A3']
                    val['pt3_V4'+mo] = temp['pt3_V4']

            return val
        return fcn


    def fit(self, data_dic, pt2_t, pt3_A3, pt3_V4, mom_ls, best_p0=None, corr=False, priors=None, save=False):
        if priors == None:
            priors = self.prior

        t_tsep_tau = {}
        Amp = {}

        if self.include_2pt == True:
            t_tsep_tau['proton'] = pt2_t['proton']

            pt2_amp = []
            for t in pt2_t['proton']:
                pt2_amp.append(data_dic['proton'][t])

            Amp['proton'] = np.array(pt2_amp)

            #* mom != 0


        if self.include_3pt == True:
            pass

        if best_p0 == None:
            fit_result = lsf.nonlinear_fit(data=(t_tsep_tau, Amp), prior=priors, fcn=self.fit_func(mom_ls), maxit=10000, fitter='scipy_least_squares') # scipy_least_squares   # gsl_multifit

        else:
            fit_result = lsf.nonlinear_fit(data=(t_tsep_tau, Amp), prior=priors, fcn=self.fit_func(mom_ls), maxit=10000, fitter='scipy_least_squares', p0=best_p0)

        print(t_tsep_tau)


        if corr == True:
            corr = gv.evalcorr(Amp)
        else:  
            corr = 0      
        
        return fit_result, corr