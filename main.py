'''
This is the main function of this program. It will demonstrate the whole process of TMDPDF analysis.

The main steps are:
1. g.s. fit
2. renormalization
3. extrapolation in the coordinate space and FT
4. matching
5. large Pz and physical pion mass extrapolation
6. combine the gamma t and gamma z to get the final result

'''

# %%
#* do the g.s. fit to all sets, save the log and figures, output the chi2_ls, p_value_ls, pdf_re_ls, pdf_im_ls
import numpy as np
import gvar as gv
import os
import multiprocessing as mp

from read_raw_module import Read_Raw
from gs_fit_module import Gs_Fit
from prior_setting import two_state_fit


def read_and_fit(loop_paras):
    gamma, mass, mom, ll, b, z = loop_paras
    fit_id='{}{}_P{}_L{}_b{}_z{}'.format(mass, gamma, mom, ll, b, z)

    # if the fit result already exists, skip it
    if os.path.exists('dump/gs_fit/{}_Q_chi_re_im'.format(fit_id)):
        return

    read_raw = Read_Raw('data_raw/')

    data_dic = {}
    temp_2pt = read_raw.read_2pt_bs(mass, mom)
    data_dic['2pt_re'] = np.real( temp_2pt )
    data_dic['2pt_im'] = np.imag( temp_2pt )


    for tseq in range(4, 9): # this tseq range comes from the raw data
        temp_ra_re, temp_ra_im = read_raw.read_ratio_bs(gamma, mass, mom, ll, b, z, tseq=tseq)

        data_dic['ra_re_tseq_{}'.format(tseq)] = temp_ra_re[:, 1:tseq]
        data_dic['ra_im_tseq_{}'.format(tseq)] = temp_ra_im[:, 1:tseq]

    gs_fit = Gs_Fit(two_state_fit(), fit_id)
    gs_fit.para_set(pt2_tmin=3, pt2_tmax=10, ra_tmin=4, ra_tmax=8, tau_cut=0)

    p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls = gs_fit.main(data_dic)

    Q_chi_re_im = {}
    Q_chi_re_im['Q'] = np.array(p_value_ls)
    Q_chi_re_im['chi2'] = np.array(chi2_ls)
    Q_chi_re_im['re'] = np.array(pdf_re_ls)
    Q_chi_re_im['im'] = np.array(pdf_im_ls)

    gv.dump(Q_chi_re_im, 'dump/gs_fit/{}_Q_chi_re_im'.format(fit_id))

    return

#* parallel processing
with mp.Pool() as pool: 
    mass = '310'
    ll = 6
    loop_paras_ls = [(gamma, mass, mom, ll, b, z) for gamma in ['t', 'z'] for mom in [8, 10, 12] for b in range(1, 6) for z in range(13)]

    pool.map(read_and_fit, loop_paras_ls)







# %%
'''
Select out the bad fits
'''

# bad_fit_bs_id = {}

# for file in os.listdir('dump/gs_fit'):
#     temp = [x for x in file.split('_')[0:5]]
#     b = int(temp[3][1:])
#     z = int(temp[4][1:])

#     if b < 3 and z < 13:
#         load = gv.load('dump/gs_fit/'+file)['chi2']
#         for idx in range(800):
#             if load[idx] > 2:
#                 if str(idx) not in bad_fit_bs_id:
#                     bad_fit_bs_id[str(idx)] = 0
#                 bad_fit_bs_id[str(idx)] += 1

# print(sum([bad_fit_bs_id[key] for key in bad_fit_bs_id]))


# %%
'''
check z dependence plot
'''

# import gvar as gv
# import numpy as np

# #* read the gs fit result
# re_ls = []
# im_ls = []

# for z in range(13):
#     test = gv.load('dump/gs_fit/310t_P8_L6_b2_z{}_Q_chi_re_im'.format(z))

#     re = gv.dataset.avg_data(test['re'], bstrap=True)
#     im = gv.dataset.avg_data(test['im'], bstrap=True)

#     re_ls.append(re)
#     im_ls.append(im)

# from funcs import errorbar_plot
# errorbar_plot(np.arange(13), [v.mean for v in re_ls], [v.sdev for v in re_ls], 'z dependence at b=2, mom=8, real')

# errorbar_plot(np.arange(13), [v.mean for v in im_ls], [v.sdev for v in im_ls], 'z dependence at b=2, mom=8, imag')

# %%
