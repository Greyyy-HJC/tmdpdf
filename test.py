#%% 
'''
generate bs list
'''

import numpy as np
import gvar as gv

# define a function to generate a random list for bootstrap sampling with max = N_conf and times = N_bs
def gen_bs_ls(N_conf, N_bs):
    bs_ls = []
    for i in range(N_bs):
        exclude_conf = np.array([192, 207, 382, 523, 775, 920])

        random_ls = np.random.randint(N_conf, size=N_conf)

        while np.in1d(random_ls,exclude_conf).any():
            random_ls = np.random.randint(N_conf, size=N_conf)

        bs_ls.append(np.random.randint(N_conf, size=N_conf))
    return bs_ls

# bs_ls = gen_bs_ls(1000, 800)
# gv.dump(bs_ls, 'data_raw/bs_ls_800_clean')


# exclude_conf = np.array([3, 6, 9, 8])

# random_ls = np.random.randint(10, size=10)

# while np.in1d(random_ls,exclude_conf).any():
#     random_ls = np.random.randint(10, size=10)

# print(random_ls)


# %%
'''
select bad fits
'''

import numpy as np
import gvar as gv
import os

bad_fit_bs_id = {}
total = 0

for file in os.listdir('dump/gs_fit_bs'):
    temp = [x for x in file.split('_')[0:5]]
    gamma = temp[0][-1]
    mass = int(temp[0][:-1])
    mom = int(temp[1][1:])
    b = int(temp[3][1:])
    z = int(temp[4][1:])
    

    if b < 3 and z < 9 and mom == 8:
        load = gv.load('dump/gs_fit_bs/'+file)['Q']
        for idx in range(800):
            if load[idx] < 0.05:
                if str(idx) not in bad_fit_bs_id:
                    bad_fit_bs_id[str(idx)] = 0
                bad_fit_bs_id[str(idx)] += 1
            total += 1

bad_total = sum([bad_fit_bs_id[key] for key in bad_fit_bs_id])

print(bad_total/total)
print(total)


# very_bad_fit_bs_id = [] # bad fit with chi2 > 2 for more than 50 times
# for key in bad_fit_bs_id:
#     if bad_fit_bs_id[key] > 50:
#         # print('>>> '+key)
#         # print(bad_fit_bs_id[key])
#         very_bad_fit_bs_id.append(key)

# bs_ls_800 = gv.load('data_raw/bs_ls_800')
# very_bad_fit_conf_id = {}
# for idx in very_bad_fit_bs_id:
#     for conf_id in bs_ls_800[int(idx)]:
#         if str(conf_id) not in very_bad_fit_conf_id:
#             very_bad_fit_conf_id[str(conf_id)] = 0
#         very_bad_fit_conf_id[str(conf_id)] += 1

# print(np.mean([very_bad_fit_conf_id[key] for key in very_bad_fit_conf_id]))


# for key in very_bad_fit_conf_id:
#     if very_bad_fit_conf_id[key] > 40:
#         print('>>> '+key)
#         print(very_bad_fit_conf_id[key])



# %%
'''
throw away the bad configs
'''

import numpy as np
from funcs import *
from read_raw_module import Read_Raw
read_raw = Read_Raw('data_raw/')

mass = 220
mom = 8

pt2 = read_raw.read_2pt(mass, mom)
meff_ls = []
for idx in range(len(pt2)):
    meff = pt2_to_meff(pt2[idx])
    meff_ls.append(meff)

print(np.shape(meff_ls))

# %%
'''
adjust the bad fits
'''
import numpy as np
import gvar as gv

from read_raw_module import Read_Raw
from gs_fit_module import Gs_Fit
from prior_setting import two_state_fit


def read_and_fit(loop_paras):
    #! here is the fitting parameter setting
    ra_tmax = 8
    tau_cut = 1

    gamma, mass, mom, ll, b, z = loop_paras
    fit_id='{}{}_P{}_L{}_b{}_z{}_tmax{}_cut{}'.format(mass, gamma, mom, ll, b, z, ra_tmax, tau_cut)

    #* if the fit result already exists, skip it
    # if os.path.exists('dump/gs_fit_bs/{}_Q_chi_re_im'.format(fit_id)):
    #     return

    read_raw = Read_Raw('data_raw/')

    data_dic = {}
    temp_2pt = read_raw.read_2pt_bs(mass, mom)
    data_dic['2pt_re'] = np.real( temp_2pt )
    data_dic['2pt_im'] = np.imag( temp_2pt )



    for tseq in range(4, 9): # this tseq range comes from the raw data
        temp_ra_re, temp_ra_im = read_raw.read_ratio_bs(gamma, mass, mom, ll, b, z, tseq=tseq)

        data_dic['ra_re_tseq_{}'.format(tseq)] = temp_ra_re[:, 1:tseq]
        data_dic['ra_im_tseq_{}'.format(tseq)] = temp_ra_im[:, 1:tseq]

    #todo
    # for key in data_dic:
    #     temp = gv.dataset.avg_data(data_dic[key])
    #     print('>>>'+key)
    #     print(temp)


    gs_fit = Gs_Fit(two_state_fit(), fit_id)
    gs_fit.para_set(pt2_tmin=3, pt2_tmax=9, ra_tmin=4, ra_tmax=ra_tmax, tau_cut=tau_cut) #! here is the fitting parameter setting

    p_value_ls, chi2_ls, pdf_re_ls, pdf_im_ls = gs_fit.main(data_dic)

    Q_chi_re_im = {}
    Q_chi_re_im['Q'] = np.array(p_value_ls)
    Q_chi_re_im['chi2'] = np.array(chi2_ls)
    Q_chi_re_im['re'] = np.array(pdf_re_ls)
    Q_chi_re_im['im'] = np.array(pdf_im_ls)

    gv.dump(Q_chi_re_im, 'dump/gs_fit_bs/{}_Q_chi_re_im'.format(fit_id))


    bad_count = len( [v for v in p_value_ls if v < 0.05] )
    print('bad_count: {}'.format(bad_count))

    return

fit_id = '220z_P10_L6_b5_z3_Q_chi_re_im'
temp = [x for x in fit_id.split('_')[0:5]]
gamma = temp[0][-1]
mass = int(temp[0][:-1])
mom = int(temp[1][1:])
ll = int(temp[2][1:])
b = int(temp[3][1:])
z = int(temp[4][1:])

read_and_fit((gamma, mass, mom, ll, b, z))

# %%
import gvar as gv
import numpy as np

test = gv.load('dump/gs_fit_bs_chi2_dic')
print([key for key in test])
print(np.shape( test['total']) )

# %%
import gvar as gv

test = gv.load('dump/gs_fit_gvar/all_after_gs_fit')
print([key for key in test])
# %%
