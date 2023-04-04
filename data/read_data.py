# %%
import sys
import os
import h5py as h5
sys.path.append(os.path.abspath("../"))

from funcs import *

#!# A3: np.imag( UU - DD )
#!# V4: np.real( UU - DD )

#* hadron: 'proton', 'proton_np'
#* flavor: 'UU', 'DD'
#* current: 'A3', 'V4'

def find_key(dict, key_word):
    for key in dict:
        if key_word in key:
            return key

def find_data_2pt(ens, hadron_ls):
    #* 2pt file name 
    fname = 'callat_test.h5'

    basekey = None
    myfile = h5.File(fname, 'r')[ens]

    data = {}
    for hadron in hadron_ls:
        data[hadron] = myfile[hadron]

    return data


def find_data_3pt(flavor, current, tsep):
    #* 3pt file name 
    file_name = 'formfac_4D_a12m130_a_proton_{}_{}_cfgs_300-5295_srcs_0-31_fft_n6.h5'

    basekey = 'gf1p0_w3p0_n30_M51p2_L520_a3p0/formfac_4D/ml0p00195' 
    basekey_2 = 'momentum_current'

    fname = file_name.format(flavor, current)
    myfile = h5.File(fname, 'r')[basekey]

    tsep_key = 'proton_{}_tsep_{}_sink_mom_px0_py0_pz0'.format(flavor, tsep)

    data = myfile[tsep_key][current][basekey_2]

    return data

def test():
    #* take a look at the C2pt shape for check

    ens = 'a09m310'
    hadron_ls = ['proton']
    test_2pt = find_data_2pt(ens, hadron_ls)['proton']

    print(np.shape(test_2pt))

    series1 = test_2pt[:,:,0,0] # symmetric, so same source and sink overlap
    series2 = test_2pt[:,:,1,0] # different sink but same source as the (0,0)

    gv1 = gv.dataset.avg_data(series1)
    gv2 = gv.dataset.avg_data(series2)

    errorbar_plot(np.arange(len(gv1)), [v.mean for v in gv1], [v.sdev for v in gv1], title='series 1 test', ylim=None, save=False)

    errorbar_plot(np.arange(len(gv2)), [v.mean for v in gv2], [v.sdev for v in gv2], title='series 2 test', ylim=None, save=False)

    plt.show()

def main():
    ens = 'a09m310'
    hadron_ls = ['proton', 'piplus']

    #!# 2pt
    temp = find_data_2pt(ens, hadron_ls)

    #* average sym sink-source and non-sys sink-source, but actually they can be fitted separately
    p_data = (temp['proton'][:,:,0,0] + temp['proton'][:,:,1,0]) / 2
    pi_data = (temp['piplus'][:,:,0,0] + temp['piplus'][:,:,1,0]) / 2

    #* take the real part of 2pt data
    p_data = np.real(p_data)
    pi_data = np.real(pi_data)


    #!# collect data
    data_set = {}
    key_ls = ['proton', 'piplus'] # indicate the key list of the out put file if necessary
    data_set['proton'] = bootstrap(p_data, times=500)
    data_set['piplus'] = bootstrap(pi_data, times=500)

    data_set_avg = gv.dataset.avg_data(data_set, bstrap=True)

    errorbar_plot(np.arange(len(data_set_avg['proton'])), [v.mean for v in data_set_avg['proton']], [v.sdev for v in data_set_avg['proton']], title='C2pt for proton', ylim=None, save=False)
    plt.show()

    return data_set_avg

#* do the test
# test()

data_set_tidy = main()
print(data_set_tidy)

#!# ATTENTION: dump will overwrite the file with the same name
# gv.dump(data_set_tidy, '../dump/data_set_tidy')


# %%
