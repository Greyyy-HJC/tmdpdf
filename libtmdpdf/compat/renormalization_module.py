'''
This module is about the renormalization process, including Wilson loop and ZO factor.
'''


# %%
import h5py as h5
import numpy as np


class Renormalization():
    def __init__(self, folder_path):
        self.wloop_path = 'data_raw/record_extro_wloop.h5'

    def readin_wloop(self, ll, b, z):
        ll = 6

        h5_file = h5.File('data/record_extro_wloop.h5', "r")['/b='+str(b)+'/cv'][:]
        
        select = 2*ll+z
        wloop = h5_file[select]
        wloop = np.mean(wloop, dtype=np.float64)

        return wloop
