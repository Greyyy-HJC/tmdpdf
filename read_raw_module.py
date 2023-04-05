'''
This is the module to read from the raw data file. 

The input will be the path of the raw data folder.
The output will be a list with 1000 configurations. 

For 2pt, the parameters are pion mass, hadron momentum. 
For 3pt, the parameters are gamma structure, pion mass, hadron momentum, L, b, z and tseq.

gamma structure: 't', 'z'
pion mass: 220, 310
hadron momentum(mom): 8, 10, 12 
L(ll): 6, 8, 10

With all the parameters, you can access to the specific data.
'''

#! do bootstrap here, include correlation between all sets by using the same bootstrap seed
#! average both Ldir and bdir here
#! average zdir here

#%% 
import h5py as h5
from funcs import *

class Read_Raw():
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.bs_seed_path = self.folder_path + '/bs_ls_800'

    def get_pt2_path(self, mass):
        return self.folder_path + '/a12m130p_tmdpdf_m{}_2pt.h5'.format(mass)
    
    def get_pt3_path(self, gamma, mass):
        return self.folder_path + '/a12m130p_tmdpdf_3pt_m{}_{}.h5'.format(mass, gamma)
    
    def read_2pt(self, mass, mom):
        pt2_path = self.get_pt2_path(mass)
        myfile = h5.File(pt2_path, 'r')

        base_key = 'hadron_121050'
        mom_key = 'mom_{}'.format(mom)

        
        data_complex = myfile[base_key][mom_key][:,:,1] + 1j * myfile[base_key][mom_key][:,:,2]

        #* return a numpy array with shape (1000 conf, 16 t)
        return data_complex

    def read_3pt(self, gamma, mass, mom, ll, b, z, tseq):
        ll_idx = int(ll / 2 - 3) # ll = 6, 8, 10

        pt3_path = self.get_pt3_path(gamma, mass)
        myfile = h5.File(pt3_path, 'r')

        mom_key = 'mom_{}'.format(mom)
        tseq_key = 'tseq_{}'.format(tseq)

        avg_ls1 = ['Ldir_-1', 'Ldir_1']
        avg_ls2 = ['bdir_-1', 'bdir_1']

        temp1 = myfile[mom_key][tseq_key]

        #* average over Ldir and bdir
        temp2 = ( temp1[avg_ls1[0]][avg_ls2[0]][:] + temp1[avg_ls1[0]][avg_ls2[1]][:] + temp1[avg_ls1[1]][avg_ls2[0]][:] + temp1[avg_ls1[1]][avg_ls2[1]][:] ) / 4

        #* average over zdir
        temp3 = ( temp2[:,z,ll_idx,b,0,:,:] + temp2[:,z,ll_idx,b,1,:,:] ) / 2


        data_complex = temp3[:,:,1] + 1j * temp3[:,:,2]

        #* return a numpy array with shape (1000 conf, 16 tau)
        return data_complex
    
    def read_2pt_bs(self, mass, mom):
        return bootstrap(self.read_2pt(mass, mom), seed_path=self.bs_seed_path)

    def read_3pt_bs(self, gamma, mass, mom, ll, b, z, tseq):
        return bootstrap(self.read_3pt(gamma, mass, mom, ll, b, z, tseq), seed_path=self.bs_seed_path)

    def read_ratio_bs(self, gamma, mass, mom, ll, b, z, tseq):
        # here we have 3pt with shape (800 bs, 16 tau) and 2pt with shape (800 bs), so we need to reshape 2pt to (800 bs, 1) to make the division work

        ratio = self.read_3pt_bs(gamma, mass, mom, ll, b, z, tseq) / self.read_2pt_bs(mass, mom)[:,tseq].reshape(-1,1)

        # only return the real part, the shape is (800 bs, 16 tau)
        return np.real(ratio)
    

if __name__ == "__main__":
    read_raw = Read_Raw('data_raw/')

    # test = read_raw.read_2pt_bs(220, 8)
    test = read_raw.read_ratio_bs('z', 310, 8, 6, 1, 0, 4)

    print(np.shape(test))
    print(test[0])



# %%
