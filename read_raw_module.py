'''
This is the module to read from the raw data file. 

The input will be the path of the raw data folder.
The output will be a list with N_bs samples, they are complex. 

For 2pt, the parameters are pion mass, hadron momentum. 
For 3pt, the parameters are gamma structure, pion mass, hadron momentum, L, b, z and tseq.

gamma structure: 't', 'z'
pion mass: 220, 310
hadron momentum(mom): 8, 10, 12 
L(ll): 6, 8, 10

With all the parameters, you can access to the specific data.

Example usage can be found in the end. 
'''

#! do bootstrap here, include correlation between all sets by using the same bootstrap seed
#! average both Ldir and bdir here
#! average zdir here, but it should be noted that we need to get ratio and separate real part and imag part first
#! the real and imag of gamma_z should be exchanged, because i gamma_z corresponds to gamma_t


#%% 
import h5py as h5
from funcs import *

class Read_Raw():
    def __init__(self, folder_path):
        self.folder_path = folder_path
        self.bs_seed_path = self.folder_path + '/bs_ls_800_clean' #todo

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
        #* shape = (1000, 13 z, 3 L, 9 b, 2 zdir, 16 t, 3 id_re_im)


        data_complex = temp2[:,z,ll_idx,b,:,:,1] + 1j * temp2[:,z,ll_idx,b,:,:,2]


        #* return a numpy array with shape (1000 conf, 2 zdir, 16 tau)
        return data_complex
    
    def read_ratio(self, gamma, mass, mom, ll, b, z, tseq):
        # this is used to check

        ratio = self.read_3pt(gamma, mass, mom, ll, b, z, tseq) / self.read_2pt(mass, mom)[:,tseq].reshape(-1,1,1)


        #* separate real and imag part then do the zdir average
        ra_pz_re = np.real( ratio[:,0,:] ) # zdir is postive
        ra_pz_im = np.imag( ratio[:,0,:] ) # zdir is postive
        ra_nz_re = np.real( ratio[:,1,:] ) # zdir is negative
        ra_nz_im = np.imag( ratio[:,1,:] ) # zdir is negative

        if z == 0:
            if gamma == 't':
                return ra_pz_re, ra_pz_im
            elif gamma == 'z':
                return -ra_nz_im, ra_pz_re #! exchange real and imag part here


        if gamma == 't':
            ra_re = ( ra_pz_re + ra_nz_re ) / 2 
            ra_im = ( ra_pz_im - ra_nz_im ) / 2 

        elif gamma == 'z': #! exchange real and imag part here
            ra_im = ( ra_pz_re - ra_nz_re ) / 2
            ra_re = ( - ra_pz_im - ra_nz_im ) / 2 

        # return a numpy array with shape (1000 conf, 16 tau)
        return ra_re, ra_im
    
    
    def read_2pt_bs(self, mass, mom):
        return bootstrap(self.read_2pt(mass, mom), seed_path=self.bs_seed_path)

    def read_3pt_bs(self, gamma, mass, mom, ll, b, z, tseq):
        return bootstrap(self.read_3pt(gamma, mass, mom, ll, b, z, tseq), seed_path=self.bs_seed_path)

    def read_ratio_bs(self, gamma, mass, mom, ll, b, z, tseq):
        # here we have 3pt with shape (800 bs, 2 zdir, 16 tau) and 2pt with shape (800 bs, 16 t), so we need to reshape 2pt[:,tseq] to (800 bs, 1, 1) to make the division work

        ratio = self.read_3pt_bs(gamma, mass, mom, ll, b, z, tseq) / self.read_2pt_bs(mass, mom)[:,tseq].reshape(-1,1,1)


        #* separate real and imag part then do the zdir average
        ra_pz_re = np.real( ratio[:,0,:] ) # zdir is postive
        ra_pz_im = np.imag( ratio[:,0,:] ) # zdir is postive
        ra_nz_re = np.real( ratio[:,1,:] ) # zdir is negative
        ra_nz_im = np.imag( ratio[:,1,:] ) # zdir is negative

        if z == 0:
            if gamma == 't':
                return ra_pz_re, ra_pz_im
            elif gamma == 'z':
                return -ra_nz_im, ra_pz_re #! exchange real and imag part here


        if gamma == 't':
            ra_re = ( ra_pz_re + ra_nz_re ) / 2 
            ra_im = ( ra_pz_im - ra_nz_im ) / 2 

        elif gamma == 'z': #! exchange real and imag part here
            ra_im = ( ra_pz_re - ra_nz_re ) / 2
            ra_re = ( - ra_pz_im - ra_nz_im ) / 2 

        # return a numpy array with shape (800 bs, 16 tau)
        return ra_re, ra_im
    

if __name__ == "__main__":
    read_raw = Read_Raw('data_raw/')

    test = np.real( read_raw.read_2pt_bs(220, 8) )
    test_meff = [pt2_to_meff(test[i])[1:10] for i in range(len(test))]

    m_avg = gv.dataset.avg_data(test_meff, bstrap=True)

    errorbar_ls_plot([np.arange(len(m_avg))], [gv.mean(m_avg)], [gv.sdev(m_avg)], title='meff', save=False)


    mean_re = []
    sdev_re = []
    mean_im = []
    sdev_im = []

    for tseq in range(4, 9):
        temp_re, temp_im = read_raw.read_ratio_bs(gamma='z', mass=220, mom=8, ll=6, b=0, z=0, tseq=tseq)

        avg_re = gv.dataset.avg_data(temp_re[:,1:tseq], bstrap=True)
        avg_im = gv.dataset.avg_data(temp_im[:,1:tseq], bstrap=True)

        mean_re.append(gv.mean(avg_re))
        sdev_re.append(gv.sdev(avg_re))
        mean_im.append(gv.mean(avg_im))
        sdev_im.append(gv.sdev(avg_im))


    x_ls = [np.arange(1, tseq) - tseq/2 for tseq in range(4, 9)]

    errorbar_ls_plot(x_ls, mean_re, sdev_re, title='real', save=False)

    errorbar_ls_plot(x_ls, mean_im, sdev_im, title='imag', save=False)

# %%