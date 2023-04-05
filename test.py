#%% 
import numpy as np
import gvar as gv

# define a function to generate a random list for bootstrap sampling with max = N_conf and times = N_bs
def gen_bs_ls(N_conf, N_bs):
    bs_ls = []
    for i in range(N_bs):
        bs_ls.append(np.random.randint(N_conf, size=N_conf))
    return bs_ls

# bs_ls = gen_bs_ls(1000, 800)
# gv.dump(bs_ls, 'data_raw/bs_ls_800')
# %%
