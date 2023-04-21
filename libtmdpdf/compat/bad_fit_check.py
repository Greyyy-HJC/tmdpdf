# %%
'''
check the fits quality in the dump folder
'''

import os

import gvar as gv
import numpy as np

bad_fit_bs_id = {}
total = 0

for file in os.listdir('dump/gs_fit_bs'):
    temp = [x for x in file.split('_')[0:7]]
    gamma = temp[0][-1]
    mass = int(temp[0][:-1])
    mom = int(temp[1][1:])
    ll = int(temp[2][1:])
    b = int(temp[3][1:])
    z = int(temp[4][1:])
    tmax = int(temp[5][4:])
    tau_cut = int(temp[6][3:])
    

    if tmax == 8 and tau_cut == 1 and ll == 6: #b < 3 and z < 7 and mom == 10 and tmax == 9 and tau_cut == 0:
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


# %%
'''
move the good quality fits to the read_from_here folder
'''

import shutil

# specify the source and destination paths
src_path = 'dump/gs_fit_bs/'
dest_path = 'read_from_here/gs_fit/'

# use the shutil.copy() function to copy the file
for file in os.listdir(src_path):
    temp = [x for x in file.split('_')[0:7]]
    gamma = temp[0][-1]
    mass = int(temp[0][:-1])
    mom = int(temp[1][1:])
    ll = int(temp[2][1:])
    b = int(temp[3][1:])
    z = int(temp[4][1:])
    tmax = int(temp[5][4:])
    tau_cut = int(temp[6][3:])

    # if tmax == 8 and tau_cut == 1 and ll == 6:
    #     shutil.copy(src_path+file, dest_path+file)


# %%
'''
make sure the fits in the read_from_here folder are good
'''
# import gvar as gv

# bad_fit_bs_id = {}
# total = 0

# for file in os.listdir('read_from_here/gs_fit'):
#     temp = [x for x in file.split('_')[0:7]]
#     gamma = temp[0][-1]
#     mass = int(temp[0][:-1])
#     mom = int(temp[1][1:])
#     b = int(temp[3][1:])
#     z = int(temp[4][1:])
#     tmax = int(temp[5][4:])
#     tau_cut = int(temp[6][3:])
    

#     if True:
#         load = gv.load('read_from_here/gs_fit/'+file)['Q']
#         for idx in range(800):
#             if load[idx] < 0.05:
#                 if str(idx) not in bad_fit_bs_id:
#                     bad_fit_bs_id[str(idx)] = 0
#                 bad_fit_bs_id[str(idx)] += 1
#             total += 1

# bad_total = sum([bad_fit_bs_id[key] for key in bad_fit_bs_id])

# print(bad_total/total)
# print(total)
# %%
