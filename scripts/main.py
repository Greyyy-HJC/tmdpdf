"""
This is the main function of this program. It will demonstrate the whole process of TMDPDF analysis.

The main steps are:
1. g.s. fit
2. renormalization
3. extrapolation in the coordinate space and FT
4. matching
5. large Pz and physical pion mass extrapolation
6. combine the gamma t and gamma z to get the final result

"""

# %%
import multiprocessing as mp
import os
import gvar as gv
import numpy as np
from tqdm import tqdm

from libtmdpdf.compat.funcs import *
from libtmdpdf.compat.gs_fit_module import Gs_Fit
from libtmdpdf.compat.prior_setting import two_state_fit
from libtmdpdf.compat.read_raw_module import Read_Raw

# %%
#############################################################
"""
do the g.s. fit to all sets, save the log and figures, output the chi2_ls, p_value_ls, pdf_re_ls, pdf_im_ls
"""
#############################################################

#! bs fits
if False:
    from gs_fit_module import read_and_fit_bs

    # * parallel processing
    with mp.Pool() as pool:
        import warnings

        # Filter out RuntimeWarning messages
        warnings.filterwarnings("ignore", category=RuntimeWarning)

        loop_paras_ls = [
            (gamma, mass, mom, ll, b, z)
            for ll in [6]
            for gamma in ["t", "z"]
            for mass in [220, 310]
            for mom in [8, 10, 12]
            for b in range(1, 6)
            for z in range(13)
        ]

        # pool.map(read_and_fit, loop_paras_ls)

        results = list(
            tqdm(pool.imap(read_and_fit_bs, loop_paras_ls), total=len(loop_paras_ls))
        )  # * use tqdm to show the progress

#! gvar fit
if False:
    from gs_fit_module import read_and_fit_gvar

    read_and_fit_gvar()

#############################################################
#############################################################


# %%
#############################################################
"""
do the renormalization, make z dependence plots
"""
#############################################################
if True:
    # * read the gs fit result
    load = gv.load("data/read_from_here/all_after_gs_fit")

    mass = 220
    gamma = "t"
    mom = 8
    ll = 6

    re_dic = {}
    im_dic = {}

    for b in range(1, 6):
        re_dic["b{}".format(b)] = []
        im_dic["b{}".format(b)] = []
        for z in range(13):
            fit_id = "{}{}_P{}_L{}_b{}_z{}".format(mass, gamma, mom, ll, b, z)

            re_dic["b{}".format(b)].append(load[fit_id]["re"])
            im_dic["b{}".format(b)].append(load[fit_id]["im"])

    # todo renormalization

    # * make z dependence plots of different b in one figure
    label_ls = ["b{}".format(b) for b in range(1, 6)]

    errorbar_ls_plot(
        [np.arange(13) for i in range(5)],
        [gv.mean(re_dic["b{}".format(b)]) for b in range(1, 6)],
        [gv.sdev(re_dic["b{}".format(b)]) for b in range(1, 6)],
        label_ls,
        "z dependence mix b at {}{}, mom={}, L={}, real".format(mass, gamma, mom, ll),
    )

    errorbar_ls_plot(
        [np.arange(13) for i in range(5)],
        [gv.mean(im_dic["b{}".format(b)]) for b in range(1, 6)],
        [gv.sdev(im_dic["b{}".format(b)]) for b in range(1, 6)],
        label_ls,
        "z dependence mix b at {}{}, mom={}, L={}, imag".format(mass, gamma, mom, ll),
    )

#############################################################
#############################################################


# %%
