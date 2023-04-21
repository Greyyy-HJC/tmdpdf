# %%
import matplotlib.pyplot as plt
import numpy as np

from .funcs import *
from .head import *

t_label = r"$\rm{t (a) }$"
meff_label = r"$m_{eff}$"
chi2_label = r"$\chi^2 / d.o.f.$"
Q_label = r"$Q$"


def meff_plot(pt2_ls, ti, tf, fit_res, title):
    meff_ls = pt2_to_meff(pt2_ls)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    ax.errorbar(
        np.arange(len(meff_ls)),
        [val.mean for val in meff_ls],
        [val.sdev for val in meff_ls],
        color=orange,
        marker="D",
        label="meff",
        **errorb
    )

    t1_fit = np.linspace(ti, tf - 1, 100)
    t2_fit = np.linspace(ti + 1, tf, 100)

    x = {}
    x["proton"] = t1_fit

    c1_fit = fit_res.fcn(x, fit_res.p)["proton"]

    x = {}
    x["proton"] = t2_fit

    c2_fit = fit_res.fcn(x, fit_res.p)["proton"]

    meff_fit = []
    for i in range(100):
        meff = np.log(c1_fit[i] / c2_fit[i])
        meff_fit.append(meff)

    ax.fill_between(
        t1_fit,
        [v.mean + v.sdev for v in meff_fit],
        [v.mean - v.sdev for v in meff_fit],
        color=blue,
        alpha=0.4,
        label="fit",
    )

    ax.set_ylim([0, 1])
    # ax.set_xlim([-0.5, 1.5])
    ax.set_xlabel(t_label, **fs_p)
    ax.set_ylabel(meff_label, **fs_p)
    ax.set_title(title, font)
    ax.legend(loc="upper right")
    ax.tick_params(direction="in", **ls_p)
    ax.grid(linestyle=":")
    plt.savefig("fig/" + title + ".pdf", transparent=True)
    plt.show()


def fit_on_data_plot_2pt(
    x, gv_y, fit_res, key, title, log_folder, ylim=None, save=True
):
    """this is a general plot function to make effective mass plot with fit band on the data

    Args:
        x (_type_): _description_
        gv_y (_type_): _description_
        fit_res (_type_): _description_
        key (_type_): _description_
        title (_type_): _description_
        log_folder (_type_): _description_
        ylim (_type_, optional): _description_. Defaults to None.
        save (bool, optional): _description_. Defaults to True.
    """

    # * data part
    mx = x[:-1]
    gv_my = pt2_to_meff(gv_y)
    my = [v.mean for v in gv_my]
    myerr = [v.sdev for v in gv_my]

    # * fit part
    fit_x = np.linspace(mx[0], mx[-1], 100)
    input_x = {}
    input_x[key] = fit_x
    fit_y = fit_res.fcn(input_x, fit_res.p)[key]
    fit_mx = []
    fit_my = []
    fit_myerr = []

    for i in range(len(fit_x) - 1):
        val = np.log(fit_y[i]) - np.log(fit_y[i + 1])
        val = val / (fit_x[i + 1] - fit_x[i])
        fit_mx.append(fit_x[i])
        fit_my.append(val.mean)
        fit_myerr.append(val.sdev)

    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    ax.errorbar(mx, my, myerr, label="data", **errorb)
    ax.fill_between(
        fit_mx,
        [fit_my[i] + fit_myerr[i] for i in range(len(fit_my))],
        [fit_my[i] - fit_myerr[i] for i in range(len(fit_my))],
        alpha=0.4,
        label="fit",
    )
    ax.tick_params(direction="in", top="on", right="on", **ls_p)
    ax.grid(linestyle=":")
    ax.set_ylim(ylim)
    plt.title(title, font)
    plt.legend()
    if save == True:
        plt.savefig(log_folder + title + ".pdf", transparent=True)
    # plt.show()


def fit_on_data_plot_ratio(
    ra_t, ra_tau, ra_re_gv, ra_im_gv, fit_res, title, log_folder
):
    """
    This function is used to make a plot of the 3pt / 2pt ratio with fit results on the data points
    Plot both real and imag parts
    ra_t and ra_tau are just t_ls and tau_ls used for fits, ra_re_gv and ra_im_gv are the y values for fits

    Args:
        ra_t (_type_): _description_
        ra_tau (_type_): _description_
        ra_re_gv (_type_): _description_
        ra_im_gv (_type_): _description_
        fit_res (_type_): _description_
        title (_type_): _description_
        log_folder (_type_): _description_
    """

    tmin = min(ra_t)
    tmax = max(ra_t) + 1

    input_x = {}
    input_x["2pt_re"] = np.arange(3, 10)
    input_x["2pt_im"] = np.arange(3, 10)
    input_x["ra_re"] = [ra_t, ra_tau]
    input_x["ra_im"] = [ra_t, ra_tau]

    fit_re_val = fit_res.fcn(input_x, fit_res.p)["ra_re"]
    fit_im_val = fit_res.fcn(input_x, fit_res.p)["ra_im"]

    tau_dic = {}
    ra_re_dic = {}
    ra_im_dic = {}
    fit_re_dic = {}
    fit_im_dic = {}

    for idx in range(len(ra_t)):
        key = "tseq_{}".format(ra_t[idx])
        if key not in tau_dic:
            tau_dic[key] = []
            ra_re_dic[key] = []
            ra_im_dic[key] = []
            fit_re_dic[key] = []
            fit_im_dic[key] = []

        tau_dic[key].append(ra_tau[idx])
        ra_re_dic[key].append(ra_re_gv[idx])
        ra_im_dic[key].append(ra_im_gv[idx])
        fit_re_dic[key].append(fit_re_val[idx])
        fit_im_dic[key].append(fit_im_val[idx])

    # * plot real part
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for tseq in range(tmin, tmax):
        key = "tseq_{}".format(tseq)
        ax.errorbar(
            np.array(tau_dic[key]) - tseq / 2,
            [v.mean for v in ra_re_dic[key]],
            [v.sdev for v in ra_re_dic[key]],
            label="tseq = {}".format(tseq),
            color=color_ls[tseq - tmin],
            **errorb
        )

        ax.fill_between(
            np.array(tau_dic[key]) - tseq / 2,
            [v.mean + v.sdev for v in fit_re_dic[key]],
            [v.mean - v.sdev for v in fit_re_dic[key]],
            color=color_ls[tseq - tmin],
            alpha=0.4,
        )

    ax.tick_params(direction="in", top="on", right="on", **ls_p)
    ax.grid(linestyle=":")

    plt.title(title + "_real", font)
    plt.legend(ncol=3)
    plt.xlabel(r"$\tau - t/2$", font)
    plt.ylabel(r"g.s.", font)
    plt.savefig(log_folder + title + "_real.pdf", transparent=True)
    # plt.show()

    # * plot imag part
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)

    for tseq in range(tmin, tmax):
        key = "tseq_{}".format(tseq)
        ax.errorbar(
            np.array(tau_dic[key]) - tseq / 2,
            [v.mean for v in ra_im_dic[key]],
            [v.sdev for v in ra_im_dic[key]],
            label="tseq = {}".format(tseq),
            color=color_ls[tseq - tmin],
            **errorb
        )

        ax.fill_between(
            np.array(tau_dic[key]) - tseq / 2,
            [v.mean + v.sdev for v in fit_im_dic[key]],
            [v.mean - v.sdev for v in fit_im_dic[key]],
            color=color_ls[tseq - tmin],
            alpha=0.4,
        )

    ax.tick_params(direction="in", top="on", right="on", **ls_p)
    ax.grid(linestyle=":")

    plt.title(title + "_imag", font)
    plt.legend(ncol=3)
    plt.xlabel(r"$\tau - t/2$", font)
    plt.ylabel(r"g.s.", font)
    plt.savefig(log_folder + title + "_imag.pdf", transparent=True)
    # plt.show()


def hist_plot(
    dic,
    xlabel,
    title,
    xlim=None,
    ylim=None,
    accumulate=False,
    save=True,
    figsize=(6, 6),
):
    """
    define a function to make a histogram plot of a list of values.

    Args:
        dic (_type_): _description_
        xlabel (_type_): _description_
        title (_type_): _description_
        xlim (_type_, optional): _description_. Defaults to None.
        ylim (_type_, optional): _description_. Defaults to None.
        accumulate (bool, optional): _description_. Defaults to False.
        save (bool, optional): _description_. Defaults to True.
        figsize (tuple, optional): _description_. Defaults to (6, 6).
    """
    all = np.concatenate(list(dic.values()))
    all = all.ravel()

    fig = plt.figure(figsize=figsize)
    ax = plt.axes(plt_axes)

    if accumulate == False:
        # * for chi2 plot

        for key in dic:
            lis = dic[key]
            ax.hist(lis, bins=150, density=True, rwidth=1, alpha=0.3, label=key)

        ax.hist(
            all,
            bins=200,
            density=True,
            color=blue,
            rwidth=1,
            alpha=0.8,
            edgecolor="black",
            linewidth=1,
            label="Total",
        )

        check = np.percentile(all, 95)

        ax.axvline(check, color="red", linestyle="dashed", linewidth=1.5)
        ax.text(
            check + 0.05,
            0.8,
            "95% percentile: {}={:.2f}".format(chi2_label, check),
            color="red",
            fontname="Times New Roman",
            fontsize=14,
        )

    elif accumulate == True:
        # * for p value plot

        for key in dic:
            lis = dic[key]
            ax.hist(
                lis,
                bins=100,
                density=True,
                rwidth=1,
                alpha=1,
                histtype="step",
                cumulative=True,
                label=key,
            )

            print(
                ">>> {}: {:.2f}%".format(
                    key, len([v for v in lis if v < 0.05]) / len(lis) * 100
                )
            )

        ax.hist(
            all,
            bins=200,
            color=blue,
            density=True,
            rwidth=1,
            alpha=1,
            edgecolor="black",
            linewidth=1.5,
            histtype="step",
            cumulative=True,
            label="Total",
        )

        check = len([v for v in all if v < 0.05]) / len(all)

        ax.axvline(0.05, color="red", linestyle="dashed", linewidth=1.5)
        ax.text(
            0.05 + 0.05,
            0.7,
            "Proportion of Q < 0.05: {:.2f}%".format(check * 100),
            color="red",
            fontname="Times New Roman",
            fontsize=14,
        )

    ax.tick_params(direction="in", top="on", right="on", **ls_p)
    ax.legend(ncol=3, prop=small_font)
    ax.grid(linestyle=":")
    ax.set_xlabel(xlabel, font)
    ax.set_ylabel("Density", font)
    plt.title(title, font)
    plt.xlim(xlim)
    plt.ylim(ylim)
    if save == True:
        plt.savefig("fig/" + title + ".pdf", transparent=True)
    plt.show()


if __name__ == "__main__":
    import os

    import gvar as gv

    """
    make the chi2 distribution plot for bs gs fits
    """

    collect_chi2 = {"220t": [], "220z": [], "310t": [], "310z": []}
    collect_Q = {"220t": [], "220z": [], "310t": [], "310z": []}

    #! bs fits
    # #* iterate over all files in the folder
    # for file in os.listdir('dump/gs_fit_bs'):
    #     temp = [x for x in file.split('_')[0:7]]
    #     gamma = temp[0][-1]
    #     mass = int(temp[0][:-1])
    #     mom = int(temp[1][1:])
    #     ll = int(temp[2][1:])
    #     b = int(temp[3][1:])
    #     z = int(temp[4][1:])
    #     tmax = int(temp[5][4:])
    #     tau_cut = int(temp[6][3:])

    #     if tmax == 8 and tau_cut == 1 and ll == 6:
    #         collect_chi2['{}{}'.format(mass, gamma)].append(gv.load('dump/gs_fit_bs/'+file)['chi2'])
    #         collect_Q['{}{}'.format(mass, gamma)].append(gv.load('dump/gs_fit_bs/'+file)['Q'])

    #! gvar fits
    load = gv.load("read_from_here/all_after_gs_fit")
    for key in load.keys():
        temp = [x for x in key.split("_")]
        gamma = temp[0][-1]
        mass = int(temp[0][:-1])
        mom = int(temp[1][1:])
        ll = int(temp[2][1:])
        b = int(temp[3][1:])
        z = int(temp[4][1:])

        collect_chi2["{}{}".format(mass, gamma)].append(load[key]["chi2"])
        collect_Q["{}{}".format(mass, gamma)].append(load[key]["Q"])

    # flatten the list
    for key in collect_chi2.keys():
        collect_chi2[key] = np.array(collect_chi2[key]).reshape(-1)
        collect_Q[key] = np.array(collect_Q[key]).reshape(-1)

    hist_plot(
        collect_chi2,
        chi2_label,
        "chi2_distribution_hist".format(chi2_label),
        xlim=(0.2, 1.5),
        ylim=(-0.05, 5),
        save=True,
    )

    hist_plot(
        collect_Q,
        Q_label,
        "Q_distribution_hist",
        xlim=(0, 1),
        ylim=(-0.05, 1.05),
        accumulate=True,
        save=True,
    )


# %%
