import gvar as gv


def two_state_fit():
    priors = gv.BufferDict()
    priors['E0'] = gv.gvar(1, 10)
    priors['log(dE1)'] = gv.gvar(0, 10)

    priors['pdf_re'] = gv.gvar(1, 10)
    priors['pdf_im'] = gv.gvar(1, 10)

    priors['re_c0'] = gv.gvar(1, 10)
    priors['re_c1'] = gv.gvar(1, 10)
    priors['re_c2'] = gv.gvar(1, 10)
    priors['re_c3'] = gv.gvar(1, 10)
    priors['im_c0'] = gv.gvar(1, 10)
    priors['im_c1'] = gv.gvar(1, 10)
    priors['im_c2'] = gv.gvar(1, 10)
    priors['im_c3'] = gv.gvar(1, 10)

    return priors