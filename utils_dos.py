import numpy as np


def thermal_voltage(flip_v, temp=293):
    return [((-1 * flip_v_elem * e) / (k * temp)) for flip_v_elem in flip_v]


def dub_crop(i, i_err, therm_v, low_v_lim, upp_v_lim):
    dub_i, dub_i_err, dub_v = [], [], []
    for v_ind, v_elem in enumerate(therm_v):
        if low_v_lim < v_elem < upp_v_lim:
            dub_i.append(i[v_ind])
            dub_i_err.append(i_err[v_ind])
            dub_v.append(v_elem)
    return dub_i, dub_i_err, dub_v


def log_i_err(i, i_err, mul_fac):
    log_i_err = []
    for num, i_err_elem in enumerate(i_err):
        log_i_err.append(np.log(1 + (i_err_elem * mul_fac / i[num])))
    return log_i_err


def a_to_cutoff(a, temp=293):
    return (k*temp*a) / e


e = 1.6 * 10 ** -19
k = 1.38 * 10 ** -23
