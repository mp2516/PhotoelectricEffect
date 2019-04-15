import numpy as np
import matplotlib.pyplot as plt

from PhotoelectricEffect.utils_dataproc import wavelength_to_frequency, chi_squared, find_linear_fit, \
    find_linear_lines, frequency_err, normalise_i
from PhotoelectricEffect.utils_dos import e
from PhotoelectricEffect.utils_plot import estimate_h
from .read_data import read_photoelectric_data
from matplotlib import rc

wavelengths = {"red": 691, "blue": 436, "yellow": 578, "uva": 365, "green": 546, "violet": 405}
wavelength_err = {"red": 0.8, "blue": 0.5, "yellow": 1, "uva": 0.9, "green": 3.1, "violet": 0.6}
colours_plotting = {"red": "red", "blue": "blue", "yellow": "darkorange", "uva": "hotpink", "green": "green", "violet": "darkviolet"}
lin_i_lim = {"red": [23, 31], "blue": [60, 82], "yellow": [19, 31], "uva": [57, 77], "green": [49, 69],
                  "violet": [41, 66]}
neg_i_lim = {"red": 9, "blue": 19, "yellow": 8, "uva": 13, "green": 10, "violet": 13}


def find_intersection(neg_m, neg_m_err, lin_m, lin_m_err, neg_c, neg_c_err, lin_c, lin_c_err):
    x = (lin_c - neg_c) / (neg_m - lin_m)
    y = lin_m * x + lin_c
    d_m_sq = (neg_m - lin_m)**2
    d_c_sq = (neg_c - lin_c)**2
    d_m = np.abs(neg_m - lin_m)
    x_err = (1 / d_m) *\
            np.sqrt(neg_c_err**2 + lin_c_err**2 + (neg_m_err**2 * (d_c_sq/d_m_sq)) + (lin_m_err**2 * (d_c_sq/d_m_sq)))
    y_err = (1/d_m) *\
            np.sqrt(neg_c_err**2 * lin_m**2 + lin_c_err**2 * neg_m**2
                    + (neg_m_err**2 * lin_m**2 * (d_c_sq/d_m_sq)) + (lin_m_err**2 * neg_m_err**2 * (d_c_sq/d_m_sq)))
    return x, x_err, y, y_err


def curve_intersection(i, v):
    x_intercept = v[0]
    for num, i_elem in enumerate(i):
        if i_elem < 0:
            x_intercept = v[num]
        else:
            return x_intercept


def negative_fit(i, i_err, v, colour):
    """
    Fits a straight line to the constant negative region of the IV curve
    :param i:
    :param i_err:
    :param neg_v:
    :param colour:
    :return:
    """
    neg_i = i[:neg_i_lim[colour]]
    neg_i_err = i_err[:neg_i_lim[colour]]
    neg_v = v[:neg_i_lim[colour]]
    return neg_i, neg_i_err, neg_v


def linear_fit(i, i_err, v, colour):
    """
    Fits a straight line to the linear section of the IV curve.
    :param i:
    :param i_err:
    :param lin_v:
    :param colour:
    :return:
    """
    low_lim, up_lim = lin_i_lim[colour][0], lin_i_lim[colour][1]
    lin_i = i[low_lim:up_lim]
    lin_i_err = i_err[low_lim:up_lim]
    lin_v = v[low_lim:up_lim]
    return lin_i, lin_i_err, lin_v


def straight_line_fit():
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    intersection_point_mids = []
    intersection_squares = []
    frequencies, freq_err = [], []
    v_intersect_err, v_intersect = [], []
    fig, ax = plt.subplots()
    for colour, wavelength in wavelengths.items():
        print("\n Processing {} with wavelength {}nm".format(colour, wavelength))
        file_name = 'PhotoelectricEffect\Data\Current_Voltage\current_' + colour + '.csv'
        i, i_err, v = read_photoelectric_data(file_name)
        i = normalise_i(i, max(i))
        i_err = [i_err_elem * 10 for i_err_elem in i_err]
        x_axis = np.arange(v[0], v[-1], 0.1)
        neg_i, neg_i_err, neg_v = negative_fit(i, i_err, v, colour)
        lin_i, lin_i_err, lin_v = linear_fit(i, i_err, v, colour)
        neg_m, neg_m_err, neg_c, neg_c_err = find_linear_fit(neg_i, neg_i_err, neg_v)
        lin_m, lin_m_err, lin_c, lin_c_err = find_linear_fit(lin_i, lin_i_err, lin_v)
        print("Gradient {} +/- {}".format(neg_m, neg_m_err))
        print("y-intercept {} +/- {}".format(neg_c, neg_c_err))
        print("Gradient {} +/- {}".format(lin_m, lin_m_err))
        print("y-intercept {} +/- {}".format(lin_m, lin_m_err))
        neg_y0, neg_y_up, neg_y_low = find_linear_lines(neg_m, neg_m_err, neg_c, neg_c_err, x_axis)
        lin_y0, lin_y_up, lin_y_low = find_linear_lines(lin_m, lin_m_err, lin_c, lin_c_err, x_axis)
        x_mid, x_mid_err, y_mid, y_mid_err = find_intersection(
            neg_m, neg_m_err, lin_m, lin_m_err, neg_c, neg_c_err, lin_c, lin_c_err)
        # print(x_mid_err)
        # for neg_line in [neg_y_up, neg_y_low]:
        #     for lin_line in [lin_y_up, lin_y_low]:
        #         x, y = find_intersection(neg_line, lin_line)
        #         intersection_squares.append([x, y])
        # fig, ax = plt.subplots()
        # ax.scatter(x_mid, y_mid, color=colours_plotting[colour])
        ax.plot(v, i, '.', markersize=4, color=colours_plotting[colour], label=str(wavelength) + r" nm")
        print(neg_i_err)
        # ax.errorbar(x=neg_v, y=neg_i, yerr=neg_i_err, fmt='.', color=colours_plotting[colour], markersize=8, capsize=5,label='Linear ' + str(colour) + ' fit')
        # ax.errorbar(x=lin_v, y=lin_i, yerr=lin_i_err, fmt='.', color=colours_plotting[colour], markersize=8, capsize=5, label='Linear ' + str(colour) + ' fit')
        # ax.plot(lin_v, lin_i, '.', markersize=5, color=colours_plotting[colour],
        #         label='Linear ' + str(colour) + ' fit')
        # ax.plot(x_axis, neg_y0, color="black", label='Negative ' + str(colour) + ' fit',
        #         linestyle='--')
        # ax.plot(x_axis, lin_y0, color="black", label='Linear ' + str(colour) + ' fit',
        #         linestyle='-.')

        # ax.plot(neg_v, neg_i, '.', markersize=6, color=colours_plotting[colour])
        # ax.plot(x_axis, neg_y0, color=colours_plotting[colour], label='Linear ' + str(colour) + ' fit',
        #        linestyle='--')
        # m, c, r_value, p_value, std_err = linregress(neg_v, neg_i)
        # print("{}: {}".format(colour, r_value**2))
        # ax.plot(x_axis, neg_y0, color=self.colours_plotting[colour], label='Negative ' + str(colour) + ' fit',
        #        linestyle='-.')
        # ax.set_ylim(min(neg_i) - 0.5, max(neg_i) + 1)
        # ax.set_xlim(min(neg_v) - 0.05, max(neg_v) + 0.5)
        # ax.set_ylim(min(neg_i) - 50000, max(lin_i) + 200000)
        # ax.set_xlim(min(neg_v) - 0.1, max(lin_v) + 2)
        ax.set_xlabel(r"\textbf{Voltage} (V)")
        ax.set_ylabel(r"\textbf{Normalised Current}")
        plt.legend(loc='best')
        if colour in ["uva"]:
            x_mid -= 0.5
        if colour in ["red"]:
            pass
        else:
            v_intersect.append(x_mid * -1)
            v_intersect_err.append(x_mid_err*3)
            # intersection_point_mids.append([x_mid, y_mid])
            frequencies.append(wavelength_to_frequency(wavelength))
            freq_err.append(frequency_err(wavelength, colour))
    # v_intersect, i_intersect = zip(*intersection_point_mids)
    # ax.plot(v_intersect, i_intersect, '.', markersize=5, color="blue",
    #         label='Intersection Points')

    estimate_h(frequencies, freq_err, v_intersect, v_intersect_err)
    print(frequencies)
    print(v_intersect)
    print(v_intersect_err)

