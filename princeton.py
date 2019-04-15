import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .read_data import read_photoelectric_data, moving_average
from .straightline import lin_i_lim, neg_i_lim, wavelengths, colours_plotting
from PhotoelectricEffect.utils_dos import e
from PhotoelectricEffect.utils_dataproc import wavelength_to_frequency, chi_squared, find_linear_fit, find_linear_lines, frequency_err, normalise_i
from PhotoelectricEffect.utils_plot import estimate_h
from matplotlib import rc

fermi_energy = 2.12


def dos_equation(v, a, c):
    return np.piecewise(v, [v < -1*c, v >= -1*c],
                        [0, lambda v: a * (0 - (v + c)**(3./2.))])


def dosfd_integrand(v, mu, t):
    return v**0.5 / (np.exp((v - mu)/t) + 1)


def fit_dos(v, i, a_est, c_est):
    [a, c], err = curve_fit(dos_equation, p0=[a_est, c_est], xdata=v, ydata=i)
    a_err, c_err = np.sqrt(np.diag(err))
    return a, a_err, c, c_err


def fit_dosfd(v, i_dif, mu_est, t_est):
    [mu, t], err = curve_fit(dosfd_integrand, p0=[mu_est, t_est], xdata=v, ydata=i_dif)
    mu_err, t_err = np.sqrt(np.diag(err))
    return mu, mu_err, t, t_err


def dos_crop(i, i_err, v, colour):
    up_lim = lin_i_lim[colour][0]
    low_lim = neg_i_lim[colour]
    dos_i = i[low_lim:up_lim]
    dos_i_err = i_err[low_lim:up_lim]
    dos_v = v[low_lim:up_lim]
    return dos_i, dos_i_err, dos_v


def find_fitted_dos(a, a_err, c, c_err, x):
    y0 = dos_equation(x, a, c)
    # y_err = np.sqrt((m_err * m) ** 2 + c_err ** 2)
    # y_up = m * x + (c + y_err)
    # y_low = m * x + (c - y_err)
    return y0


def dosfd():
    v_cut, v_cut_err = [], []
    frequencies = []
    # fig, ax = plt.subplots()
    for colour, wavelength in wavelengths.items():
        fig, ax = plt.subplots()
        print("\n Processing {} with wavelength {}nm".format(colour, wavelength))
        file_name = 'PhotoelectricEffect\Data\Current_Voltage\current_' + colour + '.csv'
        i, i_err, v = read_photoelectric_data(file_name)

        temporal_window = 1
        if colour == "green":
            i_smooth = i
            v_smooth = v
            i_smooth_err = i_err
        else:
            i_smooth = moving_average(i, temporal_window)
            v_smooth = [(v_elem / 10) + v[0] for v_elem in np.arange(0, len(i_smooth))]
            i_smooth_err = moving_average(i_err, temporal_window)

        i_smooth_dif = np.gradient(i_smooth)
        ax.set_title(colour)
        x_axis = np.arange(v[0], v[-1], 0.01)
        dos_i, dos_i_err, dos_v = dos_crop(i_smooth_dif, i_smooth_err, v_smooth, colour)
        ax.plot(dos_v, dos_i, '.', markersize=5, color=colours_plotting[colour], label='Data ' + str(colour))
        a_est = -1 * max(dos_i) / 10.
        c_est = -1 * min(dos_v) - 0.5
        a, a_err, c, c_err = fit_dos(dos_v, dos_i, a_est, c_est)
        print("A // Guess: {} // Calculated: {} +/- {}".format(a_est, a, a_err))
        print("C // Guess: {} // Calculated: {} +/- {}".format(c_est, c, c_err))
        x_axis = np.arange(dos_v[0], dos_v[-1], 0.01)
        y0 = find_fitted_dos(a, a_err, c, c_err, x_axis)
        y_guess = find_fitted_dos(a_est, a_err, c_est, c_err, x_axis)
        # ax.plot(x_axis, y_guess, '-.', color=colours_plotting[colour], label='Initial ' + str(colour) + ' fit')
        #
        ax.plot(x_axis, y0, '--', color=colours_plotting[colour], label='Final ' + str(colour) + ' fit')
        plt.legend(loc='best')
        v_cut.append(c - fermi_energy)
        v_cut_err.append(c_err)
        frequencies.append(wavelength_to_frequency(wavelength))
        # plt.savefig("C:\\Users\\18072\\OneDrive\\Documents\\Academic\\University\\Physics\\3rd_Year\\Labs\\PhotoelectricEffect\\Graphs\\DosFitting" + str(colour) + "_dos_fit_.png")


    m, m_err, c, c_err = find_linear_fit(v_cut, v_cut_err, frequencies)
    print(c)
    freq_axis = np.arange(min(frequencies), max(frequencies), 100000000)
    y0, y1, y2 = find_linear_lines(m, m_err, c, c_err, freq_axis)
    # ax.plot(freq_axis, y0)
    # ax.plot(freq_axis, y1)
    # ax.plot(freq_axis, y2)
    print(m)
    print(m*1.6*10**-19)
    print(m_err*1.6*10**-19)

    # ax.errorbar(x=frequencies, y=v_cut, yerr=v_cut_err, fmt='.')
    plt.show()


def dos():
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    v_cut, v_cut_err = [], []
    frequencies, freq_err = [], []
    fig, ax = plt.subplots()
    for colour, wavelength in wavelengths.items():
        # fig, ax = plt.subplots()
        print("\n Processing {} with wavelength {}nm".format(colour, wavelength))
        file_name = 'PhotoelectricEffect\Data\Current_Voltage\current_' + colour + '.csv'
        i, i_err, v = read_photoelectric_data(file_name)
        i_err = [i_err_elem * 10 for i_err_elem in i_err]
        dos_i, dos_i_err, dos_v = dos_crop(i, i_err, v, colour)
        dos_i_err = normalise_i(dos_i_err, max(dos_i))
        dos_i = normalise_i(dos_i, max(dos_i))

        # ax.plot(dos_v, dos_i, '.', markersize=5, color=colours_plotting[colour], label='Data ' + str(colour))
        a_est = -1 * max(dos_i) / 100.
        c_est = -1 * min(dos_v) - 0.5
        a, a_err, c, c_err = fit_dos(dos_v, dos_i, a_est, c_est)
        print("A // Guess: {} // Calculated: {} +/- {}".format(a_est, a, a_err))
        print("C // Guess: {} // Calculated: {} +/- {}".format(c_est, c, c_err))
        x_axis = np.arange(dos_v[0], dos_v[-1], 0.01)
        y0 = find_fitted_dos(a, a_err, c, c_err, x_axis)
        y_guess = find_fitted_dos(a_est, a_err, c_est, c_err, x_axis)
        # ax.plot(x_axis, y_guess, '-.', color=colours_plotting[colour], label='Initial ' + str(colour) + ' fit')
        if colour == "green":
            ax.errorbar(x=dos_v, y=dos_i, yerr=dos_i_err, fmt='.', marker='s', markersize=5, color=colours_plotting[colour], capsize=5,
                        label='Data ' + str(colour))
            ax.plot(x_axis, y0, '-.', color='black', label='Fit ' + str(colour) + ' fit')
        elif colour == "violet":
            ax.errorbar(x=dos_v, y=dos_i, yerr=dos_i_err, fmt='.', marker='d', markersize=5, color=colours_plotting[colour],
                        capsize=5, label='Data ' + str(colour))
            ax.plot(x_axis, y0, '--', color='black', label='Fit ' + str(colour) + ' fit')
        # else:
            # ax.errorbar(x=dos_v, y=dos_i, yerr=dos_i_err, fmt='.', marker='o', markersize=5,
            #             color=colours_plotting[colour], capsize=5, label='Data ' + str(colour))
            # ax.plot(x_axis, y0, '-', color='black', label='Fit ' + str(colour) + ' fit')
        # plt.legend(loc='best')
        v_cut.append(c - fermi_energy)
        v_cut_err.append(c_err*7)
        frequencies.append(wavelength_to_frequency(wavelength))
        freq_err.append(frequency_err(wavelength, colour))
        ax.set_xlabel(r"\textbf{Voltage} (V)")
        ax.set_ylabel(r"\textbf{Normalised Current}")

    estimate_h(frequencies, freq_err, v_cut, v_cut_err)
    print(frequencies)
    print(v_cut)
    print(v_cut_err)
