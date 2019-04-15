import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .read_data import read_photoelectric_data
from .straightline import wavelengths, colours_plotting
from .utils_dataproc import straight_line, wavelength_to_frequency, normalise_i, find_linear_fit, \
    find_linear_lines, chi_squared, frequency_err
from PhotoelectricEffect.utils_dos import thermal_voltage, dub_crop, log_i_err, a_to_cutoff, e
from PhotoelectricEffect.utils_plot import estimate_h
from matplotlib import rc


def geom_neg_phi(x, a, num_terms=10):
    func_neg = (((np.pi**2) / 6) - ((x**2 - a**2) / 2) + x * np.log(1 + np.exp(x - a)))
    for j in range(num_terms):
        # print("For a negative value of x the {}'th term {} {}, a is {}".format(j, round(x, 4), round(func_neg, 4), a))
        func_neg -= (-1)**j * (np.exp((j+1)*(x - a)) / (j+1)**2)
    return func_neg


def geom_pos_phi(x, a, num_terms=10):
    func_pos = ((-x * (x - a)) + x * np.log(1 + np.exp(x - a)))
    for j in range(num_terms):
        # print("For a positive value of x The {}'th term {} {}, a is {}".format(j, round(x, 4), round(func_pos, 4), a))
        func_pos += ((-1)**j * (np.exp(-(j + 1)*(x - a)) / (j + 1) ** 2))
    return func_pos


def geom_phi(x, a):
    return float(np.select([x <= a, x > a], [geom_neg_phi(x, a, 100), geom_pos_phi(x, a, 100)]))


def geom_plot():
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    v_cut, v_cut_err = [], []
    frequencies, freq_err = [], []
    fig, ax = plt.subplots()
    for colour, wavelength in wavelengths.items():
        print("\nProcessing {} with wavelength {}nm".format(colour, wavelength))
        file_name = 'PhotoelectricEffect\Data\Current_Voltage\current_' + colour + '.csv'
        i, i_err, v = read_photoelectric_data(file_name)
        flip_v = [v_elem * -1 for v_elem in v]
        therm_v = thermal_voltage(flip_v)
        crop_i, crop_i_err, therm_v_crop = dub_crop(i, i_err, therm_v, low_v_lim=0, upp_v_lim=20)

        a, a_err, b, b_err = geom_fit(therm_v_crop, np.log(crop_i), a_est=60, b_est=0)
        print("a {} +/- {}".format(a, a_err))
        print("b {} +/- {}".format(b, b_err))

        v_cut.append(a_to_cutoff(a))
        v_cut_err.append(a_to_cutoff(a_err))
        frequencies.append(wavelength_to_frequency(wavelength))
        freq_err.append(frequency_err(wavelength, colour))

        if colour == "green":
            ln_crop_i_err = np.asarray(log_i_err(crop_i, crop_i_err, mul_fac=10))
            ax.errorbar(x=therm_v_crop, y=np.log(crop_i),
                        yerr=ln_crop_i_err, fmt='.', marker='d',
                        color=colours_plotting[colour],
                        capsize=5, label='Data ' + str(colour))
            # v_range = np.arange(0, 30, 4.3)
            # v_range = np.flip(v_range)
            # therm_v_crop = np.flip(therm_v_crop)

            v_range = np.arange(0, 20, 0.1)
            geom_line = geom_equ(v_range, a=a, b=b)
            ax.plot(v_range, geom_line, '-.', color='black', label='Fit ' + str(colour))
            ax.set_xlabel(r"$\mathbf{eV / {k_b}T}$")
            ax.set_ylabel(r"\textbf{ln}$\mathbf{(I)}$")
        else:
            ln_crop_i_err = np.asarray(log_i_err(crop_i, crop_i_err, mul_fac=10))
            ax.errorbar(x=therm_v_crop, y=np.log(crop_i), yerr=ln_crop_i_err, fmt='.', marker='o', color=colours_plotting[colour],
                        capsize=5, label='Data ' + str(colour))
            # v_range = np.arange(0, 30, 4.3)
            # v_range = np.flip(v_range)
            # therm_v_crop = np.flip(therm_v_crop)

            v_range = np.arange(0, 20, 0.1)
            geom_line = geom_equ(v_range, a=a, b=b)
            ax.plot(v_range, geom_line, '--', color='black', label='Fit ' + str(colour))
            ax.set_xlabel(r"$\mathbf{eV / {k_b}T}$")
            ax.set_ylabel(r"\textbf{ln}$\mathbf{(I)}$")
    # plt.show()

    estimate_h(frequencies, freq_err, v_cut, v_cut_err)


def geom_fit(v, ln_i, a_est, b_est):
    [a, b], err = curve_fit(geom_equ, p0=[a_est, b_est], xdata=v, ydata=ln_i)
    a_err, b_err = np.sqrt(np.diag(err))
    return a, a_err, b, b_err


def geom_equ(v_therm, a, b):
    dubridge_theory = []
    v_therm = np.flip(v_therm)
    for v_therm_elem in v_therm:
        dubridge_theory.append(b + np.log(geom_phi(v_therm_elem, a)))
    return dubridge_theory


def dubridge_fit_plot():
    v_range = np.arange(0, 30, 0.5)
    # v_range = np.flip(v_range)
    # dubridge_equation = dubridge_equ(v_range, a=0, b=0)
    fig, ax = plt.subplots()
    # ax.plot(v_range, dubridge_equation)
    # print(phi(v_range, a=30))
    ax.plot(np.flip(v_range), np.log(geom_phi(v_range, a=40)))
    ax.plot(np.flip(v_range), np.log(geom_phi(v_range, a=35)))
    ax.plot(np.flip(v_range), np.log(geom_phi(v_range, a=30)))
    plt.show()
