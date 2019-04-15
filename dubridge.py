import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from PhotoelectricEffect.utils_dos import thermal_voltage, dub_crop, log_i_err, a_to_cutoff, e
from .read_data import read_photoelectric_data
from .straightline import wavelengths, colours_plotting
from PhotoelectricEffect.utils_dataproc import wavelength_to_frequency, chi_squared, find_linear_fit, find_linear_lines, frequency_err
from PhotoelectricEffect.utils_plot import estimate_h
from matplotlib import rc


def dub_neg_phi(x, num_terms=10):
    func_neg = 0
    for j in range(num_terms):
        # print("For a negative value of x the {}'th term {} {}".format(j, round(x, 4), round(func_neg, 4)))
        func_neg += (-1)**j * (np.exp((j+1)*x) / (j+1)**2)
    return func_neg


def dub_pos_phi(x, num_terms=10):
    func_pos = (((np.pi**2) / 6) + ((x**2) / 2))
    for j in range(num_terms):
        # print("For a positive value of x The {}'th term {} {}".format(j, round(x, 4), round(func_pos, 4)))
        func_pos -= ((-1)**j * (np.exp(-(j + 1)*x) / (j + 1) ** 2))
    return func_pos


def dub_phi(x):
    return np.piecewise(float(x), [x <= 0., x > 0.], [dub_neg_phi(x, 10), dub_pos_phi(x, 10)])


def dub_plot():
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    v_cut = []
    frequencies, freq_err = [], []
    v_cut_err = []
    fig, ax = plt.subplots()
    for colour, wavelength in wavelengths.items():
        print("\nProcessing {} with wavelength {}nm".format(colour, wavelength))
        file_name = 'PhotoelectricEffect\Data\Current_Voltage\current_' + colour + '.csv'
        i, i_err, v = read_photoelectric_data(file_name)
        flip_v = [v_elem * -1 for v_elem in v]
        # fig, ax = plt.subplots()

        therm_v = thermal_voltage(flip_v)
        crop_i, crop_i_err, therm_v_crop = dub_crop(i, i_err, therm_v, low_v_lim=-20, upp_v_lim=60)
        # crop_i = normalise_i(crop_i, max(crop_i))

        ln_crop_i_err = np.asarray(log_i_err(crop_i, crop_i_err, mul_fac=10))
        ax.errorbar(x=therm_v_crop, y=np.log(crop_i), yerr=ln_crop_i_err, fmt='.', color=colours_plotting[colour], capsize=5,
                    label='Data ' + str(colour))

        a, a_err, b, b_err = dub_fit(therm_v_crop, np.log(crop_i), 40, 0)
        print("a {} +/- {}".format(a, a_err))
        print("b {} +/- {}".format(b, b_err))

        v_range = np.arange(-20, 60, 0.5)
        dub_line = dub_equ(v_range, a=a, b=b)
        ax.plot(v_range, dub_line, color='black', linestyle='--')
        ax.set_xlabel(r"$\mathbf{eV / {k_b}T}$")
        ax.set_ylabel(r"\textbf{ln}$\mathbf{(I)}$")

        v_cut.append(a_to_cutoff(a))
        v_cut_err.append(a_to_cutoff(a_err * 3))
        frequencies.append(wavelength_to_frequency(wavelength))
        freq_err.append(frequency_err(wavelength, colour))

    estimate_h(frequencies, freq_err, v_cut, v_cut_err)


def dub_fit(v, ln_i, a_est, b_est):
    [a, b], err = curve_fit(dub_equ, p0=[a_est, b_est], xdata=v, ydata=ln_i)
    a_err, b_err = np.sqrt(np.diag(err))
    return a, a_err, b, b_err


def dub_equ(v_therm, a, b):
    dubridge_theory = []
    for v_therm_elem in v_therm:
        dubridge_theory.append(b + np.log(dub_phi(a + v_therm_elem)))
    return dubridge_theory


def dubridge_fit_plot():
    v_range = np.arange(-20, 60, 0.5)
    dubridge_equation = dub_equ(v_range, a=0, b=0)
    fig, ax = plt.subplots()
    ax.plot(v_range, dubridge_equation)
    plt.show()

