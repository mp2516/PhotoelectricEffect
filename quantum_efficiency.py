import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from .read_data import read_qe
from .straightline import colours_plotting, wavelengths, wavelength_err
from PhotoelectricEffect.utils_dataproc import chi_squared
from matplotlib import rc


def cal_qe(wavelength, gradient):
    h = 6.63 * 10 ** -34
    e = 1.6 * 10 ** -19
    c = 3 * 10 ** 8
    return (h * c) / (e * wavelength * gradient * 10 ** -3)


def qe_all():
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    quantum_efficiency, quantum_efficiency_err = [], []
    for colour, wavelength in wavelengths.items():
        file_name = "C:\\Users\\18072\PycharmProjects\\3rdYearLab\\PhotoelectricEffect\\Data\\Quantum_Efficiency\\QE_" + str(colour) + '.csv'
        i, i_err, pow, pow_err = read_qe(file_name)
        pow_err = [err / 2 for err in pow_err]
        # i_err = [i_err_elem / max(i) for i_err_elem in i_err]
        # i = [i_elem / max(i) for i_elem in i]
        fig, ax = plt.subplots()
        ax.errorbar(i, pow, xerr=i_err, yerr=pow_err, fmt='.', marker='o', elinewidth=1, capsize=5,
                    color='black', markersize=5)
        p, pcov = curve_fit(lambda t, a, b: a * t + b, i, pow, sigma=pow_err,
                            absolute_sigma=True)
        straight_line = [value * p[0] + p[1] for value in i]
        chisquare = chi_squared(straight_line, pow, pow_err)
        reduced_chisquare = chisquare / (len(pow) - 2)
        perr = np.sqrt(np.diag(pcov))
        xaxis = np.arange(0, max(i), 1)
        y_error = np.sqrt((perr[0] * p[0]) ** 2 + perr[1] ** 2)
        y1 = p[0] * xaxis + (p[1] + y_error)
        y2 = p[0] * xaxis + (p[1] - y_error)
        # ax.plot(xaxis, xaxis, color='green')
        ax.plot(xaxis, p[0] * xaxis + p[1], linestyle='--', color='red')
        # plt.legend(["\nGradient: {} +/- {}"
        #             "\nOffset: {} +/- {}"
        #             "\nReduced $\chi^2$: {}".format(round(p[0], 6), round(perr[0], 6), round(p[1], 5),
        #                                             round(perr[1], 5), round(reduced_chisquare, 3))], frameon=False,
        #            handlelength=0)
        ax.fill_between(xaxis, y1, y2, color='grey', alpha=0.2)
        ax.set_ylabel(r"\textbf{Current} (pA)")
        ax.set_xlabel(r"\textbf{Power} ("+r"$\mu$"+"W)")
        # ax.errorbar(pow, i, xerr=pow_err, yerr=i_err, fmt='.', marker='o', elinewidth=1, capsize=5, color='black',
        #             markersize=5)
        # p, pcov = curve_fit(lambda t, a, b: a * t + b, pow, i, sigma=i_err, absolute_sigma=True)
        # straight_line = [value * p[0] + p[1] for value in pow]
        # chisquare = chi_squared(straight_line, i, i_err)
        # reduced_chisquare = chisquare / (len(pow) - 2)
        # perr = np.sqrt(np.diag(pcov))
        # xaxis = np.arange(0, max(pow), 1)
        # y_error = np.sqrt((perr[0] * p[0]) ** 2 + perr[1] ** 2)
        # y1 = p[0] * xaxis + (p[1] + y_error)
        # y2 = p[0] * xaxis + (p[1] - y_error)
        # # ax.plot(xaxis, xaxis, color='green')
        # ax.plot(xaxis, p[0] * xaxis + p[1], linestyle='--', color='red')
        # # plt.legend(["\nGradient: {} +/- {}"
        # #             "\nOffset: {} +/- {}"
        # #             "\nReduced $\chi^2$: {}".format(round(p[0], 6), round(perr[0], 6), round(p[1], 5),
        # #                                             round(perr[1], 5), round(reduced_chisquare, 3))], frameon=False,
        # #            handlelength=0)
        # ax.fill_between(xaxis, y1, y2, color='grey', alpha=0.2)
        # ax.set_ylabel(r"\textbf{Current} (pA)")
        # ax.set_xlabel(r"\textbf{Power} (" + r"$\mu$" + "W)")
        quantum_efficiency.append(cal_qe(wavelength, p[0]))
        quantum_efficiency_err.append(0.002)
    fig_qe, ax_qe = plt.subplots()
    ax_qe.errorbar(wavelengths.values(), quantum_efficiency, xerr=wavelength_err.values(), yerr=quantum_efficiency_err, fmt='.', marker='o', elinewidth=1, capsize=5, color='black',
                markersize=5)
    wls = list(wavelengths.values())
    ax_qe.set_ylabel(r"\textbf{Quantum Efficiency}")
    ax_qe.set_xlabel(r"\textbf{Wavelength} (nm)")
    p, pcov = curve_fit(lambda t, a, b: a * t + b, wls, quantum_efficiency, sigma=np.asarray(quantum_efficiency_err), absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    xaxis = np.arange(100, 1000, 1)
    y_error = np.sqrt((perr[0] * p[0]) ** 2 + perr[1] ** 2)
    y1 = p[0] / xaxis + (p[1] + y_error)
    y2 = p[0] / xaxis + (p[1] - y_error)
    ax_qe.fill_between(xaxis, y1, y2, color='grey', alpha=0.2)
    ax_qe.plot(xaxis, p[0] * xaxis + p[1], linestyle='--', color='red')
    ax_qe.plot(wavelengths.values(), quantum_efficiency, '.')
    ax_qe.set_xlim(300, 700)
    ax_qe.set_ylim(-0.005, 0.02)
    print(p[0])
    print(perr[0])
    plt.show()