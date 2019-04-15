import matplotlib.pyplot as plt
import numpy as np
from PhotoelectricEffect.utils_dataproc import find_linear_lines, find_linear_fit, chi_squared, frequency_err
from PhotoelectricEffect.utils_dos import e
from matplotlib import rc


def estimate_h(frequencies, freq_err, v_cut, v_cut_err):
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    # freq_err = [freq_err_elem for freq_err_elem in freq_err]

    fig_h, ax_h = plt.subplots()
    frequencies = np.asarray(frequencies)
    freq_err = [freq_err_elem * 3 for freq_err_elem in freq_err]
    freq_axis = np.arange(min(frequencies) - 1e15, max(frequencies) + 1e15, 1e12)
    m, m_err, c, c_err = find_linear_fit(v_cut, v_cut_err, frequencies)
    y0, y1, y2 = find_linear_lines(m, m_err, c, c_err, frequencies)

    print("\nR-chi: {}"
          "\nh: {} +/- {}"
          "\npsi: {} +/- {}".format(chi_squared(y0, v_cut, v_cut_err) / 4, m * e, m_err * e, c, c_err))

    y0, y1, y2 = find_linear_lines(m, m_err, c, c_err, freq_axis)
    ax_h.plot(freq_axis, y0, '--', color="red")
    ax_h.fill_between(freq_axis, y1, y2, color='grey', alpha=0.2)
    ax_h.set_xlabel(r"\textbf{Frequencies} (Hz)")
    ax_h.set_ylabel(r"\textbf{Cut-Off Voltage} (V)")
    ax_h.set_xlim(3e14, 10e14)
    ax_h.set_ylim(0.25, 1.25)

    ax_h.errorbar(x=frequencies, y=v_cut, yerr=v_cut_err, xerr=freq_err, fmt='.', color="black", capsize=3, label='Data')
    plt.show()

