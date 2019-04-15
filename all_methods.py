from PhotoelectricEffect.straightline import straight_line_fit
from PhotoelectricEffect.princeton import dos
from PhotoelectricEffect.dubridge import dub_plot
from PhotoelectricEffect.geometry import geom_plot
import numpy as np
import matplotlib.pyplot as plt
from PhotoelectricEffect.utils_dataproc import find_linear_lines, find_linear_fit, chi_squared
from PhotoelectricEffect.utils_dos import e


def all_methods_plotted():
    fig_h, ax_h = plt.subplots()
    for method in [straight_line_fit(), dos(), dub_plot(), geom_plot()]:
        frequencies, v_cut, v_cut_err = method
        frequencies = np.asarray(frequencies)
        freq_axis = np.arange(min(frequencies) - 1e14, max(frequencies) + 1e14, 1e12)
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

        ax_h.errorbar(x=frequencies, y=v_cut, yerr=v_cut_err, fmt='.', color="black", capsize=5, label='Data')
        plt.show()
