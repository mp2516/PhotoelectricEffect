import numpy as np
from scipy.optimize import curve_fit

c = 3 * 10 ** 8

def straight_line(x, a, b):
    return x * a + b


def wavelength_to_frequency(wavelength):
    return 3 * 10**17 / wavelength


def chi_squared(fitted, observed, observed_error):
    chisquared = []
    for num in range(len(fitted)):
        chisquared.append(((observed[num] - fitted[num]) / observed_error[num])**2)
    return np.sum(chisquared)


def normalise_i(i, max_i):
    return [i_elem / max_i for i_elem in i]


def find_linear_fit(i, i_err, v):
    [m, c], err = curve_fit(straight_line, xdata=v, ydata=i,
                            sigma=i_err, absolute_sigma=True)
    m_err, c_err = np.sqrt(np.diag(err))
    return m, m_err, c, c_err


def find_linear_lines(m, m_err, c, c_err, x):
    y0 = m * x + c
    y_err = np.sqrt((m_err * m) ** 2 + c_err ** 2)
    y_up = m * x + (c + y_err)
    y_low = m * x + (c - y_err)
    return y0, y_up, y_low


def frequency_err(wavelength, colour):
    wavelength_err = {"red": 0.8, "blue": 0.5, "yellow": 1, "uva": 0.9, "green": 3.1, "violet": 0.6}
    # return np.sqrt((((c * wavelength_err[colour] * 10**9) / wavelength**2) ** 2 +
    #                 (2*((c*10**9)/wavelength))*10**-3)**2)
    return (2*((c*10**9)/wavelength))*10**-3
