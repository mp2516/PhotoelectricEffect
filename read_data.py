import csv
import numpy as np


def moving_average(data, temporal_window=2):
    """
    Produces an array of constant probability which is normalised over the whole data set.
    Mode ‘valid’ returns output of length max(M, N) - min(M, N) + 1. The convolution product is only given for points
    where the signals overlap completely. Values outside the signal boundary have no effect.
    :param data: The list of heights to average
    :param temporal_window: The size of the transient
    :return: The convolution of the data and the temporal_window
    """
    window = np.ones(temporal_window) / temporal_window
    return np.convolve(data, window, 'valid')


def read_photoelectric_data(file_name):
    """
    Opens the .csv file and parses the data converting into lists of data:
        current, voltage and current_error
    :param file_name:
    :return: voltage, current, current_error
    """
    line_count = 0
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        voltage = []
        current, current_error = [], []
        for row in csv_reader:
            if line_count == 0:
                #print(f'Column names are {", ".join(row)}')
                line_count += 1
            elif line_count == 1:
                #print(f'\t Order of magnitude {row[0]} {row[1]} {row[2]}.')
                line_count += 1
            else:
                #print(f'\tVoltage: {row[0]}, Current: {row[1]} +/- {row[2]}.')
                voltage.append(float(row[0]))
                current.append(float(row[1]))
                current_error.append(float(row[2]))
                line_count += 1
        print(f'Processed {line_count} lines.')

    current.sort()
    current_error.sort()
    voltage.sort()
    return current, current_error, voltage


def read_spectrum(file_name, low_index, high_index):
    line_count = 0
    with open(file_name) as spectrum_file:
        spectrum_reader = csv.reader(spectrum_file, delimiter=',')
        wavelength = []
        intensity = []
        for row in spectrum_reader:
            if line_count < low_index:
                line_count += 1
            elif line_count > high_index:
                line_count += 1
            else:
                wavelength.append(float(row[0]))
                intensity.append(float(row[1]))
                line_count += 1
    return wavelength, intensity


def read_qe(file_name):
    line_count = 0
    with open(file_name) as qe_file:
        qe_reader = csv.reader(qe_file, delimiter=',')
        i, i_err = [], []
        p, p_err = [], []
        for row in qe_reader:
            if line_count in [0, 1]:
                line_count += 1
            else:
                i.append(float(row[0]))
                i_err.append(float(row[1]))
                p.append(float(row[2]))
                p_err.append(float(row[3]))
                line_count += 1
    return i, i_err, p, p_err


def normalise_photoelectric(i, i_err):
    """
    Normalises the photoelectric current by dividing by the maximum value of current in the list. Returns a list of the
    same length with the largest value in the list = 1.
    :param i: The raw current data
    :param i_err: The raw current error
    :return: norm_i, norm_i_error
    """
    max_i = max(i)
    norm_i = [i / max_i for i in i]
    norm_i_err = [i_e / max_i for i_e in i_err]
    return norm_i, norm_i_err
