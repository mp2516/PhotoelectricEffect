import matplotlib.pyplot as plt
from .read_data import read_spectrum
from.straightline import colours_plotting
from matplotlib import rc

colour_index = {"red": [2217, 2417], "blue": [1078, 1278], "yellow": [1720, 1920], "uva": [748, 948], "green": [1578, 1778], "violet": [935, 1135]}
colour_time_filter = {"red": 6.5, "blue": 6.5, "yellow": 6.5, "uva": 12.5, "green": 6.5, "violet": 12.5}
colour_time_source = {"red": 45.5, "blue": 2.5, "yellow": 1, "uva": 20.5, "green": 0.5, "violet": 7.5}
source_to_filter = {"red": "R", "blue": "VB", "yellow": "Y", "uva": "B", "green": "G", "violet": "VA"}


def filter_source_spectrum():
    rc('font', **{'family': 'serif', 'weight': 'bold', 'size': 16})
    rc('text', usetex=True)
    for colour, indices in colour_index.items():
        file_name_filter = "C:\\Users\\18072\PycharmProjects\\3rdYearLab\\PhotoelectricEffect\\Data\\Filter\\"\
                           + source_to_filter[colour] + '_Filter_' + str(colour_time_filter[colour]) + 'ms.csv'
        file_name_source = "C:\\Users\\18072\PycharmProjects\\3rdYearLab\\PhotoelectricEffect\\Data\\Source\\"\
                           + colour + '_' + str(colour_time_source[colour]) + 'ms.csv'
        f_nm, f_inten = read_spectrum(file_name_filter, colour_index[colour][0], colour_index[colour][1])
        s_nm, s_inten = read_spectrum(file_name_source, colour_index[colour][0], colour_index[colour][1])

        fig, ax = plt.subplots()
        f_inten = [f / max(f_inten) for f in f_inten]
        s_inten = [s / max(s_inten) for s in s_inten]
        ax.plot(f_nm, f_inten, '--', color='black', linewidth=2)
        ax.plot(s_nm, s_inten, color=colours_plotting[colour], linewidth=2)
        ax.set_xlabel(r"\textbf{Wavelength} (nm)")
        ax.set_ylabel(r"\textbf{Normalised Intensity}")
        # plt.savefig("C:\\Users\\18072\\OneDrive\\Documents\\Academic\\University\\Physics\\3rd_Year\\Labs\\PhotoelectricEffect\\Graphs\\Spectrum\\" + str(colour) + "_spectrum.png")
    plt.show()

