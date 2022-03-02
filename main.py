from file_operations import cut_data_wl, load_file
from gui import Gui, set_wavelength_window
import sys

# set the limits for the gaussian parameters
depth_min, depth_max = 0.01, 1
fwhm_min, fwhm_max = 0.0001, 0.03
xaxis_step = 0.1

def main():
    # load the spectrum
    try:
        fpath, fname, data = load_file()
    except AttributeError:
        print('No file selected.')
        sys.exit(0)

    # select the wavelength window
    xmin, xmax = set_wavelength_window(data['wl'])
    if xmin == -1 or xmax == -1 or (xmax - xmin <= 0):
        print('Error: wrong or no limits selected.')
        sys.exit(1)

    # cut the spectrum to the desired wavelength window
    data = cut_data_wl(data, xmin, xmax)
    if len(data['wl']) == 0:
        print('Error: no data left after cutting.')
        sys.exit(1)

    # run the GUI
    gui = Gui(fpath, fname, data, xmin, xmax, xaxis_step, depth_min, depth_max, fwhm_min, fwhm_max)

if __name__ == "__main__":
    main()
