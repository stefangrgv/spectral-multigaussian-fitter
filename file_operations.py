from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfile
import tkinter as tk
import numpy as np
import os


def load_file():
    '''
    Opens a file dialog and loads a spectrum

    Returns (filepath, filename, data)
    '''

    tk.Tk().withdraw()

    fpath = askopenfilename(filetypes=[('Narval files','*.s *.snr *.spect'), ('LSD files', '*.lsd *.avg')])

    chunks = fpath.split('/')[-1].split('.')
    fname = ''
    for i in range(len(chunks)-1):
        fname += chunks[i] + '.'
    fname = fname[:-1]

    wl, i = np.loadtxt(fpath, usecols=(0,1), skiprows=2, unpack=True)
    data = {'wl':wl, 'i':i}

    return (fpath, fname, data)

def cut_data_wl(data, xmin, xmax):
    '''
    Restricts the spectrum to the desired wavelength window [xmin, xmax]

    Returns the restricted data
    '''

    cut_wl = data['wl'][ (data['wl'] >= xmin) & (data['wl'] <= xmax) ]
    cut_i  = data['i'][ (data['wl'] >= xmin) & (data['wl'] <= xmax) ]

    cut_data = {'wl': cut_wl, 'i': cut_i}
    return cut_data

def save_file(fname, fit):
    '''
    Saves the fit result to an ASCII file
    '''

    f = asksaveasfile(filetypes=[('.fit files', '*.fit')])
    if f is not None:
        f.write('Fit results:\n')
        for n in range(fit['n']):
            f.write('#%i\n'%(n+1))
            f.write('Line depth\t=\t% 6.5f\n'%(fit['a'][n]))
            f.write('Line peak\t=\t% 10.6f\n'%(fit['x0'][n]))
            f.write('Line FWMH\t=\t% 7.5f\n\n'%np.abs(2*fit['sigma'][n]))
        f.write('Continuum\t=\t% 4.3f\n'%(fit['continuum']))
        f.close()
