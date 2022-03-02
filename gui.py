import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import tkinter as tk
import math
import sys
from gauss_fit import Multigaussian
from file_operations import save_file


listening_for = 'nothing' # contains info if we're listening for mouse clicks or not

xticks = [] # x-axis ticks
bounds_left, bounds_right = [], [] # the boundaires of the spectrum which will contain the fit
fit = {} # will contain information about the fit
plot_fit = False


def set_wavelength_window(wl):
    """
    Restricts the wavelength range

    Returns (wl_min: float, wl_max: float)
    """
    spectrum_wl_min, spectrum_wl_max = min(wl), max(wl)

    wl_min, wl_max = -1, -1

    win = tk.Tk()
    win.wm_title('Set wavelength window')

    tk.Label(master=win, text='Wavelength min').grid(row=0, column=0)
    wl_min_entry = tk.Entry(master=win)
    wl_min_entry.insert(0, '%f'%spectrum_wl_min)
    wl_min_entry.grid(row=0, column=1)

    tk.Label(master=win, text='Wavelength max').grid(row=1, column=0)
    wl_max_entry = tk.Entry(master=win)
    wl_max_entry.insert(0, '%f'%spectrum_wl_max)
    wl_max_entry.grid(row=1, column=1)

    # handle the exit
    def close_window():
        win.quit()
        win.destroy()

    def set_limits():
        nonlocal wl_min, wl_max

        choice_ok = False

        try:
            wl_min = float(wl_min_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", "Wavelength min must be a valid number!")
            choice_ok = False

        try:
            wl_max = float(wl_max_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", "Wavelength max must be a valid number!")
            choice_ok = False

        if wl_min >= spectrum_wl_min and wl_max <= spectrum_wl_max:
            choice_ok = True
        else:
            tk.messagebox.showerror("Error", "Wavelength window must fit in the range %f"%spectrum_wl_min + " -- %f"%spectrum_wl_max)

        if choice_ok:
            close_window()

    set_limit_button = tk.Button(master=win, text='Set limits', command=set_limits)
    set_limit_button.grid(row=2, column=1)

    win.protocol("WM_DELETE_WINDOW", close_window)
    
    win.mainloop()

    return wl_min, wl_max


def set_bounds(peak_min, peak_max):
    global bounds_left, bounds_right

    bounds_left, bounds_right = [], []

    for n in range(number_of_gaussians):
        bounds_left.append(depth_min)
        bounds_left.append(peak_min)
        bounds_left.append(fwhm_min)

        bounds_right.append(depth_max)
        bounds_right.append(peak_max)
        bounds_right.append(fwhm_max)

class Gui:
    '''
    The main GUI window class
    '''

    def __init__(self, fpath, fname, data, xmin, xmax, xaxis_step, depth_min_input, depth_max_input, fwhm_min_input, fwhm_max_input):
        global depth_min, depth_max, fwhm_min, fwhm_max, xticks
        self.fname, self.data, self.xmin, self.xmax, self.xaxis_step = fname, data, xmin, xmax, xaxis_step
        depth_min, depth_max, fwhm_min, fwhm_max = depth_min_input, depth_max_input, fwhm_min_input, fwhm_max_input

        self.root = tk.Tk()

        #COSMETICS
        self.root.wm_title(self.fname)
        # get screen width and height
        ws = self.root.winfo_screenwidth()
        hs = self.root.winfo_screenheight()
        w = 1300 # width for the Tk root
        h = 850 # height for the Tk root
        # calculate x and y coordinates for the Tk root window
        x = ws - w
        y = (hs/2) - (h/2)
        self.root.geometry('%dx%d+%d+%d' % (w, h, x, y))

        #CANVAS
        self.canvas_frame = tk.Frame(master=self.root)
        self.canvas_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.buttons_frame = tk.Frame(master=self.root)
        self.buttons_frame.pack()

        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.canvas_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        self.subplot = self.fig.add_subplot(111)
        self.subplot.set_xlim([self.xmin, self.xmax])
        self.subplot.set_xticks(np.arange(self.xmin, self.xmax, self.xaxis_step))

        xticks = []
        for a in np.arange(self.xmin, self.xmax, self.xaxis_step):
            xticks.append('%6.2f'%a)

        self.subplot.set_xticklabels(xticks)

        self.subplot.plot(self.data['wl'], self.data['i'], marker='.')

        #initialize variables that track data selection
        self.fit_start, self.fit_end = -1, -1

        #listen for mouse
        mouse_click = self.canvas.mpl_connect('button_press_event', self.onclick)
        mouse_release = self.canvas.mpl_connect('button_release_event', self.onrelease)

        #generate buttons and labels
        self.limit_button = tk.Button(master=self.buttons_frame, text='Select fit range', underline=11, command=self.fit_range)
        self.limit_button.grid(row=0, column=0)
        self.root.bind('<Alt-r>', self.fit_range)

        self.fit_button = tk.Button(master=self.buttons_frame, text='Fit', underline=0, command=self.do_fit)
        self.fit_button.grid(row=0, column=1)
        self.root.bind('<Alt-f>', self.do_fit)

        self.save_fit_button = tk.Button(master=self.buttons_frame, text='Save fit', underline=0, command=self.save_fit)
        self.save_fit_button.grid(row=0, column=2)
        self.root.bind('<Alt-s>', self.save_fit)

        tk.Label(master=self.buttons_frame, text='\t').grid(row=0, column=3)#cosmetics

        continuum_label = tk.Label(master=self.buttons_frame, text='Continuum level: ').grid(row=0, column=4)
        self.continuum_entry = tk.Entry(master=self.buttons_frame)
        self.continuum_entry.insert(0, '1')
        self.continuum_entry.grid(row=0, column=5)

        tk.Label(master=self.buttons_frame, text='\t').grid(row=0, column=6)#cosmetics

        ngaussians_label = tk.Label(master=self.buttons_frame, text='Number of gaussians: ').grid(row=0, column=7)
        self.number_of_gaussians_entry = tk.Entry(master=self.buttons_frame)
        self.number_of_gaussians_entry.insert(0, '1')
        self.number_of_gaussians_entry.grid(row=0, column=8)

        self.limit_button = tk.Button(master=self.buttons_frame, text='Set limits', underline=4, command=self.set_limits)
        self.limit_button.grid(row=0, column=9)
        self.root.bind('<Alt-l>', self.set_limits)

        #handle closing
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        self.root.mainloop()


    def mouse(self, event):
        '''
        Get the position of the mouse in data coordinates
        '''
        if event.xdata is not None:  # None would correspond to a mouseclick outside the canvas
            delta_wl = list(map(lambda x: math.fabs(x - event.xdata), self.data['wl']))
            mouse_index = delta_wl.index(min(delta_wl))

            return mouse_index


    def onclick(self, event):
        '''
        Handle mouse click
        '''

        self.clicked = self.mouse(event)


    def onrelease(self, event):
        '''
        Handle mouse release
        '''

        global listening_for

        self.released = self.mouse(event)

        if listening_for == 'fit_range':
            self.fit_start = self.clicked
            self.fit_end = self.released
            
            self.fit_range_data = {'wl':self.data['wl'][self.fit_start:self.fit_end], 'i':self.data['i'][self.fit_start:self.fit_end]}

            self.draw_all()


    def fit_range(self, event=None):
        '''
        Select the fit range using the mouse
        '''

        global listening_for, plot_fit

        self.fit_start, self.fit_end = -1, -1
        plot_fit = False
        listening_for = 'fit_range'


    def draw_all(self):
        '''
        Redraw the plot
        '''

        global number_of_gaussians, plot_fit

        self.subplot.cla()
        self.subplot.set_xlim([self.xmin, self.xmax])
        self.subplot.set_xticks(np.arange(self.xmin, self.xmax, self.xaxis_step))
        self.subplot.set_xticklabels(xticks)
        self.subplot.plot(self.data['wl'], self.data['i'], marker='.')

        if self.fit_start > 0 and self.fit_end > 0:
            self.subplot.plot(self.fit_range_data['wl'], self.fit_range_data['i'], marker='.', color='purple') #plot limits of fit area
        
        if plot_fit:
            for n in range(number_of_gaussians):
                self.subplot.plot(self.fit['wl'], self.fit['individual_gaussians'][n], marker='.', color='red') #plot the individual gaussians

            self.subplot.plot(self.fit['wl'], self.fit['fit'], marker='.', color='green') #plot the whole fit

        self.canvas.draw()


    def do_fit(self, event=None):
        ''' 
        Perform the gauss fit
        '''
        global number_of_gaussians, listening_for, bounds_left, bounds_right, plot_fit

        #set up fit
        try:
            number_of_gaussians = int(self.number_of_gaussians_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", 'The number of gaussians must be valid integer number!')
            return None
        if number_of_gaussians <= 0:
            tk.messagebox.showerror("Error", 'The number of gaussians must be above 0!')
            return None

        try:
            continuum = float(self.continuum_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", "The continuum level must be a valid number!")
            return None
        if continuum <= 0:
            tk.messagebox.showerror("Error", 'The continuum level must be above 0!')
            return None

        set_bounds(np.min(self.fit_range_data['wl']), np.max(self.fit_range_data['wl']))

        #fit
        self.multigaussian_fit = Multigaussian(number_of_gaussians=number_of_gaussians, continuum=continuum, wl_fit=self.fit_range_data['wl'], i_fit=self.fit_range_data['i'], bounds_left=bounds_left, bounds_right=bounds_right)
        self.fit = self.multigaussian_fit.gauss_fit()

        #set canvas
        plot_fit = True
        listening_for = 'nothing'
        self.draw_all()


    def save_fit(self, event=None):
        '''
        Save the fit result to a file using the function in file_operations
        '''

        save_file(self.fname, self.fit)


    def set_limits(self, event=None):
        '''
        Set the limits for the gaussian parameters: depth, peak position, full width at half-maximum
        '''
        global number_of_gaussians

        try:
            number_of_gaussians = int(self.number_of_gaussians_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", 'The number of gaussians must be valid integer number!')
            return None
        if number_of_gaussians <= 0:
            tk.messagebox.showerror("Error", 'The number of gaussians must be above 0!')
            return None

        try:
            self.continuum = float(self.continuum_entry.get())
        except ValueError:
            tk.messagebox.showerror("Error", "The continuum level must be a valid number!")
        if self.continuum <= 0:
            tk.messagebox.showerror("Error", 'The continuum level must be above 0!')
            return None

        try:
            peak_min, peak_max = np.min(self.fit_range_data['wl']), np.max(self.fit_range_data['wl'])
        except AttributeError:
            tk.messagebox.showerror("Error", 'You must first select which line(s) you want to fit using the "Select fit range" button!')
            return None

        Limiter_window = Limiter(self, self.fit_range_data['wl'], self.fit_range_data['i'], number_of_gaussians, self.continuum, depth_min, depth_max, peak_min, peak_max, fwhm_min, fwhm_max)

    def on_closing(self):
        '''
        Properly handle the exit
        '''
        self.root.destroy()
        sys.exit(0)


class Limiter:
    '''
    GUI for setting limits on the gaussian profiles
    '''

    def __init__(self, master_gui, wl_fit, i_fit, number_of_gaussians_input, continuum, depth_min, depth_max, peak_min, peak_max, fwhm_min, fwhm_max):
        global number_of_gaussians

        self.master_gui = master_gui

        self.limiter = tk.Tk()
        self.limiter.wm_title('Gaussian configuration')

        self.wl_fit, self.i_fit = wl_fit, i_fit
        number_of_gaussians = number_of_gaussians_input
        self.continuum = continuum

        #set the window position

        #interface        
        self.frames = []
        self.depth_min_entry, self.depth_max_entry, self.peak_min_entry, self.peak_max_entry, self.fwhm_min_entry, self.fwhm_max_entry = [], [], [], [], [], []
        for n in range(number_of_gaussians):
            self.frames.append(tk.Frame(master=self.limiter))
            self.frames[-1].pack()
            tk.Label(self.frames[-1], text='min').grid(row=0, column=1)
            tk.Label(self.frames[-1], text='max').grid(row=0, column=2)        
            tk.Label(self.frames[-1], text='Depth range:').grid(row=1, column=0)
            self.depth_min_entry.append(tk.Entry(self.frames[-1]))
            self.depth_min_entry[-1].insert(0, str(depth_min))
            self.depth_min_entry[-1].grid(row=1, column=1)
            self.depth_max_entry.append(tk.Entry(self.frames[-1]))
            self.depth_max_entry[-1].insert(0, str(depth_max))
            self.depth_max_entry[-1].grid(row=1, column=2)

            tk.Label(self.frames[-1], text='Peak range:').grid(row=2, column=0)
            self.peak_min_entry.append(tk.Entry(self.frames[-1]))
            self.peak_min_entry[-1].insert(0, str(peak_min))
            self.peak_min_entry[-1].grid(row=2, column=1)
            self.peak_max_entry.append(tk.Entry(self.frames[-1]))
            self.peak_max_entry[-1].insert(0, str(peak_max))
            self.peak_max_entry[-1].grid(row=2, column=2)

            tk.Label(self.frames[-1], text='HWHM range:').grid(row=3, column=0)
            self.fwhm_min_entry.append(tk.Entry(self.frames[-1]))
            self.fwhm_min_entry[-1].grid(row=3, column=1)
            self.fwhm_min_entry[-1].insert(0, str(fwhm_min))
            self.fwhm_max_entry.append(tk.Entry(self.frames[-1]))
            self.fwhm_max_entry[-1].grid(row=3, column=2)
            self.fwhm_max_entry[-1].insert(0, str(fwhm_max))


        self.fit_button = tk.Button(self.limiter, text='Fit', command=self.fit_with_limits)
        self.fit_button.pack()


        self.limiter.mainloop()


    def fit_with_limits(self):
        '''
        Use the limits that the user inputs in the Limiter window
        '''
    
        global plot_fit

        bounds_left, bounds_right = [], []
        for n in range(number_of_gaussians):
            #check if min < max for all the entries
            if float(self.depth_min_entry[n].get()) >= float(self.depth_max_entry[n].get()):
                tk.messagebox.showerror("Error", "The lower limit for the depth of gaussian #%i must be strictly lower than the upper one"%(n+1))
                return None
            if float(self.peak_min_entry[n].get()) >= float(self.peak_max_entry[n].get()):
                tk.messagebox.showerror("Error", "The lower limit for the peak poition of gaussian #%i must be strictly lower than the upper one"%(n+1))
                return None
            if float(self.fwhm_min_entry[n].get()) >= float(self.fwhm_max_entry[n].get()):
                tk.messagebox.showerror("Error", "The lower limit for the FWHM of gaussian #%i must be strictly lower than the upper one"%(n+1))
                return None

            #check if the limits of the fit are inside the selected window
            if float(self.peak_min_entry[n].get()) >= self.master_gui.fit_range_data['wl'][-1] or float(self.peak_max_entry[n].get()) <= self.master_gui.fit_range_data['wl'][0]:
                tk.messagebox.showerror("Error", "The limits of the fit must contain the selected wavelength range (%4.3f - %4.3f)"%(self.master_gui.fit_range_data['wl'][0], self.master_gui.fit_range_data['wl'][-1]))
                return None

            bounds_left.append(float(self.depth_min_entry[n].get()))
            bounds_left.append(float(self.peak_min_entry[n].get()))
            bounds_left.append(float(self.fwhm_min_entry[n].get()))

            bounds_right.append(float(self.depth_max_entry[n].get()))
            bounds_right.append(float(self.peak_max_entry[n].get()))
            bounds_right.append(float(self.fwhm_max_entry[n].get()))

        set_bounds(bounds_left, bounds_right)

        # perform the fit
        self.multigaussian_fit = Multigaussian(number_of_gaussians=number_of_gaussians, continuum=self.continuum, wl_fit=self.wl_fit, i_fit=self.i_fit, bounds_left=bounds_left, bounds_right=bounds_right)
        self.master_gui.fit = self.multigaussian_fit.gauss_fit()
        plot_fit = True
        self.master_gui.draw_all()
