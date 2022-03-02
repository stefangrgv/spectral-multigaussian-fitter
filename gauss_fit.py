import numpy as np
from scipy.optimize import curve_fit


class Multigaussian:
    '''
    This class contains the parameters of the multigaussians, along with the fitting function
    '''

    def __init__(self, number_of_gaussians, continuum, wl_fit, i_fit, bounds_left, bounds_right):
        self.number_of_gaussians = number_of_gaussians
        self.continuum = continuum
        self.wl_fit = wl_fit
        self.i_fit = i_fit

        self.init_guess = []

        self.bounds = (bounds_left, bounds_right)

        for n in range(number_of_gaussians):
            self.init_guess.append((bounds_left[3*n]+bounds_right[3*n])/2)
            self.init_guess.append((bounds_left[3*n+1]+bounds_right[3*n+1])/2)
            self.init_guess.append((bounds_left[3*n+2]+bounds_right[3*n+2])/2)

    def set_initial(self, guess):
        '''
        Set the initial guess values for the parameters in the fitting process
        '''
        for n in range(self.number_of_gaussians):
            self.init_guess.append(guess[n][0])
            self.init_guess.append(guess[n][1])
            self.init_guess.append(guess[n][2])


    def gauss(self, x, a, x0, sigma):
        return -abs(a)*np.exp(-(x-x0)**2/(2*sigma**2)) #there's no continuum here


    def multigauss(self, x, *args):
        '''
        The sum of gaussians
        '''
        mg = self.continuum
        values = np.array(args)

        if values.shape == (1, self.number_of_gaussians*3):
            values = values.reshape(self.number_of_gaussians*3)
        valcount = self.number_of_gaussians*3

        for valnum in range(valcount):
            if valnum%3 == 0:
                mg += self.gauss(self.wl_fit, values[valnum], values[valnum+1], values[valnum+2])

        return mg


    def gauss_fit(self):
        self.popt, self.pcov = curve_fit(self.multigauss, self.wl_fit, self.i_fit, p0=self.init_guess, bounds=self.bounds)

        self.gaussians, self.parameters = [], []
        self.fit_a, self.fit_x0, self.fit_sigma = [], [], []
        for n in range(self.number_of_gaussians):
            self.fit_a.append(self.popt[3*n])
            self.fit_x0.append(self.popt[3*n+1])
            self.fit_sigma.append(self.popt[3*n+2])

            self.gaussians.append(self.continuum + self.gauss(self.wl_fit, self.fit_a[-1], self.fit_x0[-1], self.fit_sigma[-1]))

            self.parameters.append(self.popt[3*n])
            self.parameters.append(self.popt[3*n+1])
            self.parameters.append(self.popt[3*n+2])
            
            print('Gaussian #' + str(n+1))
            print('a = %2.2f'%self.fit_a[-1])
            print('x0= %2.2f'%self.fit_x0[-1])
            print('s = %2.2f'%self.fit_sigma[-1])

        print()

        return({'a':self.fit_a, 'x0':self.fit_x0, 'sigma':self.fit_sigma, 'wl':self.wl_fit, 'n':self.number_of_gaussians, 'individual_gaussians':self.gaussians, 'fit':self.multigauss(self.wl_fit, self.parameters), 'continuum':self.continuum, 'covariance':self.pcov})



