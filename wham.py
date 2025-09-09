# Author: Julius Neumann <julius-neumann@gmx.net>
# 2025-09-09
###############################################

import numpy as np
import os

###############################################
# Class to calculate PMF (Potential of Mean Force) using WHAM (Weighted Histogram Analysis Method)
# Based on Roux, B. (1995). The calculation of the potential of mean force using computer simulations.
# Computer Physics Communications, 91(1), 275–282. https://doi.org/10.1016/0010-4655(95)00053-I

kB = 0.001987191 # Boltzmann constant in kcal/(mol*K)
T_default = 310 # default temperature in Kelvin (approximate body temp in humans)

class WHAM:
    def __init__(self, windows_file, kT=kB*T_default, print_progress=False):
        self.windows_file = windows_file
        self.kT = kT
        
        if print_progress:
            print('Loading index file `'+windows_file+'`')
        self.windows = np.loadtxt(windows_file, dtype=str) # list of file names for each window
        
        self.window_centers = np.array(self.windows[:,1],dtype=float) # window centers (A)
        self.k_spring = np.array(self.windows[:,4],dtype=float) # spring constant for each window (kcal/mol/A^2)
        
        self.Nw = self.window_centers.size # number of windows
        if print_progress:
            print('Found '+str(self.Nw)+' windows in simulation index file `'+windows_file+'`')

        self.xi_data = [] # list of lists of xi values over time for each window
        self.xi_data_selected = [] # selection of xi values to be analyzed

        # time values for loaded data
        self.arr_t_start = []
        self.arr_t_end = []
        self.arr_delta_t = []

        # bins on which data will be histogrammed
        self.xi_binEdges = [] # shape (n_bins+1)
        self.xi_binCenters = [] # shape (n_bins)

        # window restraining potentials, shape (Nw, n_bins)
        self.w_i = []

        # step 1: directly histogram the simulation data from each window
        # result: biased histogram and biased pdf for each window, shape (Nw, n_bins)
        self.xi_hist = []
        self.pdf_biased = []

        # step 2: combine windows using WHAM
        # results: free energy constants F_i, shape (Nw)
        # and the unbiased pdf pmf over all windows, shape (n_bins)
        self.exp_F_i = []
        self.pdf = []

        # step 3: convert unbiased pdf to unbiased pmf
        self.pmf = []

    # helper function for load_data
    def calc_skiprows_maxrows(self, t_min, t_max, t_start, delta_t):
        if t_min!=None and t_max!=None:
            t_min = float(t_min) # convert string from argparse to float
            t_max = float(t_max) # convert string from argparse to float

            skiprows = int((t_min - t_start) / delta_t) + 1
            max_rows = int((t_max - t_start) / delta_t) + 1 - skiprows
        
        elif t_min!=None and t_max==None:
            t_min = float(t_min) # convert string from argparse to float
            skiprows = int((t_min - t_start) / delta_t) + 1
            max_rows = None
        
        elif t_min==None and t_max!=None:
            t_max = float(t_max) # convert string from argparse to float
            skiprows = 0
            max_rows = int((t_max - t_start) / delta_t) + 1
        else:
            skiprows = 0
            max_rows = None
        
        return skiprows, max_rows
    
    # load data (list of values for coordinate (e.g. RMSD) over the entire simulation)
    # t_min, t_max specify the data selection to be LOADED
    def load_data(self, t_min=None, t_max=None, print_progress=False):
        if (t_min is not None and t_max is not None) and t_min > t_max:
            raise Exception('t_min > t_max')
                
        if print_progress:
            print('Loading coordinate data from data files')
        
        for i in range(self.Nw):
            # load the first 2 time values in the file efficiently
            t_firsttwo = np.loadtxt(self.windows[i,0], usecols=(0), max_rows=2)
            t_start = t_firsttwo[0]

            # assume that time values are evenly spaced
            delta_t = t_firsttwo[1] - t_firsttwo[0] # time between data points
            self.arr_delta_t.append(delta_t)

            # load the values for the coordinate xi (e.g. RMSD)
            skiprows, max_rows = self.calc_skiprows_maxrows(t_min, t_max, t_start, delta_t)
            arr_xi = np.loadtxt(self.windows[i,0], usecols=(1), skiprows=skiprows, max_rows=max_rows)

            t_start_loaded = np.loadtxt(self.windows[i,0], usecols=(0), skiprows=skiprows, max_rows=1)
            self.arr_t_start.append(t_start_loaded)

            if max_rows is not None:
                t_end_loaded = np.loadtxt(self.windows[i,0], usecols=(0), skiprows=skiprows+max_rows-1, max_rows=1)
            else:
                # load the last line of the file efficiently
                with open(self.windows[i,0], 'rb') as f:
                    try:  # catch OSError in case of a one line file
                        f.seek(-2, os.SEEK_END)
                        while f.read(1) != b'\n':
                            f.seek(-2, os.SEEK_CUR)
                    except OSError:
                        f.seek(0)
                    last_line = f.readline().decode()

                # store the last time value
                t_end_loaded = np.array([float(x) for x in last_line.split(' ')])[0]

            self.arr_t_end.append(t_end_loaded)

            # add the array of xi values to our list
            self.xi_data.append(arr_xi)
            
            if print_progress:
                print('Loaded '+str(arr_xi.size)+' from `'+self.windows[i,0]+'`')

            # raise error if empty
            if self.xi_data[i].size==0:
                raise Exception('Loaded data from `'+self.windows[i,0]+'` is empty. '
                                +'Please select t_min, t_max so that each window contains'
                                +'at least one datapoint with t_min <= t <= t_max')

    # t_min, t_max specify the data selection to be ANALYZED
    # so that a subset of the loaded data may be analyzed
    def select_data(self, t_min=None, t_max=None):
        if (t_min is not None and t_max is not None) and t_min > t_max:
            raise Exception('t_min > t_max')

        # clear previous selection
        self.xi_data_selected = []

        for i in range(self.Nw):
            if t_min is not None:
                index_start = max(0, int((t_min - self.arr_t_start[i])/self.arr_delta_t[i]))
            else:
                index_start = 0
            if t_max is not None:
                index_stop = min(self.xi_data[i].size,
                                 int((t_max - self.arr_t_start[i])/self.arr_delta_t[i]))
            else:
                index_stop = self.xi_data[i].size

            self.xi_data_selected.append(self.xi_data[i][index_start:index_stop])
                    
    # step 1: directly histogram the simulation data from each window
    # result: biased pdf for each window, shape (Nw, n_bins)
    def hist_data(self, xi_min=None, xi_max=None, n_bins=200, print_progress=False):
        if (xi_min is not None and xi_max is not None) and xi_min > xi_max:
            raise Exception('xi_min > xi_max')

        if xi_min is None:
            xi_min = min([np.min(self.xi_data_selected[i]) for i in range(self.Nw)])
        if xi_max is None:
            xi_max = max([np.max(self.xi_data_selected[i]) for i in range(self.Nw)])

        self.xi_binEdges = np.linspace(xi_min, xi_max, n_bins+1)
        self.xi_binCenters = (self.xi_binEdges[1:] + self.xi_binEdges[:-1])/2

        # window restraining potentials
        self.w_i = np.zeros((self.Nw, n_bins))
        for i in range(self.Nw):
            self.w_i[i,:] = 0.5*self.k_spring[i]*(self.xi_binCenters - self.window_centers[i])**2

        if print_progress:
            print('Histogramming the data with '+str(n_bins)
                  +' bin centers ranging from '+str(self.xi_binCenters[0])
                  +' to '+str(self.xi_binCenters[-1]))
        
        # histogram the data into the bins defined by xi_binEdges
        self.xi_hist = np.zeros((self.Nw, n_bins))
        self.pdf_biased = np.zeros((self.Nw, n_bins))
        total_n_datapoints = 0
        for i in range(self.Nw):
            self.xi_hist[i,:] = np.histogram(self.xi_data_selected[i], bins=self.xi_binEdges)[0]
            
            n_i = np.sum(self.xi_hist[i,:]) # number of independent data points
            total_n_datapoints += n_i

            if n_i > 0:
                # convert counts to probabilities
                self.pdf_biased[i,:] = self.xi_hist[i,:]/n_i

                # convert probabilities to probability densities
                self.pdf_biased[i,:] /= (self.xi_binEdges[1:] - self.xi_binEdges[:-1])

        if total_n_datapoints==0:
            raise Exception('Selected data is empty. '
                            +'Please select t_min, t_max, xi_min, xi_max so that there is '
                            +'at least one datapoint with t_min <= t <= t_max, xi_min <= xi <= xi_max')

    # step 2: combine the results from each window using WHAM
    # results: free energy constants F_i, shape (Nw)
    # and the unbiased pdf over all windows, shape (n_bins)
    def combine_windows(self, tolerance=0.0001, print_progress=False):
        if print_progress:
            print('Precomputing the numerator of Roux Eq. 8')
        
        arr_n_i = np.array([np.sum(self.xi_hist[i,:]) for i in range(self.Nw)])

        # numerator of Roux Eq 8 (constant with respect to F_i)
        numerator = np.einsum('i,ij->j', arr_n_i, self.pdf_biased)

        if print_progress:
            print('Precomputing useful constants for Roux Eq. 8 and 9')
        # useful constants (with respect to F_i) for Roux Eq 8 and 9
        n_bins = self.xi_binEdges.size-1
        exp_w_i = np.zeros((self.Nw, n_bins))
        denomFactor_i = np.zeros((self.Nw, n_bins))
        for i in range(self.Nw):
            exp_w_i[i,:] = np.exp(-self.w_i[i,:]/self.kT)
            denomFactor_i[i,:] = arr_n_i[i]*exp_w_i[i,:]

        exp_F_i = np.ones(self.Nw) # undetermined free energy constants F_i for each window
        
        if print_progress:
            print('Beginning WHAM iterations')
        # WHAM iterations to determine the F_i
        step = 0
        while True:
            exp_F_i_old = np.copy(exp_F_i)
            pdf = numerator/np.einsum('ij,i->j',denomFactor_i,1/(exp_F_i)) # Roux Eq 8
            exp_F_i = np.trapezoid(exp_w_i*pdf, x=self.xi_binCenters) # Roux Eq 9
            error = np.max(np.abs(exp_F_i-exp_F_i_old))
            step+=1
            if step%100==0 and print_progress:
                print('WHAM iteration '+str(step)+': '+str(error))
            if (error <= tolerance):
                if print_progress:
                    print('WHAM iteration '+str(step)+': '+str(error))
                break

        self.exp_F_i = exp_F_i # free energy constants F_i; shape (Nw)
        self.pdf = pdf # unbiased pdf; shape (n_bins)
        if print_progress:
            print('WHAM iterations finished')

    # convert unbiased pdf to unbiased pmf over all windows, shape (n_bins)
    def pdf_to_pmf(self):
        # make sure that we do not take log(0)
        mask_nonzero = (self.pdf!=0)
        mask_zero = (self.pdf==0)

        # definition of PMF (Potential of Mean Force) from PDF (Probability Density Function)
        pmf = np.copy(self.pdf)
        pmf[mask_nonzero] = -self.kT*np.log(self.pdf[mask_nonzero]) # kcal/mol

        # set values of pmf where pdf is zero to the max value of the pmf where the pdf is nonzero
        pmf[mask_zero] = np.max(pmf[mask_nonzero])
        
        pmf -= np.min(pmf) # move the minimum of the pmf to zero
        self.pmf = pmf

