#!/usr/bin/env python3

import warnings
import psd_utils

#for debugging purposes, treat warnings as exceptions
warnings.simplefilter('error')

min_freq = 0.1 #Hz, the cutoff frequency for the fitting
white_fit_freq = 50 #Hz, cutoff for the white frequency guess value
        #bins above this frequency is averaged to get an initial guess of the
        #white noise value for the fit
over_f_fit_freq = 1 #similarly, frequency bins below this frequency are fit
        #to a line in log scale to give an initial guess of the 1/f parameters
        #for the fit


if __name__ == '__main__':
    temps = [252, 275, 300, 325, 350, 375] #temperature values in mK
    temp_dict = {temp:'{0}_tsweep_psddata.sav'.format(temp) for temp in temps}
            #this should map the temps values to the psddata.sav file produced
            #by Steve's IDL code
    psd_utils.temp_data(temp_dict, white_fit_freq, min_freq, over_f_fit_freq)
