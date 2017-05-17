#!/usr/bin/env python3

import pdb
import warnings

import numpy as np
from scipy.io import readsav
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib import pyplot as plt

#for debugging purposes, treat warnings as exceptions
warnings.simplefilter('error')

min_freq = 0.1 #Hz
white_fit_freq = 50 #Hz
over_f_fit_freq = 1

def flip_freqs(freq_array):
    split_point = len(freq_array) // 2
    split_point += 1
    
    return np.hstack((freq_array[split_point:], freq_array[:split_point]))
        

def is_symmetric(array):
    value_range = len(array) // 2
    for i in range(value_range):
        if array[i] != -1 * array[-1 - i]:
            return False

    return True


def read_data(file_dict):
    data = {key:readsav(f_name) for key, f_name in file_dict.items()}
    data = {key:[data[key][data_key] for data_key in data[key].keys()][0]
            for key in data.keys()} #readsav returns a 1 item dict
                    #because...IDL?

    return data


def get_sxx(sav_data):
    #slicing does not work on these memory map objects
    resonances = [tone for i, tone in enumerate(sav_data[0]) if i != 0]
            #first element is a string with the time
    sxx_sweeps = [tone.s_xx[0] for tone in resonances]
    freq_list = [tone.freq_vec[0] for tone in sxx_sweeps]
    Sxx_list = [tone.csd_xy[0] for tone in sxx_sweeps]

    return freq_list, Sxx_list


def get_S_par_per(sav_data):
    resonances = [tone for i, tone in enumerate(sav_data[0]) if i != 0]
    spar_sweeps = [tone.s_par[0] for tone in resonances]
    spar_freq = [tone.freq_vec[0] for tone in spar_sweeps]
    spar_csd = [tone.csd_xy[0] for tone in spar_sweeps]
    sper_sweeps = [tone.s_per[0] for tone in resonances]
    sper_freq = [tone.freq_vec[0] for tone in sper_sweeps]
    sper_csd = [tone.csd_xy[0] for tone in sper_sweeps]
    
    return sper_freq, sper_csd, spar_freq, spar_csd


def amp_noise_comp(sav_data):
    sper_freq, sper_csd, spar_freq, spar_csd = get_S_par_per(sav_data)
    #pull freqs over the cutoff
    sper_masks = [freq >= min_freq for freq in sper_freq]
    spar_masks = [freq >= min_freq for freq in spar_freq]
    sper_freq = [freq[mask] for freq, mask in zip(sper_freq, sper_masks)]
    sper_csd = [freq[mask] for freq, mask in zip(sper_csd, sper_masks)]
    spar_freq = [freq[mask] for freq, mask in zip(spar_freq, spar_masks)]
    spar_csd = [freq[mask] for freq, mask in zip(spar_csd, spar_masks)]
    #ratio of per/par
    ratio = [np.abs(sper_array) / np.abs(spar_array)
            for sper_array, spar_array in zip(sper_csd, spar_csd)]
    inv_ratio = [1. - tone_ratio for tone_ratio in ratio]

    return inv_ratio


def noise_fn(freq, log_white_level, exponent, amplitude):
    #1 over f
    over_f = amplitude * (freq**exponent)
    ret_val = over_f + np.exp(log_white_level)

    return ret_val


def log_noise_fn(freq, log_white_level, exponent, amplitude):
    #1 over f
    over_f = amplitude * (freq**exponent)
    ret_val = over_f + np.exp(log_white_level)

    return np.log(ret_val)


def starting_guess(freq, sxx):
    #just drop any invalid values
    freq = freq[sxx >= 0]
    sxx = sxx[sxx >= 0]
    log_sxx = np.log(sxx)
    log_freq = np.log(freq)
    exponent, log_amp, _, __, ___ = linregress(log_freq, log_sxx)

    return exponent, np.exp(log_amp)

    
def plot_fit(freq, sxx, filename, fit_values):
    binning_num = 10
    def bin_sxx(freq, sxx, bin_num):
        new_array = np.zeros(len(sxx) // bin_num)
        new_freq = np.zeros(len(sxx) // bin_num)
        for i in range(len(new_array)):
            new_array[i] = np.average(sxx[i*bin_num:(i+1)*bin_num])
            new_freq[i] = np.average(freq[i*bin_num:(i+1)*bin_num])
        return new_freq, new_array

    noise_fit = noise_fn(freq, *fit_values)
    fig = plt.figure()
    ax = plt.axes()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.plot(freq, noise_fit, 'g-')
    ax.plot(*bin_sxx(freq, sxx, binning_num), 'b-')
    ax.set_xlabel('Frequency(Hz)')
    ax.set_ylabel('S_xx(Hz^2)')
    fig.savefig(filename, format='png')
    plt.close(fig)

    
def fit_sxx(freq, sxx, temp, res_id):
    filename_template = 'fit_plots/temp_{0}_resonator_{1}_Sxx_fit.png'
    guess_white = sxx[freq > white_fit_freq]
    guess_white = np.log(np.average(guess_white[guess_white > 0.]))
    #maybe help numerically to bring things to order unity
    scale_factor = np.exp(guess_white)
    sxx /= scale_factor
    guess_white = 0.
    guess_exponent, guess_amp = starting_guess(freq[freq < over_f_fit_freq],
        sxx[freq < over_f_fit_freq])
    if (guess_exponent > 0) or np.isnan(guess_amp):
        guess_exponent = -0.01
    if (guess_amp < 0) or np.isnan(guess_amp):
        guess_amp = 0.
    sigma = np.sqrt(1. / freq) #I think that this is the appropriate weighting
            #since you will have 1/freq amount of cycles averaged
    fit_vals, _ = curve_fit(noise_fn, freq, sxx, p0=(guess_white,
            guess_exponent, guess_amp), bounds=((-np.inf, -np.inf, 0.),
                    (np.inf, 0, np.inf)), max_nfev=2000)
    #rescale back to actual values
    sxx *= scale_factor
    fit_vals[0] = fit_vals[0] + np.log(scale_factor)
    fit_vals[2] = fit_vals[2] * scale_factor
    plot_fit(freq, sxx, filename_template.format(temp, res_id), fit_vals)

    return fit_vals[0] #currently only looking at the white component


def log_fit(freq, sxx, temp, res_id):
    def log_fit_fn(freq, log_white_level, exponent, amplitude):
        return np.log(noise_fn(freq, log_white_level, exponent, amplitude))

    sxx = np.abs(sxx)
    filename_template = 'fit_plots/temp_{0}_resonator_{1}_Sxx_fit.png'
    guess_white = sxx[freq > white_fit_freq]
    guess_white = np.log(np.average(guess_white[guess_white > 0.]))
    #maybe help numerically?
    scale_factor = np.exp(guess_white)
    sxx /= scale_factor
    guess_white = 0.
    guess_exponent, guess_amp = starting_guess(freq[freq < over_f_fit_freq],
        sxx[freq < over_f_fit_freq])
    #sigma = np.sqrt(1. / freq) #I think that this is the appropriate weighting
            #since you will have 1/freq amount of cycles averaged
    fit_vals, _ = curve_fit(log_noise_fn, freq, np.log(sxx), p0=(guess_white,
            guess_exponent, guess_amp), maxfev=2000)
    sxx *= scale_factor
    fit_vals[0] = fit_vals[0] + np.log(scale_factor)
    fit_vals[2] = fit_vals[2] * scale_factor
    #fit_vals = (guess_white, guess_exponent, guess_amp)
    plot_fit(freq, sxx, filename_template.format(temp, res_id), fit_vals)

    return fit_vals[0] #currently only looking at the white component

def temp_data(temp_dict):
    data = read_data(temp_dict)
    data_sxx = [(key, get_sxx(temp_data)) for key, temp_data in data.items()]
    noise_ratio = {key: amp_noise_comp(temp_data)
            for key, temp_data in data.items()}
    temperatures = [temp for temp, data in data_sxx]
    noise_ratio = [noise_ratio[key] for key in temperatures]
    frequencies = [data[0] for temp, data in data_sxx]
    Sxx_list = [data[1] for temp, data in data_sxx]
    #change index order, so freq[resonator][temperature] instead of
    #[temperature][resonator]
    frequencies = [[frequencies[j][i] for j in range(len(frequencies))]
            for i in range(len(frequencies[0]))]
    Sxx_list = [[Sxx_list[j][i] for j in range(len(Sxx_list))]
            for i in range(len(Sxx_list[0]))]
    noise_ratio = [[noise_ratio[j][i] for j in range(len(noise_ratio))]
            for i in range(len(noise_ratio[0]))]
    #for right now, just use positive frequencies
    array_masks = [[freq >= min_freq for freq in tone] for tone in frequencies]
    frequencies = [[freq[mask] for freq, mask in zip(tone_freqs, tone_masks)]
            for tone_freqs, tone_masks in zip(frequencies, array_masks)]
    Sxx_list = [[np.abs(sxx[mask]) for sxx, mask in zip(tone_sxx, tone_masks)]
            for tone_sxx, tone_masks in zip(Sxx_list, array_masks)]
    Sxx_list = [[elem * ratio for elem, ratio in zip(tone_sxx, tone_ratios)]
            for tone_sxx, tone_ratios in zip(Sxx_list, noise_ratio)]
    white_noise_levels = [[fit_sxx(freq, sxx, temperatures[j], i+1)
            for j, (freq, sxx)
            in enumerate(zip(tone_frequencies, tone_Sxx_list))]
            for i, (tone_frequencies, tone_Sxx_list)
            in enumerate(zip(frequencies, Sxx_list))]
    with open('fit_vals.csv', 'w') as f:
        import csv
        writer = csv.writer(f)
        for resonator in white_noise_levels:
            writer.writerow(resonator)

    white_noise_levels = [[np.exp(value) for value in resonator]
            for resonator in white_noise_levels]
    plot_temp_white(temperatures, white_noise_levels)


def plot_temp_white(temps, white_noise_levels):
    filename_template = 'noise_vs_temp_resonator_{0}.png'
    for i, resonator in enumerate(white_noise_levels):
        fig = plt.figure()
        ax = plt.axes()
        ax.cla()
        ax.set_xlabel('Temperature(K)')
        ax.set_ylabel('White Noise Level(Hz)')
        ax.plot(temps, resonator)
        fig.savefig(filename_template.format(i+1), format='png')



