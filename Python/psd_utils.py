#!/usr/bin/env python3

import warnings
import csv

import numpy as np
from scipy.io import readsav
from scipy.optimize import curve_fit
from scipy.stats import linregress
from matplotlib import pyplot as plt

#for debugging purposes, treat warnings as exceptions
warnings.simplefilter('error')


def noise_fn(freq, log_white_level, exponent, amplitude):
    #1 over f
    over_f = amplitude * (freq**exponent)
    ret_val = over_f + np.exp(log_white_level)

    return ret_val


def over_f_starting_guess(freq, sxx):
    #just drop any invalid values
    freq = freq[sxx > 0]
    sxx = sxx[sxx > 0]
    log_sxx = np.log(sxx)
    log_freq = np.log(freq)
    try:
        exponent, log_amp, _, __, ___ = linregress(log_freq, log_sxx)
    except (ValueError, RuntimeWarning):
        #sxx > 0 is null set
        exponent, log_amp = 0., 0. #hopefully error catching will drop these
            #later

    return exponent, np.exp(log_amp)


def fit_data(freq, csd, white_freq, min_freq, over_f_freq):
    guess_white = csd[freq > white_freq]
    if len(guess_white[guess_white > 0.]) > 0:
        guess_white = np.log(np.average(guess_white[guess_white > 0.]))
    else:
        guess_white = 0
    #improve numerical stability
    scale_factor = np.exp(guess_white)
    csd /= scale_factor
    guess_white = 0.
    guess_exponent, guess_amp = over_f_starting_guess(
            freq[freq < over_f_freq], csd[freq < over_f_freq])
    #reset invalid values
    if (guess_exponent > 0) or np.isnan(guess_amp):
        guess_exponent = -0.01
    if (guess_amp < 0) or np.isnan(guess_amp):
        guess_amp = 0.
    if len(csd[csd > 0]) == 0:
        return 0., 0.
    fit_vals, covariance = curve_fit(noise_fn, freq[freq > min_freq],
            csd[freq > min_freq], p0=(guess_white, guess_exponent, guess_amp),
            bounds=((-np.inf, -np.inf, 0.),
                    (np.inf, 0, np.inf)), max_nfev=2000)

    #rescale back to actual values
    csd *= scale_factor
    fit_vals[0] = fit_vals[0] + np.log(scale_factor)
    fit_vals[2] = fit_vals[2] * scale_factor

    return fit_vals, covariance


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


def amp_noise_comp(sav_data, white_freq, min_freq, over_f_freq):
    sper_freq, sper_csd, spar_freq, spar_csd = get_S_par_per(sav_data)
    #pull freqs over the cutoff
    sper_masks = [freq >= min_freq for freq in sper_freq]
    spar_masks = [freq >= min_freq for freq in spar_freq]
    sper_freq = [freq[mask] for freq, mask in zip(sper_freq, sper_masks)]
    sper_csd = [freq[mask] for freq, mask in zip(sper_csd, sper_masks)]
    spar_freq = [freq[mask] for freq, mask in zip(spar_freq, spar_masks)]
    spar_csd = [freq[mask] for freq, mask in zip(spar_csd, spar_masks)]
    #fit values
    sper_fit = [fit_data(freq, np.abs(csd), white_freq, min_freq, over_f_freq)[0]
            for freq, csd in zip(sper_freq, sper_csd)]
    spar_fit = [fit_data(freq, np.abs(csd), white_freq, min_freq, over_f_freq)[0]
            for freq, csd in zip(spar_freq, spar_csd)]
    #ratio of per/par
    sper_fit_csd = [noise_fn(freq, *sper_fit_vals)
            for freq, sper_fit_vals in zip(sper_freq, sper_fit)]
    spar_fit_csd = [noise_fn(freq, *spar_fit_vals)
            for freq, spar_fit_vals in zip(spar_freq, spar_fit)]
    ratio = [sper_array / spar_array
            for sper_array, spar_array in zip(sper_fit_csd, spar_fit_csd)]
    inv_ratio = [1. - tone_ratio for tone_ratio in ratio]
    for elem in inv_ratio:
        elem[elem < 0.] = 0.
                #zero points will be dropped later

    return inv_ratio


def read_data(f_name, white_freq, min_freq, over_f_freq):
    data = readsav(f_name)
    data = [data[key] for key in data.keys()][0]
    data_sxx = get_sxx(data)
    noise_ratio = amp_noise_comp(data, white_freq, min_freq, over_f_freq)
    frequencies = data_sxx[0]
    Sxx_list = data_sxx[1]

    return frequencies, Sxx_list, noise_ratio


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

    
def fit_sxx(freq, sxx, white_freq, min_freq, over_f_freq, temp, res_id):
    filename_template = 'fit_plots/temp_{0}_resonator_{1}_Sxx_fit.png'
    fit_vals, _ = fit_data(freq, sxx, white_freq, min_freq, over_f_freq)
    #rescale back to actual values
    plot_fit(freq, sxx, filename_template.format(temp, res_id), fit_vals)

    return fit_vals[0] #currently only looking at the white component


def plot_temp_white(temps, white_noise_levels):
    filename_template = 'noise_vs_temp_resonator_{0}.'
    for i, resonator in enumerate(white_noise_levels):
        fig = plt.figure()
        ax = plt.axes()
        ax.cla()
        ax.set_xlabel('Temperature(K)')
        ax.set_ylabel('White Noise Level(Hz)')
        #filter out defaulted white values when there were no points to fit
        plot_temps = [temp for temp, white in zip(temps, resonator)
                if white != 0. and white < 1e-8]
        resonator = [white for white in resonator if white != 0. and white < 1e-8]
        ax.plot(plot_temps, resonator, 'bo-')
        fig.savefig(filename_template.format(i+1) + 'png', format='png')
        with open(filename_template.format(i+1) + 'csv', 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['#temp', 'Sxx'])
            for temp, res in zip(plot_temps, resonator):
                writer.writerow([temp, res])

    colors = ['ro-', 'bo-', 'go-', 'yo-', 'ko-',
            'r*-', 'b*-', 'g*-', 'y*-', 'k*-',
            'rx-', 'bx-', 'gx-', 'yx-', 'kx-']
    fig = plt.figure()
    ax = plt.axes()
    ax.cla()
    ax.set_xlabel('Temperature(K)')
    ax.set_ylabel('White Noise Level(Hz)')
    for i, resonator in enumerate(white_noise_levels):
        #filter out defaulted white values when there were no points to fit
        plot_temps = [temp for temp, white in zip(temps, resonator)
                if white != 0. and white < 1e-8]
        resonator = [white for white in resonator if white != 0. and white < 1e-8]
        ax.plot(plot_temps, resonator, colors[i])
        ax.set_yscale('log')
    fig.savefig('all_resonators.png', format='png')

def temp_data(temp_dict, white_freq, min_freq, over_f_freq):
    data = [(key, f_name) for key, f_name in temp_dict.items()]
    temperatures = [data_elem[0] for data_elem in data]
    data = [read_data(data_elem[1], white_freq, min_freq, over_f_freq)
            for data_elem in data]
    frequencies = [data_tuple[0] for data_tuple in data]
    Sxx_list = [data_tuple[1] for data_tuple in data]
    noise_ratio = [data_tuple[2] for data_tuple in data]
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
    white_noise_levels = [[fit_sxx(freq, sxx, white_freq, min_freq,
            over_f_freq, temperatures[j], i+1)
            for j, (freq, sxx)
            in enumerate(zip(tone_frequencies, tone_Sxx_list))]
            for i, (tone_frequencies, tone_Sxx_list)
            in enumerate(zip(frequencies, Sxx_list))]

    white_noise_levels = [[np.exp(value) for value in resonator]
            for resonator in white_noise_levels]
    plot_temp_white(temperatures, white_noise_levels)



