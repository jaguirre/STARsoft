#!/usr/bin/env python3

import csv
import glob
import sys

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def dB_to_lin(vals):
    return np.power(10., vals / 10.)


def lin_to_dB(vals):
    return 10. * np.log10(vals)


def read_file(filename):
    N_header_lines = 6
    freqs = []
    S_vals = []
    with open(filename, 'r') as f:
        #this throws out the header lines
        for i in range(N_header_lines):
            f.readline()
        #reading values
        reader = csv.reader(f)
        for row in reader:
            freqs.append(float(row[0]))
            S_21 = float(row[1]) + 1.j * float(row[2])
            S_vals.append(S_21)
    #convert to numpy arrays
    freqs = np.array(freqs)
    S_vals = np.array(S_vals, dtype=np.complex128)
    
    return freqs, S_vals


def read_dir(dirname):
    name_template = '/VNAslice*.txt'
    files = glob.glob(dirname + name_template)
    data = [read_file(filename) for filename in files]
    freqs = np.concatenate([freq_sub for freq_sub, S_sub in data])
    S_vals = np.concatenate([S_sub for freq_sub, S_sub in data])
    sort_indices = freqs.argsort()
    S_vals = S_vals[sort_indices]
    freqs = freqs[sort_indices]

    return freqs, S_vals


def non_lin_res(freqs, amp, f_0, a, Qc, Qi):
    angle_cutoff = 0.05
    Qr = 1. / ((1. / Qc) + (1. / Qi))
    y_0 = Qr * (freqs / f_0 - 1.)
    try:
        freqs[0]
    except IndexError:
        freqs = np.array([freqs])
    for i, detune in enumerate(y_0):
        roots = np.roots([4., -4. * detune, 1., -a - detune])
        
        real_roots = [root for root in roots
                if np.abs(np.imag(root) / np.real(root)) <= np.sin(
                angle_cutoff)]
        try:
            y_0[i] = np.abs(min(real_roots))
        except ValueError:
            print(roots)
            raise
    values = 1. - Qr / (Qc * (1. + 2.j * y_0))

    return amp * np.abs(values)


def remove_cryostat_atten(freqs, S_vals):
    def cable_delay(freq, cable_len, zero_phase):
        return (cable_len * freq + zero_phase + np.pi) % (2. * np.pi) - np.pi
    S_mag = np.abs(S_vals)
    N_poly = 21 #order of polynomial to fit to
    coeffs = np.polyfit(freqs, S_mag, N_poly)
    guess_vals = [-2 * np.pi / 1e8, 0.] #3m of cable
    fit_vals, _ = curve_fit(cable_delay, freqs, np.angle(S_vals), guess_vals,
            bounds=([-np.inf, 0.], [np.inf, 2 * np.pi]))
    cable_len, zero_phase = fit_vals
    
    return S_mag / (np.polyval(coeffs, freqs)
            * np.exp(1.j * cable_len * freqs + zero_phase))


def find_res(S_mag):
    N_window = 80
    min_depth = 0.8 #dB
    indices = []
    bot = N_window // 2
    for i in range(len(S_mag) - N_window):
        if S_mag[i+bot] == np.min(S_mag[i:i+N_window]):
            if lin_to_dB(np.max(S_mag[i:i+N_window]) /
                    np.min(S_mag[i:i+N_window])) >= min_depth:
                indices.append(i+bot)

    return indices
        

def fit_res(freqs, S_mag, f_min):
    amp = np.mean(S_mag)
    guesses = [amp, f_min, 0.1, 1e5, 1e4] #guess values
    try:
        fit_vals, _ = curve_fit(non_lin_res, freqs, S_mag, p0 = guesses,
                bounds=([0., 0., 0., 100., 100.], [np.inf, np.inf, 5., np.inf, np.inf]),
                maxfev=4000)
        #if freqs[0] > 1e8:
        #    plt.plot(freqs, S_mag, 'b.')
        #    plt.plot(freqs, non_lin_res(freqs, *fit_vals), 'r.')
        #    plt.show()
    except RuntimeError:
        fit_vals = [f_min, -1., 0., 0., 0.]
    return fit_vals


def write_fit_results(file_name, fit_list):
    with open(file_name, 'w') as f:
        f.write('#f_0, a, Qc, Qi\n')
        writer = csv.writer(f)
        for fit in fit_list:
            writer.writerow(fit)


def write_S_mag(freqs, S_vals, filename):
    with open(filename, 'w') as f:
        writer = csv.writer(f)
        for data_pt in zip(freqs, S_vals):
            writer.writerow(data_pt)
            

def main(directory):
    fit_window = 200
    freqs, S_vals = read_dir(directory)
    norm_S = remove_cryostat_atten(freqs, S_vals)
    norm_S = np.abs(norm_S)
    indices = find_res(norm_S)
    #fit_vals = [fit_res(freqs[i-fit_window:i+fit_window],
    #        norm_S[i-fit_window:i+fit_window], freqs[i]) for i in indices
    #        if (i>fit_window) and (i<len(freqs)-fit_window)]
    #fit_vals = [params[1:] for params in fit_vals] #don't care about amplitude
    #write_fit_results('fit_results.csv', fit_vals)
    res_freq = [freqs[i] for i in indices]
    with open('found_frequencies.csv', 'w') as f: 
        for freq in res_freq:
            f.write('{:1.4e}\n'.format(freq))
    res_S = [np.abs(S_vals[i]) for i in indices]
    write_S_mag(freqs, S_vals, 'sweep_data.csv')
    plt.plot(freqs, lin_to_dB(np.abs(S_vals)), 'b.')
    plt.plot(res_freq, res_S, 'r*')
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) == 1:
        print('Give a data directory')
        sys.exit(1)
    main(sys.argv[1])
