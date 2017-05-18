;+
; NAME:
;       ACQUIRE_AND_CALIBRATE_MULTITONE_DATA
;
; PURPOSE:
;       For a set of streaming data and associated calibration scans:
;        - Read in the raw data, optionally filter the time stream, and save (both as a .sav file, and selected fields as a .fits file)
;        - Fit the data to extract basic resonator parameters, and save results
;        - Calibrate the streaming data to directions parallel and perpendicular to the resonance circle, both in raw units (S21) and resonator units (f0_resonator), and save
;          --> We eventually want to extend this calibration to a delta(1/Q) as well, but this is not yet implemented
;        - Compute the psd of all time streams, and save
;        - Save postscript versions of plots tracking the calibration and psd computation
;        
;       The frequency, I, and Q for the gain calibration data, for channel n are obtained as (a similar structure for med and fine cal data):
;         f_gain = (multitone_data_raw.gain_caldata.data_calibration.tones)[n,*]*multitone_data_raw.bininhz
;         ivec = (multitone_data_raw.gain_caldata.data_calibration.i)[n,*]
;         qvec = (multitone_data_raw.gain_caldata.data_calibration.q)[n,*]
;       The frequency, I, and Q for the streaming data, for channel n, are obtained as:
;         f_stream = (multitone_data_raw.streamdata.stream_data_concat.tones)[n,0]*multitone_data_raw.streamdata.bininhz_stream
;         ivec = (multitone_data_raw.streamdata.stream_data_concat.s21i)[n,*]
;         qvec = (multitone_data_raw.streamdata.stream_data_concat.s21q)[n,*]
;       
;       A number of parameters in multitone_data_raw (nfreqs, bininhz, ntperf, tonebins, bins, blindbin, frlist) specify the placement of tones. It is assumed that this tone placement
;         is the same for all calibration datasets (gain, med, fine). The streaming data is allowed to have a different waveform length, and we record the (correspondingly different) bininhz
;         parameter as 'multitone_data_raw.streamdata.bininhz_stream'.
;         
;       By default, the conversion from IQ to resonator frequency is done by fitting a polynomial to the fine calibration IQ loop. The range of the IQ loop used for this calibration is the minimum satisfying the follow requirements:
;        - It contains the full range of streaming data falling sufficiently close to the IQ loop (as goverened by the rad_norm_range parameter)
;        - It is larger than a user-specified minimum phase (min_delta_phase) [radians]
;        - It encloses a minimum number of calibration IQ points (min_npoints)
;        
;       Additional sets of IQ to resonator frequency conversions are done using the fitted resonator models, rather than a polynomial fit to the fine calibration IQ loop.
;        - This is done if 'do_model_cal = 1'
;        - By default, all models aside from 'coarse' are used. The user can instead specify a list of model names with 'model_cal_names'.
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       multitone_data_raw -- a structure containing the raw, uncalibrated data
;       .date
;       .gain_caldata -- a structure containing the gain calibration data
;         .header_primary
;         .header_calsetup
;         .data_calsetup
;         .header_calibration
;         .data_calibration -- this is the field containing the calibration data; the frequency, I, and Q of channel n, are obtained as:
;          f_gain = (multitone_data_raw.gain_caldata.data_calibration.tones)[n,*]*multitone_data_raw.bininhz
;          ivec = (multitone_data_raw.gain_caldata.data_calibration.i)[n,*]
;          qvec = (multitone_data_raw.gain_caldata.data_calibration.q)[n,*] 
;       .rough_caldata --  a structure containing the rough calibration data (same fields as gain_caldata)
;       .med_caldata --  a structure containing the med calibration data (same fields as gain_caldata)
;       .fine_caldata -- a structure containing the fine calibration data (same fields as gain caldata)
;       .streamdata --   a structure containing the streaming data
;         .bininhz_stream -- bin of FFT, which sets the tones available for readout [Hz]
;         .averages_stream -- averaging number
;         .stream_data_sep
;         .stream_data_concat -- this is the field containing the streaming data for the full integration; the frequency, I, and Q of channel n are obtained as:
;          f_stream = (multitone_data_raw.streamdata.stream_data_concat.tones)[n,0]*multitone_data_raw.streamdata.bininhz_stream
;          ivec = (multitone_data_raw.streamdata.stream_data_concat.s21i)[n,*]
;          qvec = (multitone_data_raw.streamdata.stream_data_concat.s21q)[n,*]
;       .nfreqs -- number of tones used, including blind tones
;       .bininhz -- bin of FFT, which sets the tones available for readout [Hz]
;       .ntperf - number of comb tones per primary tone
;       .tonebins - spacing of tones (in bins) in a comb 
;       .bins -- list of bins corresponding to the primary tones; the frequency of each primary tone is bins*bininhz
;       .blindbin -- list of bins corresponding to blind primary tones
;       .frlist -- list of primary tones (should be equal to bins*bininhz)
;       multitone_data_calibrated -- a structure containing the calibrated data
;       .[comb]master_s21_stream_corr_rot -- series of S21 measurements, corrected for cable delay and amplitude wiggle, and shifted/rotated to put readout point on the positive, real axis, with center of IQ loop on the origin
;       .[comb]master_f_stream -- array of readout frequencies for on-resonance tones (blind tones are 0.), also the frequency used for the streaming data
;       .[comb]master_stream_phi_LS -- array of phi_LS for on-resonance tones
;       .[comb]master_f0_calculated -- array of inferred f_resonator measurements from the streaming data
;       .[comb]master_s21_stream_corr_rot_resampled - master_s21_stream_corr_rot resampled onto a regular time array (.time_vec_resampled)
;       .[comb]master_f0_calculated_resampled - master_f0_calculated resampled onto a regular time array (.time_vec_resampled)
;       .[comb]master_reso_params - fitted resonator parameters
;       .time_vec - time values [s] of the data, starting with 0
;       .time_vec_resampled - time_vec resampled to a regular grid and purged of any trailing zero-padding [s]
;
; RESTRICTIONS:
;       There are some upgrades/modifications I'd like to make to this:
;         - Include a calibration of the time stream into delta(1/Q), as well as delta(f0)
;         - Check that all data is recorded in double precision
;
; EXAMPLE:
;
; MODIFICATION HISTORY:       
;       10/09/2013 SHD
;       04/18/2014 SHD -- updated the documentation
;       05/16/2014 SHD -- added the psd calculation module
;       10/07/2015 SHD -- added 'dofit_reso_params' in call to calibrate_streaming_multitone
;       10/21/2015 SHD -- added 'maxnpts_stream_plot' in call to calibrate_streaming_multitone
;       10/21/2015 SHD -- modified call to 'compute_psd_multitone' to compute PSD for all channels
;       03/10/2016 SHD -- save the raw data as a .fits file, so it can read in with matlab
;       05/11/2016 SHD -- added 'rad_norm_range' keyword
;       06/07/2016 SHD -- various modifications to allow multiple tones per channel
;       08/19/2016 SHD -- modified handling of multiple tones per channel
;       08/29/2016 SHD -- added 'psd_use_unresampled' parameter
;       09/06/2016 SHD -- added 'bad_pts_indx_fine' and same for gain, rough, med, and modified call to 'cal_control_params'
;       09/27/2016 SHD -- added call to filter_streaming_multitone.pro to optionally filter data
;       04/23/2017 SHD -- replace the raw data with its complex conjugate
;       04/23/2017 SHD -- create outfolder and outfolder/plots/ immediately
;       04/23/2017 SHD -- added 'autoflag', 'do_flag_plot', 'do_ps_plots' keywords
;       04/23/2017 SHD -- define colorlist_x and colorlist_ps, colors for x windows and ps plots, include in 'colorlists' structure
;       04/23/2017 SHD -- add 'autoflag', 'do_flag_plot', 'do_ps_plots', 'outfolder', 'colorlists' in call to calibrate_streaming_multitone
;       04/23/2017 SHD -- add 'do_ps_plots', 'outfolder' in call to compute_psd_multitone
;       05/05/2017 SHD -- add 'do_model_cal' and 'model_cal_names', modify cal_control_params in call to calibrate_streaming_multitone
;
;-
;***********************************************************************
pro process_one_group, specfile


openr, lun, specfile, /get_lun

CDlabel = ''
devlabel = ''
datelabel_fine = ''
datelabel_gain = ''
datelabel_rough = ''
datelabel_med = ''
datelabel_stream = ''
timelabel_fine = ''
timelabel_gain = ''
timelabel_rough = ''
timelabel_med = ''
timelabel_stream = ''


readf,lun,CDlabel
readf,lun,devlabel
readf,lun,datelabel_fine
readf,lun,datelabel_gain
readf,lun,datelabel_rough
readf,lun,datelabel_med
readf,lun,datelabel_stream
readf,lun,timelabel_fine
readf,lun,timelabel_gain
readf,lun,timelabel_rough
readf,lun,timelabel_med
readf,lun,timelabel_stream
readf,lun,stage_temp
readf,lun,BB_temp

close,lun
free_lun,lun

;stop

; This block sets the main parameters that vary from run to run

; This block sets the parameters controlling the folders to read from and write to
;CDlabel = 'CD004'             ;used by this program to construct the name of the folder containing the saved data (outfolder)
;devlabel = 'atten20'             ;used by this program to construct the name of the folder containing the saved data (outfolder)
;datelabel_fine = '20170409'   ;date used by the multitone software to name the .fits file containing the fine data
;datelabel_gain = ''   ;date used by the multitone software to name the .fits file containing the gain data   ;if this is the same as 'datelabel_fine', then set this variable to a blank string
;datelabel_rough = ''  ;date used by the multitone software to name the .fits file containing the rough data  ;if this is the same as 'datelabel_fine', then set this variable to a blank string
;datelabel_med = ''    ;date used by the multitone software to name the .fits file containing the med data    ;if this is the same as 'datelabel_fine', then set this variable to a blank string
;datelabel_stream = '' ;date used by the multitone software to name the .fits file containing the stream data ;if this is the same as 'datelabel_fine', then set this variable to a blank string
;timelabel_fine = '234446'     ;timestamp used by the multitone software to name the .fits file containing the fine data;   if the cal data taken with the streaming data is to be used (rather than a separate, dedicated calibration), then set this variable to a blank string
;timelabel_gain = '235243'     ;timestamp used by the multitone software to name the .fits file containing the gain data;   if there is no gain scan, then set this variable to a blank string
;timelabel_rough = ''    ;timestamp used by the multitone software to name the .fits file containing the gain data;   if there is no rough scan, then set this variable to a blank string
;timelabel_med = '234825'      ;timestamp used by the multitone software to name the .fits file containing the med data;    if there is no medium scan, then set this variable to a blank string
;timelabel_stream = ''   ;timestamp used by the multitone software to name the .fits file containing the stream data; if there is no streaming scan, then set this variable to a blank string

;stage_temp = .257   ;stage temperature in K for these scans ; nothing is happening with this in IDL but they're included in the .csv with the fit results     
;BB_temp = 6.850      ;blackbody temperature in K for these scans ; nothing is happening with this in IDL but they're included in the .csv with the fit results

; This block sets the folders to read from and write to
home_dir = '/scr/starfire/'
data_dir = '/scr/starfire/labdata/'
infolder = data_dir + 'multitone/'                         ;folder containing the .fits files with the raw data
if (strtrim(timelabel_stream,2) ne '') then begin          ;this program saves to this folder, either based on streaming data file, or fine data file, in this order
  if (strtrim(datelabel_stream,2) ne '') then begin    
    local_outfolder = datelabel_stream + '_' + timelabel_stream 
  endif else begin
    local_outfolder = datelabel_fine + '_' + timelabel_stream 
  endelse
endif else begin  
  local_outfolder = datelabel_fine + '_' + timelabel_fine 
endelse
outfolder = home_dir + 'testdata/' + CDlabel + '/' + devlabel + '/multitone/' + local_outfolder + '/';this program saves to this folder
rawdata_filename = outfolder + 'rawdata.sav'               ;name of file in outfolder containing the raw data
rawdata_fits_filename = outfolder + 'rawdata_fits.fits'    ;name of .fits file in outfolder containing the raw data
calibrateddata_filename = outfolder + 'calibrateddata.sav' ;name of file in outfolder containing the calibrated data
psddata_filename = outfolder + 'psddata.sav'               ;name of file in outfolder containing the psd data
spawn, 'mkdir -p ' + outfolder
spawn, 'mkdir -p ' + outfolder + 'plots/'
    
; This block sets parameters for reading in the data
doreadrawdata = 1             ;need to read in the raw data? if set to 0, will read in the previously saved raw data (.save file) rather than going through the raw data files

; Optionally filter the streaming data before converting to f0
filter_streaming = 0          ;filter the streaming data?
fnotch = 60.*(dindgen(10)+1)  ;[Hz] array of frequencies to notch filter (set this to -1 to avoid notch filtering)
dfnotch = fnotch/60.          ;[Hz] the bandwidth of each notch filter
flow = 100.                   ;[Hz] cut-off frequency for a low-pass filter (set this to -1 to avoid low-pass filtering)
fhigh = 0.1                   ;[Hz] cut-on frequency for a high-pass filter (set this to -1 to avoid high-pass filtering)
order_low = -1                ;order of the low-pass filter (set this to -1 to make this a step filter)
order_high = -1               ;order of the high-pass filter (set this to -1 to make this a step filter)

; Optionally generate postscript versions of all calibration and psd plots
do_ps_plots = 1               ;make postscript plots?

; This block sets parameters for the IQ --> parallel/perpendicular (to resonator circle) calibration
docalstream = 1               ;need to calibrate the raw data? if set to 0, will read in the previously calibrated raw data rather than redoing the I/Q --> parallel/perpendicular calibration
autoflag = 1                  ;automatically flag bad phase points in the calibration scans? (requires at least one blind tone)
do_flag_plot = 0              ;make plots of bad phase point flagging in an x window? (postscript plot creation is determined by 'do_ps_plots')
do_gain_cabledelay_plots = 0  ;make plots of gain and cable delay fits in an x window? (postscript plot creation is determined by 'do_ps_plots')
dofit_reso_params = 1         ;specify whether to fit the resonator parameters
do_fit_plots = 0              ;make plots of reso fitting in an x window? (postscript plot creation is determined by 'do_ps_plots')
do_cal_plots = 0              ;make plots of IQ --> parallel/perpendicular calibration in an x window? (postscript plot creation is determined by 'do_ps_plots')
order_amp = 2                 ;order for poly fit to estimate amplitude variations in gain scan
min_delta_phase = 0.5         ;[radians] -- minimum range of theta to use for theta --> f calibration
min_npoints = 4               ;minimum number of fine scan points to use for theta --> f calibration
rad_norm_range = [0., 2.]     ;two-element array specifying minimum and maximum values; abs(S21_stream_corr_rot)/abs(S21_cal_corr_rot) must fall within this range for these streaming points to impact the theta --> f calibration
order_f_theta = 3             ;order for poly fit to frequency vs phase in fine IQ loop
maxnpts_stream_plot = 2000    ;when plotting the noise ball as check on calibration, only plot a subset of points
do_model_cal = 0              ;use resonator models to do the theta --> f calibration, in addition to the polynomial fitting?
model_cal_names = ['']        ;set of resonator models to use in the theta --> f calibration (e.g., ['fine', 'fine_seriesx']); if the first element is a blank string, all models are used aside from 'coarse'

; This block specifies a list of bad points in each of the calibration scans, to be ignored in the IQ --> frequency calibration
; If autoflag = 1, this list will be added to the list of auto-detected bad points
; If there are no bad points, the array should be [-1]
bad_pts_indx_fine = [-1]
bad_pts_indx_gain = [-1]
bad_pts_indx_rough = [-1]
bad_pts_indx_med = [-1]

; This block specifies a list of channels not to calibrate (there are problems, for example, if the resonator is too far into bifurcation)
problem_list = [192.6]*1.e6;[-1.]*1.e6 ;[Hz] -- list of channels not to calibrate (e.g., low-Q channels may be problematic)
problem_Qcut = 1.e3;50. ;tolerance for identifying problem channels (avoid channel if it's within plus/minus 1/2 of this Q from a listed problem channel)

; This block sets parameters for the tone comb analysis
tprimary = 0;8 ;specify which tone to use for fitting resonances 
tskip = [-1];[0, 1, 2, 3, 4, 5, 6, 7] ;specifies a list of tones (by index) in a comb not to calibrate individually (things break down if the tone does not catch the resonance during the fine sweep) {setting this to [-1] means no tones are discarded}
tskip_comp = [-1];[0, 1, 2, 3, 4, 5, 6] ;specifies a list of tones (by index) in a comb not to include in the simultaneous fitting {setting this to [-1] means no tones are discarded}

; This block sets the parameters for the psd calculation
dopsd = 1                ;gives the option of not analyzing the psd
psd_use_unresampled = 0  ;gives the option of working on the raw timestream, prior to resampling on a uniform time grid
do_psd_plots = 0         ;make plots of PSD computations?
psd_binfact = 10         ;factor by which the psd is binned before plotting
psd_freq_min = 0.        ;psd data at frequencies less than or equal to this value are discarded prior to binning
nosubtractedplot_psd = 0 ;if set to 1, this suppresses the plot of S_xx_bin_subtracted       
nolegend_psd = 1         ;if set to 1, this turns off the legend plotting for the psd plots  

; This block configures the plotting colors for x windows and ps plots
set_plot, 'X'
setcolors, /sys, /silent
colorlist_x = {gaincolor:!white, roughcolor:!green, medcolor:!blue, finecolor:!red, streamcolor:!cyan, offsetcolor:!orange, guesscolor:!green, fitcolor:!red}
set_plot, 'ps'
colorlist_ps = {gaincolor:cgcolor('black'), roughcolor:cgcolor('green'), medcolor:cgcolor('blue'), finecolor:cgcolor('red'), streamcolor:cgcolor('cyan'), offsetcolor:cgcolor('orange'), guesscolor:cgcolor('green'), fitcolor:cgcolor('red')}
colorlists = {colorlist_x:colorlist_x, colorlist_ps:colorlist_ps}
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This block reads in the raw calibration and streaming data, optionally filters the streaming data, and saves
if doreadrawdata then begin

; If necessary, this block sets the date labels for the calibration and streaming data to the fine date label
  if (strtrim(datelabel_gain,2) eq '')   then datelabel_gain   = datelabel_fine
  if (strtrim(datelabel_rough,2) eq '')  then datelabel_rough  = datelabel_fine
  if (strtrim(datelabel_med,2) eq '')    then datelabel_med    = datelabel_fine
  if (strtrim(datelabel_stream,2) eq '') then datelabel_stream = datelabel_fine

; This block packages up the date labels and time labels into a single structure for convenience
  date_time_labels = {datelabel_gain:datelabel_gain, $
                      datelabel_rough:datelabel_rough, $
                      datelabel_med:datelabel_med, $
                      datelabel_fine:datelabel_fine, $
                      datelabel_stream:datelabel_stream, $
                      timelabel_gain:timelabel_gain, $
                      timelabel_rough:timelabel_rough, $
                      timelabel_med:timelabel_med, $
                      timelabel_fine:timelabel_fine, $
                      timelabel_stream:timelabel_stream}

; This block reads in the raw data
  multitone_data_raw = read_rawdata_multitone(infolder, date_time_labels)

; This block replaces the raw data with its complex conjugate, correcting an error occuring somewhere in the data acquisition
  multitone_data_raw = conj_rawdata_multitone(multitone_data_raw)  

; This block filters the streaming data
  if (filter_streaming and tag_exist(multitone_data_raw, 'streamdata', /top_level)) then begin
    multitone_data_raw = filter_streaming_multitone(multitone_data_raw, fnotch=fnotch, $
                                                                        dfnotch=dfnotch, $
                                                                        flow=flow, $
                                                                        fhigh=fhigh, $
                                                                        order_low=order_low, $
                                                                        order_high=order_high)
  endif  
  
; This block saves the raw data  
  save, multitone_data_raw, filename = rawdata_filename

; This block saves some of the raw data to a .fits file  
  save_multitone_rawdata_asfits, multitone_data_raw, rawdata_fits_filename, tprimary
     
endif else begin
  restore, rawdata_filename
endelse
print, 'Done saving/restoring raw data  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This block converts the streaming data into offsets parallel and perpendicular to the resonance circle, both in raw units (Volts) and resonator units (f0_resonator)
; We eventually want to extend this calibration to a delta(1/Q) as well, but this is not yet implemented
if docalstream then begin

; This block packages up the parameters controlling the calibration for convenience
  cal_control_params = {autoflag:autoflag, $
                        do_flag_plot:do_flag_plot, $
                        do_gain_cabledelay_plots:do_gain_cabledelay_plots, $
                        dofit_reso_params:dofit_reso_params, $
                        do_fit_plots:do_fit_plots, $                        
                        do_cal_plots:do_cal_plots, $
                        order_amp:order_amp, $
                        min_delta_phase:min_delta_phase, $
                        min_npoints:min_npoints, $
                        rad_norm_range:rad_norm_range, $
                        order_f_theta:order_f_theta, $
                        maxnpts_stream_plot:maxnpts_stream_plot, $
                        do_model_cal:do_model_cal, $
                        model_cal_names:model_cal_names, $                        
                        problem_list:problem_list, $
                        problem_Qcut:problem_Qcut, $
                        tskip:tskip, $
                        tskip_comp:tskip_comp, $
                        bad_pts_indx_fine:bad_pts_indx_fine, $
                        bad_pts_indx_gain:bad_pts_indx_gain, $
                        bad_pts_indx_rough:bad_pts_indx_rough, $
                        bad_pts_indx_med:bad_pts_indx_med}
                      
; This block calibrates the data and saves                      
  multitone_data_calibrated = calibrate_streaming_multitone(multitone_data_raw, cal_control_params, do_ps_plots=do_ps_plots, outfolder=outfolder, colorlists=colorlists)
  save, multitone_data_calibrated, filename = calibrateddata_filename

; Identify real tones (not blind)
  bins = multitone_data_raw.bins
  blindbin = multitone_data_raw.blindbin
  junk = cgsetdifference(bins, blindbin, positions=reso_indx)
  non_blind = intarr(n_elements(bins))
  non_blind[reso_indx] = 1
  
  s = size(reso_indx)
  stage_temp = replicate(stage_temp,s[1],1)
  BB_temp_arr = replicate(BB_temp,s[1],1)
  write_csv, (outfolder+'fit_results_' + local_outfolder + '.csv'), $
                    (multitone_data_calibrated.master_reso_params.fitresults_fine_seriesx_nonlinear.fit_f0)[reso_indx], $
                    (multitone_data_calibrated.master_reso_params.fitresults_fine_seriesx_nonlinear.fit_QR)[reso_indx], $
                    (multitone_data_calibrated.master_reso_params.fitresults_fine_seriesx_nonlinear.fit_QC)[reso_indx], $
                    (multitone_data_calibrated.master_reso_params.fitresults_fine_seriesx_nonlinear.fit_QI)[reso_indx], $
                    (multitone_data_calibrated.master_reso_params.fitresults_fine_seriesx_nonlinear.fit_a_nl)[reso_indx], $ 
                    stage_temp, $
                    BB_temp_arr, $
                            HEADER=['f0','Qr','Qc','Qi','a_nl','T_stage','T_BB'], $
                            TABLE_HEADER=[local_outfolder]

endif else begin
  restore, calibrateddata_filename
endelse
print, 'Done saving/restoring calibrated data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; This block computes the PSDs of the noise trace projected into directions parallel and perpendicular to the resonance circle, both in raw units (Volts) and resonator units (df/f and dQ)
if (dopsd and tag_exist(multitone_data_raw, 'streamdata', /top_level)) then begin
  use_unresampled = psd_use_unresampled
  do_plots = do_psd_plots
  binfact = psd_binfact
  freq_min = psd_freq_min
  nosubtractedplot = nosubtractedplot_psd
  nolegend = nolegend_psd
  multitone_data_psd = compute_psd_multitone(multitone_data_raw, multitone_data_calibrated, do_ps_plots=do_ps_plots, $
                                                                                            problem_list = problem_list, $
                                                                                            problem_Qcut = problem_Qcut, $
                                                                                            tskip = tskip, $
                                                                                            use_unresampled = use_unresampled, $
                                                                                            do_plots = do_plots, $
                                                                                            binfact = binfact, $
                                                                                            freq_min = freq_min, $
                                                                                            nosubtractedplot = nosubtractedplot, $
                                                                                            nolegend = nolegend, $
                                                                                            outfolder = outfolder)                                                                                        
                                                                                            
  save, multitone_data_psd, filename = psddata_filename
  print, 'Done saving PSD data'
endif 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;stop
end
;***********************************************************************
