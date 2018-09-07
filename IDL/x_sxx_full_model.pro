;+
; NAME:
;       X_SXX_FULL_MODEL
;
; PURPOSE:
;       Given an input vector of incident power values, compute resonator frequency shift (x = df/f, with specified absolute reference) and Sxx [Hz^-1].
;         The expression for Sxx includes all components of GR noise (photon generation, thermal generation, and recombination), as well as a fixed, power-independent component.
;       
;       This is calculated using the 'kid_responsivity.pro' routine
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;       X -- the independent variable, here the incident power [pW]
;       P -- an array controlling the values of each component
;         P[0] - temp
;         P[1] - tc
;         P[2] - N0
;         P[3] - nstar
;         P[4] - tau_max
;         P[5] - eta_pb
;         P[6] - vol
;         P[7] - fr
;         P[8] - alpha_k
;         P[9] - gamma_t
;         P[10] - nu_opt
;         P[11] - n_gamma
;         P[12] - use_te
;         P[13] - sxx_fixed [1e-17 Hz^-1]
;         P[14] - dx_ref
;         P[15] - scale factor for dx_ref (used to keep the fitted parameters of order 1)
;         P[16] - eta_opt (system optical efficiency)
;         
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;       The function defined by X and P
;
; RESTRICTIONS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:       
;       10/19/2017 SHD
;                           
;-
;***********************************************************************
function x_sxx_full_model, X, P


; Unpack the input parameters
pabs = X*P[16]
temp = P[0]
tc = P[1]
N0 = P[2]
nstar = P[3]
tau_max = P[4]
eta_pb = P[5]
vol = P[6]
fr = P[7]
alpha_k = P[8]
gamma_t = P[9]
nu_opt = P[10]
n_gamma = P[11]
use_te = P[12]
sxx_fixed = P[13]*1.d-17
dx_ref = P[14]*P[15]

; Call kid_responsivity.pro
result = kid_responsivity(temp, pabs, tc = tc, $
                                      N0 = N0, $
                                      nstar = nstar, $
                                      tau_max = tau_max, $
                                      eta_pb = eta_pb, $
                                      vol = vol, $
                                      fr = fr, $
                                      alpha_k = alpha_k, $
                                      gamma_t = gamma_t, $
                                      nu_opt = nu_opt, $
                                      n_gamma = n_gamma, $
                                      use_te = use_te)
                  
; Concatenate xr and sum of all noise terms, and return
output = [ [result.xr + dx_ref], [result.sxx_gr + result.sxx_gamma + sxx_fixed] ]
return, output


end
;***********************************************************************