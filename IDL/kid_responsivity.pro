;+
; NAME:
;       KID_RESPONSIVITY
;
; PURPOSE:
;       Compute the responsivity of a KID at finite temperature, and under a specified optical load. A number of associated quantities are computed as well. These are:
;        - nth -- number density of quasiparticles for a superconductor in thermal equilibrium at the specified temperature [micron^-3]
;        - nqp -- number density of quasiparticles for a superconductor at the specified temperature and optical load [micron^-3]
;        - tau_qp -- quasiparticle lifetime for a superconductor at the specified temperature and optical load [microsec]
;        - s1 -- S_1(w)
;        - s2 -- S_2(w)
;        - xr -- df/f, referenced to the nqp = 0 state
;        - Qi_inv -- Qi^-1
;        - r_x -- df/f responsivity [W^-1]
;        - r_qinv -- Qi^-1 responsivity [W^-1]
;        - sxx_gr -- GR noise, excluding the photon generation noise [Hz^-1]
;        - sxx_gr0 -- alternative, simpler form of sxx_gr often found in the literature [Hz^-1]
;        - sxx_gamma -- photon generation noise [Hz^-1]
;        
;        See 'KID_derivations.pdf' for details.
;       
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;       temp -- T [K]
;       pabs -- absorbed power [pW]
;       
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       tc -- Tc [K] {default is 1}
;       N0 -- single spin density of states at the Fermi level [micron^-3 eV^-1] {default is 1.72e10 for Al}
;       nstar - n* [micron^-3] {default is 100}
;       tau_max -- tau_max [microsec] {default is 100}
;       eta_pb -- pair-breaking efficiency {default is 0.7}
;       vol -- inductor volume [micron^3] {default is 1}
;       fr -- resonator frequency [MHz] {default is 100}
;       alpha_k -- kinetic inductance fraction {default is 1}
;       gamma_t -- thin film parameter {default is 1}
;       nu_opt -- optical frequency for photon noise computation [GHz] {default is 250}
;       n_gamma -- photon occupation number in the detector, for photon noise computation {default is 0}
;       use_te -- option to use the electron temperature to compute s1 and s2 {default is 1}
;
; OUTPUTS:
;       struct -- a structure with the following fields
;       .nth
;       .nqp
;       .tau_qp
;       .s1
;       .s2
;       .xr
;       .Qi_inv
;       .r_x
;       .r_qinv
;       .sxx_gr
;       .sxx_gr0
;       .sxx_gamma
;
; RESTRICTIONS:
;       The compute_te function currently crashes for T/Tc < 0.005
;
; EXAMPLE:
;
; MODIFICATION HISTORY:  
;       07/24/17 SHD
;       07/28/17 SHD -- replaced 'n0_gamma' with 'n_gamma'
;       10/17/17 SHD -- added keyword 'use_te', and included option of using electron temperature to compute s1 and s2
;       10/26/17 SHD -- modified the computation of electron temperature to make it accurate down to T/Tc = 0.005
;       11/17/17 SHD -- use lambertw() to compute the electron temperature, which works for T/Tc > 0.006
;
;-
;***********************************************************************
function my_te_root, t

; Insert the common block
common share1, target_value

; This is the function z = sqrt(t) * exp(-1/t) with a target value subtracted off
result = sqrt(t)*exp(-1./t) - target_value
return, result


end
;***********************************************************************


;***********************************************************************
function compute_te, z, frac_tolz=frac_tolz, frac_tolt=frac_tolt

; This function solves z = sqrt(t) * exp(-1/t), under the assumption that z is positive and real.

; As I understand it, the NEWTON root finding function will iterate until z is within tolz or x is within tolx.
; I have tried to define frac_tolz and frac_tolt such that they are the fractional uncertainties on z and t,
; respectively. Empirically, I have verified that with these parameters I get very accurate results, fractional
; uncertainties on t better than 1e-5 for T/Tc > 0.005, but the function crashes for smaller values of T

; Insert the common block
common share1, target_value

; Define keywords if not specified in call
if ~keyword_set(frac_tolz) then frac_tolz = 1.e-10
if ~keyword_set(frac_tolt) then frac_tolt = 1.e-4

; The first step is to construct an initial guess
; We use different models to compute the guess for different ranges of z
tguess = 0.*z
indx = where(z ge 0.367, count)
if (count ne 0) then tguess[indx] = 10.d^(poly(alog10(z[indx]), [ 0.36886585d, 1.0912276d,   0.64289548d,   0.065642612d,  -0.16116986d  ]))
indx = where(((z ge 1.44d-5) and (z lt 0.367)), count)
if (count ne 0) then tguess[indx] = 10.d^(poly(alog10(z[indx]), [ 0.30132981d, 0.81401834d,  0.26686403d,   0.049197816d,   0.0035917102d]))
indx = where(((z ge 1.d-40) and (z lt 1.44d-5)), count)
if (count ne 0) then tguess[indx] = 10.d^(poly(alog10(z[indx]), [-0.45614331d, 0.14869495d,  0.0085570243d, 0.00029349455d, 5.1709336d-06, 3.6055103d-08]))
indx = where(z lt 1.d-40, count)
if (count ne 0) then tguess[indx] = 10.d^(poly(alog10(z[indx]), [-1.4667896d,  0.015788392d, 9.6371179d-05, 3.1821936d-07,  4.1000813d-10]))

; For each value of z, then use a root-finding algorithm to solve for t
tfit = tguess
for i=0,n_elements(tfit)-1 do begin
  target_value = z[i]
  tfit[i] = newton(tguess[i], 'my_te_root', check = check, /double, tolf=frac_tolz*z[i], tolx=frac_tolt*tguess[i])
endfor
return, tfit


end
;***********************************************************************


;***********************************************************************
function kid_responsivity, temp, pabs, tc = tc, $
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
                                       use_te = use_te


; Define keywords if not specified in call
if ~keyword_set(tc) then tc = 1.
if ~keyword_set(N0) then N0 = 1.72e10
if ~keyword_set(nstar) then nstar = 100.
if ~keyword_set(tau_max) then tau_max = 100.
if ~keyword_set(eta_pb) then eta_pb = 0.7
if ~keyword_set(vol) then vol = 1.
if ~keyword_set(fr) then fr = 100.
if ~keyword_set(alpha_k) then alpha_k = 1.
if ~keyword_set(gamma_t) then gamma_t = 1.
if ~keyword_set(nu_opt) then nu_opt = 250.
if (n_elements(n_gamma) eq 0) then n_gamma = 0.
if (n_elements(use_te) eq 0) then use_te = 1

; Make every variable double precision
temp = double(temp)
pabs = double(pabs)
tc = double(tc)
N0 = double(N0)
nstar = double(nstar)
tau_max = double(tau_max)
eta_pb = double(eta_pb)
vol = double(vol)
fr = double(fr)
alpha_k = double(alpha_k)
gamma_t = double(gamma_t)
nu_opt = double(nu_opt)
n_gamma = double(n_gamma)

; Define various constants 
k_B = 1.381d-23 ;Boltzmann constant [J K^-1]
ev_joule = 1.6022d-19 ;eV/J ratio [eV J^-1]
h_p = 6.626d-34 ;Planck constant [J s]

; Compute the ratio of Delta_0/k_B [K]
d0_kB = 1.764*tc

; Compute nth
nth = 2.*N0*k_B/ev_joule*sqrt(2.*!pi*temp*d0_kB)*exp(-1.*d0_kB/temp)

; Compute nqp
; This expression has a term of the form [sqrt(1 + eps) - 1], where eps is small when nth and pabs are small.
; When eps is small this term is not computed accurately, and here we explicitly linearize it.
;nqp = nstar*(-1. + sqrt(1. + 2.*nth/nstar + (nth/nstar)^2. + 2.*eta_pb*pabs*1.e-12*tau_max*1.e-6/(nstar*vol*d0_kB*k_B)))
eps = 2.*nth/nstar + (nth/nstar)^2. + 2.*eta_pb*pabs*1.e-12*tau_max*1.e-6/(nstar*vol*d0_kB*k_B)
term = sqrt(1. + eps) - 1.
indx = where(eps lt 1.d-8, count)
if (count ne 0) then term[indx] = 0.5*eps[indx]
nqp = nstar*term

; Compute tau_qp
tau_qp = tau_max/(1. + nqp/nstar)

; Compute S1 and S2, optionally using the electron temperature
if use_te then begin
  nqp_norm = nqp/(2.*N0*k_B/ev_joule*d0_kB*sqrt(2.*!pi))
;  temp_electron = d0_kB*compute_te(nqp_norm)
;  temp_electron = d0_kB*2./my_lamw(2.*nqp_norm^(-2.))
  temp_electron = d0_kB*2./lambertw(2.*nqp_norm^(-2.))
  temp_for_s1s2 = temp_electron
endif else begin
  temp_for_s1s2 = temp
endelse
xi = h_p*fr*1.e6/(2.*k_B*temp_for_s1s2)
s1 = (2./!pi)*sqrt(2.*d0_kB/(!pi*temp_for_s1s2))*sinh(xi)*beselk(xi, 0)
s2 = 1. + sqrt(2.*d0_kB/(!pi*temp_for_s1s2))*exp(-1.*xi)*beseli(xi, 0)

; Compute xr and Qi_inv
; Note that xr refers to the frequency shift from the nqp = 0 state
xr = -1.*alpha_k*gamma_t*s2*nqp/(4.*N0*d0_kB*k_B/ev_joule)
Qi_inv = -1.*xr*2.*s1/s2

; Compute the frequency and Qinv responsivity
r_x = -1.*alpha_k*gamma_t*s2/(4.*N0*d0_kB*k_B/ev_joule)*eta_pb*tau_qp*1.e-6/(d0_kB*k_B*vol)
r_qinv = -1.*r_x*2.*s1/s2

; Compute Sxx_gr and Sxx_gr0
tau_th = tau_max/(1. + nth/nstar) ;quasiparticle lifetime for a superconductor in thermal equilibrium at the specified temperature [microsec]
gamma_th = nth*vol/2.*(1./tau_max + 1./tau_th)*1.e6 ;quasiparticle generation rate due to thermal fluctuations ;[sec^-1]
gamma_r = nqp*vol/2.*(1./tau_max + 1./tau_qp)*1.e6 ;quasiparticle recombination rate ;[sec^-1]
sxx_gr = (alpha_k*gamma_t*s2/(4.*N0*d0_kB*k_B/ev_joule))^2.*4.*(tau_qp*1.e-6)^2./vol^2.*(gamma_th + gamma_r)
sxx_gr0 = (alpha_k*gamma_t*s2/(4.*N0*d0_kB*k_B/ev_joule))^2.*4.*nqp/vol*tau_qp*1.e-6

; Compute Sxx_gamma
sxx_gamma = (r_x)^2.*2.*h_p*nu_opt*1.e9*pabs*1.e-12*(1. + n_gamma) 

; Define the output structure and return
struct = {nth:nth, $
          nqp:nqp, $
          tau_qp:tau_qp, $
          s1:s1, $
          s2:s2, $
          xr:xr, $
          Qi_inv:Qi_inv, $
          r_x:r_x, $
          r_qinv:r_qinv, $
          sxx_gr:sxx_gr, $
          sxx_gr0:sxx_gr0, $
          sxx_gamma:sxx_gamma}
return, struct


end
;***********************************************************************