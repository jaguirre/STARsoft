;+
; NAME:
;       FIT_SXX_FULL_MODEL
;
; PURPOSE:
;       Given measurements of Sxx vs temperature, do a fit using 'sxx_full_model.pro'
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;       temp - the array of temperature points [K]
;       
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       xr - the array of measured fractional frequency shifts
;       Qi_inv - the array of measured 1/Qi values
;       sigma_xr - the array of measurement errors on xr
;       sigma_Qi_inv - the array of measurement errors on Qi_inv
;       A - the vector of guess parameters
;       pi - the structure governing the fit
;       renormalize_chisqrd - if set, renormalize the uncertainties on the fitted parameters such that reduced chisqrd = 1
;       quiet - give the option of suppresing text output
;
; OUTPUTS:
;       struct - a structure with the fitting results
;
; RESTRICTIONS:
;
; EXAMPLE:
;
; MODIFICATION HISTORY:       
;       10/18/2017 SHD
;                           
;-
;***********************************************************************
function fit_sxx_full_model, temp, sxx = sxx, $
                                   sigma_sxx = sigma_sxx, $
                                   A = A, $
                                   pi = pi, $
                                   renormalize_chisqrd = renormalize_chisqrd, $
                                   quiet = quiet


; Configure the X array, Y array, and errors array
X = temp
Y = sxx
errors = sigma_sxx

; Define pi structure if not specified
if (n_elements(pi) eq 0) then begin
  pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(A))
  pi[*].fixed = 1
  pi[3:4].fixed = 0 ;vary pabs and tau_max
endif
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize the output structure with some input parameters
struct = {X:X, Y:Y, Aguess:A, pi:pi, errors:errors}

; Do the fit
guess_model = sxx_full_model(X, A)
A = MPFITFUN('sxx_full_model', X, Y, errors, A, /quiet, PARINFO=pi, PERROR=fitsig)
model = sxx_full_model(X, A)
struct = create_struct(struct, 'modelprofile', 'sxx_full_model')

; Calculate the fit chisqrd
dof = n_elements(Y)-total(1-pi.fixed)
chisq_reduced = total(((model-Y)/errors)^2.)/dof

; If necessary, renormalize the fit chisqrd to forced the reduced chisqrd to be 1
if (renormalize_chisqrd) then fitsig = fitsig*chisq_reduced^0.5

; Print the results
if (keyword_set(quiet) ne 1) then begin
  print, '
  print, 'Input A vector:'
  print, struct.Aguess
  print, 'Fixed parameters:'
  print, pi.fixed
  print, 'Fitted A vector:'
  print, A
  print, 'Fit sigma'
  print, fitsig
endif

; Add the fit data to the output structure and return
struct = create_struct(struct, 'guess_model', guess_model, 'A', A, 'fitsig', fitsig, 'model', model, $
  'dof', dof, 'chisq_reduced', chisq_reduced, 'renormalize_chisqrd', renormalize_chisqrd)
  
return, struct
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


end
;***********************************************************************