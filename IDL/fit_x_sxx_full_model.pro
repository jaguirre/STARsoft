;+
; NAME:
;       FIT_X_SXX_FULL_MODEL
;
; PURPOSE:
;       Given measurements of resonator fractional frequency shift Sxx as a function of incident power, do a fit using 'x_sxx_full_model.pro'
;
; CATEGORY:
;
; CALLING SEQUENCE:
;
; INPUTS:
;       pinc - the array of incident power points [pW]
;       
; OPTIONAL INPUTS:
;
; KEYWORD PARAMETERS:
;       xr - the array of measured fractional frequency shifts
;       sxx - the array of measured sxx values
;       sigma_xr - the array of measurement errors on xr
;       sigma_sxx - the array of measurement errors on sxx
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
;       10/19/2017 SHD
;                           
;-
;***********************************************************************
function fit_x_sxx_full_model, pinc, xr = xr, $
                                     sxx = sxx, $
                                     sigma_xr = sigma_xr, $
                                     sigma_sxx = sigma_sxx, $
                                     A = A, $
                                     pi = pi, $
                                     renormalize_chisqrd = renormalize_chisqrd, $
                                     quiet = quiet


; Configure the X array, Y array, and errors array
X = pinc
Y = [ [reform(xr)], [reform(sxx)] ]
errors = [ [reform(sigma_xr)], [reform(sigma_sxx)] ]

; Define pi structure if not specified
if (n_elements(pi) eq 0) then pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(A))

; Make sure the scale factor for dx_ref is fixed
pi[15].fixed = 1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize the output structure with some input parameters
struct = {X:X, Y:Y, Aguess:A, pi:pi, errors:errors}

; Do the fit
guess_model = x_sxx_full_model(X, A)
A = MPFITFUN('x_sxx_full_model', X, Y, errors, A, /quiet, PARINFO=pi, PERROR=fitsig)
model = x_sxx_full_model(X, A)
struct = create_struct(struct, 'modelprofile', 'x_sxx_full_model')

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