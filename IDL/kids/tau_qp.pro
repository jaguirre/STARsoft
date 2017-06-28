function tau_qp,nqp,taumax=taumax, nstar=nstar


if ~keyword_set(taumax) then taumax = 1.e-3 ; seconds
if ~keyword_set(nstar) then nstar=1.e2 ; per microns^3

tau_qp = taumax / (1.+ nqp/nstar)
; Jonas review paper has taumax ~1e-3 and nstar ~100, which would
; extrapolate to 0.1 / n_qp in the high density limit
; Mauskopf model has tau = 10us x 1670  / n_qp = 0.0167 / n_qp
; so some discrepancy here.

return,tau_qp
end
