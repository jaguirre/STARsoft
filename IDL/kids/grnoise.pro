function grnoise,temp,tc,v,tin=tin,tau_qp = tau_qp,taumax=taumax, nstar=nstar,p_opt=p_opt,alpha=alpha,gamma=gamma
; V in cubic microns

if ~keyword_set(tin) then tin=0
if ~keyword_set(taumax) then taumax = 1.e-3 ; seconds
if ~keyword_set(nstar) then nstar=1.e2 ; per microns^3
if ~keyword_set(p_opt) then p_opt = 0. ; W
;if ~keyword_set(tau_qp) then tau_qp = 5e-6 ; sec
if ~keyword_set(nqp_min) then nqp_min=400. ; QP per cubic micron at zero Temp
if ~keyword_set(alpha) then begin
   if tin then alpha = 1.0 else alpha =0.1
endif
if ~keyword_set(gamma) then gamma = 1.

answer=fltarr(n_e(temp))

for jj=0,n_e(temp)-1 do begin
   t=temp[jj]

qpdensity = nqp_tot (t,tc,v,tin=tin,taumax=taumax,nstar=nstar) ; this is thermal QP density, without optical power

beta = dsigma2 (t,tc,nu=nu,tin=tin) / sigma2 (t,tc,nqp=qpdensity,nu=nu,tin=tin) * alpha * gamma / 2.

tau = tau_qp (nqp_tot(t,tc,v,tin=tin, taumax=taumax, nstar=nstar,p_opt=p_opt),taumax=taumax,nstar=nstar) ; tau includes the optical QPs

ef2 = 4. * beta^2 * qpdensity * tau
answer[jj]=ef2
;stop

endfor

return,answer
end
