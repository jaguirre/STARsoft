function deltaf_f2,temp,tc,p_opt=p_opt,tin=tin,tauqp=tauqp,eta=eta,nqp_min=nqp_min,taumax=taumax, nstar=nstar

; calculate quasiparticle density and fractional frequency shift due
; to bath temperature, allow optical power

if ~keyword_set(eta) then eta=0.65
if ~keyword_set(tin) then tin=0.
if ~keyword_set(nqp_min) then nqp_min=400.; QP per um^3
if ~keyword_set(p_opt) then p_opt=0. ; optical power in W
if ~keyword_set(opt_freq) then opt_freq = double(850.) ; GHz
if ~keyword_set(gamma) then gamma=1.
if ~keyword_set(alpha) then begin
   if ~tin then alpha = 0.1 else alpha=1.                 ; 0.1 for aluminum
endif

if ~keyword_set(taumax) then begin
   if ~tin then taumax = 1.e-3 else taumax = 1.e-5 
endif
                                  ; seconds
if ~keyword_set(nstar) then begin 
   if ~tin then nstar=1.e2   else nstar = 1.e4 ; not sure what to use for TiN               ; per microns^3
endif

answer = fltarr(n_e(temp))
for jj=0,n_e(temp)-1 do begin
t=double(temp[jj])
nqp = nqp_tot(t,tc,tin=tin,p_opt=p_opt,opt_freq=opt_freq,tauqp=tauqp,eta=eta,nqp_min=nqp_min,taumax=taumax,nstar=nstar)
nqp0= nqp_tot(tc/20.,tc,tin=tin,p_opt=p_opt,opt_freq=opt_freq,tauqp=tauqp,eta=eta,nqp_min=nqp_min,taumax=taumax,nstar=nstar)

sig=sigma2(t,tc,nqp,tin=tin)
sig0=sigma2(tc/20.,tc,nqp0,tin=tin)

answer [jj]= 1/2. * alpha * (sig-sig0) / sig0
endfor
stop

return, answer
end
