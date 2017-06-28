function nqp_tot,temp,tc,v,tin=tin,p_opt=p_opt,opt_freq=opt_freq,tauqp=tauqp,eta=eta,nqp_min=nqp_min,taumax=taumax,nstar=nstar
; compute qp density including both thermally and optically-generated
; QPs.  calls tau_qp and solves numerically.
;v=volume in um^3

if ~keyword_set(eta) then eta=0.65
if ~keyword_set(tin) then tin=0
if ~keyword_set(nqp_min) then nqp_min=400.; QP per um^3
if ~keyword_set(p_opt) then p_opt=0. ; optical power in W
if ~keyword_set(opt_freq) then opt_freq = double(850.) ; GHz
if ~keyword_set(taumax) then begin
   if tin then taumax = 1.e-5 else taumax=1.e-3    ; sec
endif                  
if ~keyword_set(nstar) then begin
   if tin then nstar=1.e5 else nstar = 1.e2 * (1.e-3 / taumax) ;per cubic micron, don't know what to use for TiN 
endif

if tin then N0=4.e10 else N0=1.72e10 
N0=double(N0 / 1.6e-19)         ; now in per microns^3 per Joule

delta = double(1.762*1.381e-23 * tc) ; J

answer = fltarr(n_e(temp))
 threshold = 1.e-3

for jj=0,n_e(temp) -1 do begin
   t=temp[jj]
;   qp_density_thermal = nqp(t[jj],tc,v,nqp_min=nqp_min) / v ; this is the thermal QP density
   qp_density_thermal = 2 * N0 * sqrt(2. * !pi * delta * 1.381e-23 * T) * exp(-1*Delta / (1.381e-23 * T)) + nqp_min
   done=0 & count=0
   qp_density=qp_density_thermal
   tauqp = tau_qp(qp_density,taumax=taumax,nstar=nstar) ; use dedicated function
   while done eq 0 do begin
                                ;   print,count,qp_density         
      nqp_photon = p_opt * eta / delta * tauqp
      qp_density_new=qp_density_thermal+nqp_photon
      if abs(alog10(qp_density_new/qp_density)) lt alog10(1.+threshold)then done=1
      qp_density=0.9*qp_density + 0.1 * qp_density_new
                                ; tauqp = taumax / (1.+ qp_density/nstar)
      tauqp = tau_qp(qp_density,taumax=taumax,nstar=nstar)
      count=count+1
                                ;      stop
   endwhile
  ;                                  print,count
   answer[jj]=qp_density
endfor

;stop
;return,qp_density
return,answer
end
