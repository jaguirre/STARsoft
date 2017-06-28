function rx2,temp,tc,v,tin=tin,p_opt=p_opt,opt_freq=opt_freq,tauqp=tauqp,eta=eta,nqp_min=nqp_min,taumax=taumax, nstar=nstar,gamma=gamma,alpha=alpha

; try to use d sigma 2 /d nqp

if ~keyword_set(eta) then eta=0.65
if ~keyword_set(tin) then tin=0
if ~keyword_set(nqp_min) then nqp_min=400.; QP per um^3
if ~keyword_set(p_opt) then p_opt=0. ; optical power in W
if ~keyword_set(opt_freq) then opt_freq = double(850.) ; GHz
if ~keyword_set(taumax) then taumax = 1.e-3 ; seconds
if ~keyword_set(nstar) then nstar=1.e2 ; per microns^3
if ~keyword_set(gamma) then gamma=1.
if ~keyword_set(alpha) then begin
   if tin then alpha = 1. else alpha=0.1
endif

base_power = p_opt
answer = dblarr(n_e(temp))
for jj=0,n_e(answer)-1 do begin
   t=temp[jj]
   if p_opt ne 0 then dp = p_opt / 10. else dp = 3.e-17
   nqp_base = double(nqp_tot(t,tc,v,p_opt=base_power,taumax=taumax,nstar=nstar))
   done=0 & count=0
   while not done do begin
   nqp_high = double(nqp_tot(t,tc,v,p_opt=base_power+dp,taumax=taumax,nstar=nstar))
;   nqp_low  = double(nqp_tot(t,tc,v,p_opt=base_power-dp,taumax=taumax,nstar=nstar))
   
   if alog10(nqp_high/nqp_base) lt alog10(1.02) then dp=dp * 3. else done=1
   count=count+1
   endwhile

   dnqp_dp = (nqp_high-nqp_base ) / dp
 
   dsig_dnqp = (sigma2(t,tc,nqp=nqp_high)-sigma2(t,tc,nqp=nqp_base)) / (nqp_high-nqp_base)
;   dsig_dnqp = dsigma(t,tc,nqp
   dsig_dp = dsig_dnqp * dnqp_dp
   sig=sigma2(t,tc,nqp=nqp_base) 
   dff_dp = dsig_dp / sig / 2. * gamma * alpha

   answer[jj]=dff_dp 

endfor
;stop
return,answer

end

