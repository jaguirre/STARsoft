function Rx,t,tc,v,tin=tin,popt=popt,opt_freq=opt_freq,tauqp=tauqp,eta=eta,nqp_min=nqp_min
;calculate fractional frequency responsivity in df/f per W
;v=volume in um^3

if ~keyword_set(eta) then eta=0.5
if ~keyword_set(tin) then tin=0
if ~keyword_set(nqp_min) then nqp_min=400.; QP per um^3
if ~keyword_set(popt) then popt=0. ; optical power in W
if ~keyword_set(opt_freq) then opt_freq = double(850.) ; GHz

delta = double(1.75*1.381e-23 * tc) ; J

qp_density_thermal = nqp(t,tc,v,nqp_min=nqp_min) / v ; this is the thermal QP density

if ~keyword_set(tauqp) then begin
   if tin then begin
      tauqp = 10.e-6            ; estimate for TiN
   endif else begin
      done=0 & threshold = 1e-3 & & count=0
      qp_density=qp_density_thermal
      tauqp= 8. / qp_density    ; per JZ
      while done eq 0 do begin
      ;   print,count,qp_density         
         nqp_photon = popt/(6.626e-25*opt_freq)*eta/v*tauqp
         qp_density_new=qp_density_thermal+nqp_photon
         if abs(alog10(qp_density_new/qp_density)) lt alog10(1.+threshold)then done=1
         qp_density=0.9*qp_density + 0.1 * qp_density_new
         tauqp=8./qp_density_thermal
       ;  count=count+1
       ;  stop
      end
      
      print,count
   end
   
end

stop

return,qp_density

end
