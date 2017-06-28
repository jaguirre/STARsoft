function sigma2,t,tc,nqp=nqp,nu=nu,tin=tin
; will be in units of sigma_n  Gao equation 16
; nqp in microns^-3

if ~keyword_set(nu) then nu=100. ; MHz
if ~keyword_set(tin) then tin=0
if ~keyword_set(nqp) then nqp = nqp_tot(t,tc,10.) ; calculate thermal QP density, assume arbitrary 10 micron volume 

if tin then N0=4.e10 else N0=1.72e10 ; in per microns^3 per eV

N0=double(N0 / 1.6e-19) ; now in per microns ^3 per J
D_0=double(1.762*1.381e-23* Tc) ; joules
xi = double (6.626e-34 * nu * 1.e6 / (2. *1.381e-23 * t) )

sigma2 = !Pi*D_0 / (6.626e-28 * nu) *(1. - nqp/(2.*N0*D_0)*(1.+ sqrt(2.*D_0/(!pi * 1.381e-23 * T))*exp(-1.*xi)*beseli(xi,0) ) )

return,sigma2
end
