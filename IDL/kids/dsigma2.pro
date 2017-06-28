function dsigma2,t,tc,nu=nu,tin=tin
; dsigma2 / dnqp -- will be in units of sigma_n
; nqp in microns^-3

if ~keyword_set(nu) then nu=100. ; MHz
if ~keyword_set(tin) then tin=0
if tin then N0=4.e10 else N0=1.72e10 ; in per microns^3 per eV
N0=double(N0 / 1.6e-19) ; now per microns^3 per J
D_0=double(1.762*1.381e-23* Tc) ; joules
xi = double (6.626e-34 * nu * 1.e6 / (2. *1.381e-23 * t) )

;dsigma2 = 1./(N0/1.6*6.626e-9*nu)*sqrt(2 * d_0/(!pi*1.381e-23*t))*sinh(xi)*beselk(xi,0)

dsigma2= !Pi / (2. * N0 * 6.626e-28 * nu)*(1. + sqrt(2.*D_0/(!pi*1.381e-23*t))*exp(-1.*xi)*beseli(xi,0) )


; !Pi*d_0 / (6.626e-28 * nu) *(1. - nqp/(2*N0/1.6e-19*d_0)*(1+ sqrt(2*D_0/(!pi * 1.381e-23 * T))*exp(-1*xi)*beseli(xi,0) ))


;stop

; THIS MAY NOT BE RELEVANT FOR THE CASE OF EXCESS
; QUASIPARTICLES -- don't use

return,dsigma2
end
