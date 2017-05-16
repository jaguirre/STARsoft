function nqp,t,tc,v,nqp_min=nqp_min

if ~keyword_set(nqp_min) then nqp_min=400. ; per cubic micron: zero-Temp residual QP density
; V is in cubic microns
;N0=1.72e10 for aluminum
N0=4.e10 ; for TiN
N0 = double(N0 / 1.6e-19) ; now microns^3 / Joule


;Delta = double (3.5 * 1.381e-23 * tc)
Delta = double (1.74 * 1.381e-23 * tc)
Nqp = v * 2 * N0 * sqrt(2*!pi*delta * 1.381e-23 * T) *exp(-1*Delta / (1.381e-23 * T)) + v * nqp_min

return,nqp
end
