function grnoise,t,tc, N0=N0,V,tau_qp = tau_qp
; V in cubic microns
;if ~keyword_set(N0) then N0=1.72e10 ; microns^3 / eV
if ~keyword_set(N0) then N0=4.e10 ; microns^3 / eV
if ~keyword_set(tau_qp) then tau_qp = 5e-6 ; sec
if ~keyword_set(nqp_min) then nqp_min=400. ; QP per cubic micron at zero Temp

N0 = double(N0 / 1.6e-19) ; now microns^3 / Joule

; ef^2 = 4 beta^2 Nqp tau_qp
; beta = df_0 / d Nqp , so use (df_0/ dT) (dT / dNqp)

;Delta = double (3.5 * 1.381e-23 * tc)
Delta = double (1.74 * 1.381e-23 * tc)
;Nqp = v * 2 * N0 * sqrt(2*!pi*delta * 1.381e-23 * T) *exp(-1*Delta / (1.381e-23 * T)) + v * nqp_min

dNqp_dt = v * 2. * N0 * sqrt(2*!Pi * Delta * 1.318e-23)*exp(-1*Delta / (1.381e-23 * T)) * (1./(2.*sqrt(T)) + sqrt(T) * delta / 1.381e-23/t^2)
;stop

delta_t = T/30.

dNqp_dt = (Nqp(t+delta_t, tc,v)-Nqp(t-delta_t,tc,v)) / (2 * delta_t)

beta = df_response(t,tc,100.) / dNqp_dt
; assume a frequency of 100 MHz here, shouldn't matter.

ef2 = 4. * beta^2 * Nqp(t,tc,v) * tau_qp

return,ef2
end
