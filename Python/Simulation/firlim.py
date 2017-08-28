import numpy as np
from astropy import units as u
from astropy import constants as c
from astropy.cosmology import Planck15 as cosmo

#%% The Bethermin 12 model
def sSFR_MS(z,Mstar):
    sSFR_MS0 = np.power(10.,-10.2)*np.power(u.yr,-1)
    beta_MS = -0.2
    gamma_MS = 3
    z_evo = 2.5
    zfac = z.copy() # Ugh, I don't understand Python passing
    zfac[z>z_evo] = z_evo
    zfac += 1
    return sSFR_MS0*np.power(Mstar/(1e11*u.Msun),beta_MS)*np.power(zfac,gamma_MS)

B_SB = 0.6
sigma_MS = 0.15
sigma_SB = 0.2

r_SB0 = 0.012
gamma_SB = 1.
z_SB = 1.

M_b = np.power(10.,11.2)*u.Msun
alpha = 1.3
phi_b1 = np.power(10.,-3.02)*np.power(u.Mpc,-3)
