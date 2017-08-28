import numpy as np
import astropy as ap
from astropy import cosmology
cosmo = cosmology.Planck13
from astropy import units as u
from astropy import constants as const

# Aw, fuck it.  Load up the Spinoglio data as a global variable. 
# Specify the format
datatype = {'names':('line','a','a_err','b','b_err','r','n','chisq'),'formats':('S10','f','f','f','f','f','i','f')}
spinoglio_data = np.loadtxt('spinoglio12_data.txt', dtype=datatype)

def xor(p, q):
    """ Re-inventin' the wheel """
    return ((p and not q) or (not p and q))

def line_ratio(lir, a, b):

    ratio = np.power(10.,a*np.log10(lir) + b) / lir
    return ratio

def lsol2ergs(L_solar=None,ergs=None):
    """ Convert luminosity in units of solar luminosity to erg/s.  I feel like astropy could help here ... """
    assert (xor(L_solar is None,ergs is None)),'Must set one of L_solar and erg/s'

    if (L_solar is not None):
        val = L_solar * 3.8388e26 * 1e7
    else:
        val = ergs / (3.8388e26 * 1e7)
    
    return val

def spinoglio(lines, L_IR_solar = None, L_line = None,all=False):
    """ Convert back and forth between L_IR and line luminosities using the Spinoglio et al 2011 fits """
    assert (xor(L_IR_solar is None, L_line is None)),' Must set one of L_IR_solar or L_line'
    val = {}
    if all:
        lines=spinoglio_data['line']
    
    # Ugh.  I made a bad choice in my initial configuration of this data structure
    for line in lines:
        indx = np.where(spinoglio_data['line'] == line)[0]
    
        if (L_IR_solar is not None):
            lir = lsol2ergs(L_solar=L_IR_solar) / 1e41
            log_L_line = spinoglio_data['a'][indx] * np.log10(lir) + spinoglio_data['b'][indx]
            #print 'log(L_line)', log_L_line
            L_line = lsol2ergs(ergs=np.power(10,log_L_line)*1e41)
            val[line] = L_line
        else:
            log_L_line = np.log10(lsol2ergs(L_solar=L_line) / 1e41)
            #print 'log(L_line)',log_L_line
            log_L_IR = (log_L_line - spinoglio_data['b'][indx])/spinoglio_data['a'][indx]
            #print log_L_IR
            L_IR_solar = lsol2ergs(ergs=np.power(10,log_L_IR)*1e41)
            val[line] = L_IR_solar
        
    return val

def spinoglio_line_ratios(lir_solar, upper = False, lower = False):

    """ Given an L_IR in solar luminosity, compute the line ratio relative to the given lir_solar, according to Spinoglio et al 2011 """

    # Currently only works with single lir_solar value, not vector
    assert (type(lir_solar) == float),'Does not work with vector lir input'
    # Convert to units of 10^41 erg/s
    lir = lsol2ergs(L_solar=lir_solar) / 1e41
    #lir_solar * 3.8388e26 * 1e7 / 1e41
    
    # The chi-squared gives an idea of how much of the scatter is
    # accounted for by the correlation.  Basically, want to scale up the
    # error bars (on the points) by the factor sqrt(chi-squared).  But how
    # does this propagate to the errors on fitted parameters?  Grrr ...

    if (upper):
        slope = a + a_err
        intercept = b + b_err

    if (lower):
        slope = a - a_err
        intercept = b - b_err

    if ((not (upper)) and (not(lower))):
        slope = spinoglio_data['a']
        intercept = spinoglio_data['b']

    ratios = {}
    for i in range(len(spinoglio_data['line'])):
        ratios[spinoglio_data['line'][i]] = line_ratio(lir,slope[i],intercept[i])

    return ratios

def delooze(L_CII=None,SFR=None):
    """ Calculate the De Looze et al 2011 SFR-L_CII relation, in either 
    direction.  L_CII is in solar luminosity, rather than the stupid erg/s """

    assert (xor(L_CII is None, SFR is None)),'Must set one of either L_CII or SFR'

    if (L_CII is not None):
        L_CII_ergs = lsol2ergs(L_solar=L_CII)
        SFR = np.power(L_CII_ergs,0.983)/1.028e40
        val = SFR
    else: 
        L_CII_ergs = np.power(SFR * 1.028e40,1/0.983)
        L_CII = lsol2ergs(ergs=L_CII_ergs)
        val = L_CII * u.solLum
        
    return val

def robertson_sfr(z):
    a = 0.01376 
    b = 3.26
    c = 2.59
    d = 5.68
    sfr = a * np.power(1+z,b)/(1 + np.power((1+z)/c,d))
    return sfr

def lumdens2surfbright(L_dens,z,lambda0):
    """ Convert luminosity per cosmological volume to Jy/sr """
    #Ldens = L_dens.to(u.watt/u.Mpc**3)
    Jy_per_sr = L_dens * Mpc3_per_srHz(z,lambda0) / (4 * np.pi * np.power(cosmo.luminosity_distance(z),2))
    Jy_per_sr = Jy_per_sr.to(u.Jansky)
    return Jy_per_sr

def Mpc3_per_srHz(z,lambda0):
    """ The conversion factor from Mpc^3 to sr Hz """
    # This is the comoving angular diameter distance per Bade's equation 2
    X2 = np.power(cosmo.angular_diameter_distance(z)*(1+z),2)
    Y = drdnu(z,lambda0)
    X2Y = X2 * Y
    X2Y = X2Y.to(u.Mpc**3 / u.Hertz)
    return X2Y

def drdnu(z,lambda0,approx=False):
    """ Calculate dr/dnu for a given redshift and wavelength """
    drdnu = lambda0 * np.power(1+z,2) / cosmo.H(z)
    drdnu = drdnu.to(u.Mpc/u.Hertz)
    if approx:
        O_m = cosmo.Om0
        O_l = 1 - O_m
        # Prefactor is only good for 21 cm line
        drdnu = 3.0 * np.power(1+z,2) / np.sqrt( O_m * np.power(1+z,3) + O_l )
        drdnu *= u.Mpc/u.megaHertz
    return drdnu
    
# From Furlanetto    
def r_of_z(z,delta_nu=0.1,omega_m=0.15):
    r_para = 1.7 /(u.Mpc) * (delta_nu/0.1) * np.power((1+z)/10.,0.5)*np.power(omega_m/0.15,-0.5)
    return r_para

def lfir_kennicutt(SFR):
    """ FIR luminosity given a star formation rate (M_sol/yr) """
    return 1.1e10 * SFR #* u.solLum

def V_vox(z, D, lambda0, R):
    """ Calculate the volume of a cosmological voxel at redshift z, given a solid angle, a rest wavelength, and a spectral resolution R = lambda/dlambda  """
#    assert (D.to(u.meter)),

    # Calculate observed wavelength
    lambda_obs = lambda0*(1+z)
    nu_obs = const.c/lambda_obs
    dnu = nu_obs/R
    # Calculate fwhm from telescope diameter, assuming diffraction limit
    theta = 1.22 * lambda_obs/D
    omega = np.power(theta,2)
    return (omega*Mpc3_per_srHz(z,lambda0)*dnu).to(np.power(u.Mpc,3))

def dr_los(z, lambda0, R):
    """ Calculate line of sight comoving distance """
    # Calculate observed wavelenth
    lambda_obs = lambda0*(1+z)
    nu_obs = const.c/lambda_obs
    dnu = nu_obs/R
    dr = drdnu(z,lambda0)*dnu
    return dr.to(u.Mpc)

def P_N(sigma_N,V_vox,t):
    return np.power(sigma_N,2) * V_vox / t

def nmodes(angular_extent,z0,dz,lambda0,D,R,dlnk):
    """Calculate the number of modes in a survey of linear angular extent
    theta (in degrees), centered at z0 with dz extent, observed with an
    instrument observing line lambda0 with a diffraction limited aperture
    D and resolving power R, binned into dlnk bins """
    # Calculate observed wavelength and frequency and dnu
    lambda_obs = lambda0*(1+z0)
    nu_obs = const.c/lambda_obs
    dnu = nu_obs/R
    # Calculate fwhm from telescope diameter, assuming diffraction limit
    theta = 1.22 * lambda_obs/D

    # Number of linear transverse pixels
    nkt = (np.round(angular_extent/theta))#.to(u.dimensionless_unscaled)

    return {'nkt':nkt,'theta':theta}
