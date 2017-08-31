import numpy as np
from astropy import units as u
from astropy import constants as c

per_um3 = np.power(u.micron,-3)

def hello_world():
    print('Hello, world.  Making progress with Saul.')
    return

def Rx(t,tc,v,TiN=True,popt=0.*u.W,opt_freq=850*u.GHz,tauqp=1e-6*u.s,eta=0.5,nqp_min=400*per_um3):

    """ 
    eta is optical efficiency (?)
    If TiN is set, then the quasiparticle tauqp is fixed, otherwise it is
       calculated from the quasiparticle density qp_density
    """

    delta=(1.75*c.k_B * tc).to(u.J)
    
    qp_density_thermal = nqp(t,tc,v,nqp_min=nqp_min) / v

    return delta
    
    
