import glob, os, pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splev, splrep
from scipy.integrate import odeint
import myconstants as const

from polytropes import monotrope, polytrope
from crust import SLyCrust
from eoslib import get_eos, glue_crust_and_core, eosLib

#import lalsimulation
#from pycbc import waveform
import matplotlib.pyplot as plt
plt.style.use('classic') 

c=const.CGS_C
G=const.CGS_G
Msun=const.CGS_MSUN

fm=1.e-13 #1femotometer in cm
dens_conversion=const.CGS_AMU/(fm**3)
edens_conversion=const.CGS_C**2

#--------------------------------------------------------

def save_obj(path, obj, name ):
    with open(path + name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(path, name ):
    with open(path + name + '.pkl', 'rb') as f:
        return pickle.load(f)
    
    
#---------------------------------------------------------




#-----------------------------------------------------
    



def symm_tidal_params(lambda1,lambda2,q):
    """
    Calculate best tidal parameters [Eqs. (5) and (6) in Wade et al. PRD 89, 103012 (2014)]
    Requires q <= 1
    """
    lambdap = lambda1 + lambda2
    lambdam = lambda1 - lambda2
    # Check that q <= 1, as expected
    if np.any(q > 1):
        raise ValueError("q > 1, while this function requires q <= 1.")
    dmbym = (1. - q)/(1. + q) # Equivalent to sqrt(1 - 4*eta) for q <= 1
    eta = q2eta(q)
    lam_tilde = (8./13.)*((1.+7.*eta-31.*eta*eta)*lambdap + dmbym*(1.+9.*eta-11.*eta*eta)*lambdam)
    dlam_tilde = (1./2.)*(dmbym*(1.-13272.*eta/1319.+8944.*eta*eta/1319.)*lambdap + (1.-15910.*eta/1319.+32850.*eta*eta/1319.+3380.*eta*eta*eta/1319.)*lambdam)
    return lam_tilde, dlam_tilde



################################ Ejecta Mass fits from Numerical Simulations ##################################################################
def Compactness_MBaryon(mass_input,eoskey):  #,pathEoS):

    #"""    
    mass_radius_file_loc='./m_mB_rad_data/'
    massradius_file=mass_radius_file_loc+eoskey+'.npy'
    mass_radius_rhoc_mbaryon=np.load(massradius_file)
    mass=np.array(mass_radius_rhoc_mbaryon[0])   #0:mass(msun), 1:radius(km), 2:rhoc(cgs), 3:lambda (in cgs)
    #if(mass_input>max(mass)): ##code it up as Assertion error
    #    raise ValueError('The input mass is greater than the maximum supported by the EOS:%s'%eoskey)
        
    Radius=np.array(mass_radius_rhoc_mbaryon[1])
    massBaryon=mass_radius_rhoc_mbaryon[3]
    lambda_dimensional = mass_radius_rhoc_mbaryon[4]

    sort_ind=mass.argsort()                 # sorting mass and other array; necessary for interpolation routine
    mass=mass[sort_ind[::1]]
    Radius=Radius[sort_ind[::1]]
    massBaryon=massBaryon[sort_ind[::1]]
    lambda_dimensional=lambda_dimensional[sort_ind[::1]]

    ## To check interpolation error
    #mass=mass[0::2]
    #Radius=Radius[0::2]
    #lambda_dimensional=lambda_dimensional[0::2]
    #massBaryon=massBaryon[0::2]
    
    radius_out=splev(mass_input,splrep(mass,Radius,k=3,s=0)) 
    mbaryon_out=splev(mass_input,splrep(mass,massBaryon,k=3,s=0))
    lambda_out=splev(mass_input,splrep(mass,lambda_dimensional,k=3,s=0))
    lambda_dimensionless = (lambda_out*const.CGS_G)*((const.CGS_C**2)/(const.CGS_G*mass_input*const.CGS_MSUN))**5

    return (const.CGS_G*mass_input*const.CGS_MSUN)/(radius_out*1.e5*const.CGS_C**2), lambda_dimensionless, radius_out   #from tov solver
    #return (const.CGS_G*mass_input*const.CGS_MSUN)/(radius_out*1.e5*const.CGS_C**2), mass_input+(13./200.)*mass_input**2  #from approx fit for baryon mass  #Give reference here
    