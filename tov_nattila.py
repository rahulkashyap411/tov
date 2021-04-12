#code is originally by Author: natj <nattila.joonas@gmail.com>; https://github.com/natj/tov
#Modified by Rahul.Kashyap, rahulkashyap411@gmail.com/1.Sept.2018: using this code to get mass-radius-baryon_mass relationship for all eos to be used in the calculation of the kilonovae light curve

import sys
import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import units as cgs
from math import pi
from polytropes import monotrope, polytrope
from crust import SLyCrust
from eoslib import get_eos, glue_crust_and_core, eosLib
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from label_line import label_line

import math

#from matplotlib import cm
import palettable as pal
cmap = pal.colorbrewer.qualitative.Set1_6.mpl_colormap
#cmap = pal.cmocean.sequential.Matter_8.mpl_colormap #best so far
#cmap = pal.wesanderson.Zissou_5.mpl_colormap

#--------------------------------------------------

def calc_tidal_deformability(C, Y,m):
    # """ Compute the dimensionless tidal deformability parameter Lambda from the compactness C and 
# the Lindblom y-potential at the surface of a polytropic star"""

# Eq.(C1,C2) of Lindblom & Indik 2014
    zeta = 4. * C**3 * (13. - 11.*Y + C*(3.*Y-2.) + 2.*(C**2)*(1.+Y)) + 3. * ((1.-2.*C)**2) * (2. - Y + 2.*C*(Y-1.)) *np.log(1.-2.*C) + 2. * C * (6. - 3.*Y + 3.*C*(5.*Y-8.))
    Lambda_dimensionless = (16./(15.*zeta)) * ((1.-2.*C)**2) * (2. + 2.*C*(Y-1.) - Y)  #dimensionless tidal deformability
    lambda_dimensional =   Lambda_dimensionless/cgs.G *(cgs.G*m*cgs.Msun/cgs.c**2)**5
    return lambda_dimensional

class tov:

    def __init__(self, peos):
        self.physical_eos = peos

    def tov(self, y,r):
        P, m, m_baryon, yp = y
        eden = self.physical_eos.edens_inv( P )
        rho = self.physical_eos.rho( P )
        eos_gamma = self.physical_eos.gamma_inv(P)

        dPdr = -cgs.G*(eden + P/cgs.c**2)*(m + 4.0*pi*r**3*P/cgs.c**2)
        dPdr = dPdr/(r*(r - 2.0*cgs.G*m/cgs.c**2))
        dmdr = 4.0*pi*r**2*eden
	dm_baryondr = 4.0*pi*r**2*rho*(1-2*cgs.G*m/(r*cgs.c**2))**(-0.5)
        G=cgs.G; c=cgs.c
        
        rho=eden*c**2
        dypdr= -yp**2/r -(r + (G/c**4)*4*np.pi*r**3*(P-rho))*yp/(r*(r-2*G*m/c**2)) + (G**2/c**4)*(4*(m+4*np.pi*r**3*P/c**2)**2)/(r*(r-2*G*m/c**2)**2) + 6/(r-2*G*m/c**2) - 4*np.pi*(r**2)*(5*rho+9*P+(rho+P)**2/(P*eos_gamma))*G/(c**4 * (r-2*G*m/c**2))

        return [dPdr, dmdr, dm_baryondr,dypdr]

    def tovsolve(self, rhoc):

        N = 2000
        r = np.linspace(1.e0, 1.8e6, N)
        #r = np.logspace(0.0,6.3,N)
        P = self.physical_eos.pressure( rhoc )
        eden = self.physical_eos.edens_inv( P )
        rho = self.physical_eos.rho( P )
        m = 4.0*pi*r[0]**3*eden
	m_baryon = 4.0*pi*r[0]**3*rho*(1-2*cgs.G*m/(r[0]*cgs.c**2))**(-0.5)
        yp=2.
        psol = odeint(self.tov, [P, m, m_baryon, yp], r, rtol=1.0e-6, atol=1.0e-4)
        #psol = solve_ivp(self.tov, [1.e0,5.8e6] ,[P, m, m_baryon, yp], method='RK45',t_eval=r )
        #print m, m_baryon, rhoc
        return r, psol[:,0], psol[:,1], psol[:,2], psol[:,3]
        #return psol.t, psol.y[0], psol.y[1], psol.y[2], psol.y[3] 

    def mass_radius(self):
        N = 800
        mcurve = np.zeros(N)
        rcurve = np.zeros(N)
	mbcurve = np.zeros(N)
        ypcurve = np.zeros(N)
        rhocs = np.logspace(12.5, 20.0, N)

        mass_max = 0.0
        j = 0
        for rhoc in rhocs:
            rad, press, mass, mass_baryon, yp_lambda = self.tovsolve(rhoc)

            rad  /= 1.0e5 #cm to km
            mass /= cgs.Msun
	    mass_baryon /=cgs.Msun
	    #print self.peos
	    #print rad
	    #print "Central Density=%f, Gravitational Mass=%f, Baryonic Mass=%f, Radius=%f"%(rhoc,mass.max(),mass_baryon.max(),rad.max())

            mstar = mass[-1]
            rstar = rad[-1]
            ypstar = yp_lambda[-1]
            for i, p in enumerate(press):
                if p > 0.0:
                    mstar = mass[i]
                    rstar = rad[i]
		    mbaryonStar = mass_baryon[i]
                    ypstar = yp_lambda[i]
            mcurve[j] = mstar
            rcurve[j] = rstar
	    mbcurve[j] = mbaryonStar
            
            C=(cgs.G/cgs.c**2)*(mstar*cgs.Msun)/(rstar*1.e5)
            ypcurve[j] = ypstar
            ypcurve[j] = calc_tidal_deformability(C,ypstar,mstar) 
	    


            j += 1
            if mass_max < mstar:
                mass_max = mstar
            else:
                break

        return mcurve[:j], rcurve[:j], rhocs[:j], mbcurve[:j], ypcurve[:j]


#--------------------------------------------------
def main(argv):

    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('axes', labelsize=7)
    

    fig = plt.figure(figsize=(3.54, 2.19)) #single column fig
    #fig = plt.figure(figsize=(7.48, 4.0))  #two column figure
    gs = plt.GridSpec(1, 1)

    ax = plt.subplot(gs[0, 0])
    ax.minorticks_on()
    ax.set_xlim(9.0, 16.0)
    ax.set_ylim(0.0, 3.0)

    ax.set_xlabel(r'Radius $R$ (km)')
    ax.set_ylabel(r'Mass $M$ (M$_{\odot}$)')

    
    # test single eos separately
    if False: 
        key = 'DD2'
        #key = 'BGN1H1'
        #key = 'ALF2'
        #key = 'ENG'
        #key = 'PAL6'

        dense_eos = get_eos(key)
        eos = glue_crust_and_core( SLyCrust, dense_eos )
        t = tov(eos)

        mass, rad, rhoc, massBaryon, lambda_dimensionless = t.mass_radius()
        #print mass
        #print rad
	#print rhoc
        #print massBaryon
        

        np.save('m_mB_rad_data/'+key,[mass,rad,rhoc,massBaryon,lambda_dimensionless])
	ax.plot(rad, mass)


    if True:
        i = 0
        for key, value in eosLib.iteritems():

            print "Solving TOVs for:", key, "(",i, ")"

            dense_eos = get_eos(key)
            eos = glue_crust_and_core( SLyCrust, dense_eos )
            t = tov(eos)
            mass, rad, rhoc, massBaryon, lambda_dimensionless = t.mass_radius()
            
            linestyle='solid'
            col = 'k'
            if value[4] == 'npem':
                col = 'k'
            if value[4] == 'meson':
                col = 'b'
            if value[4] == 'hyperon':
                col = 'g'
            if value[4] == 'quark':
                col = 'r'

            l, = ax.plot(rad, mass, color=col, linestyle=linestyle, alpha = 0.9)


            # labels for lines
            near_y = 1.45
            near_x = None
            rotation_offset=180.0
            if key == 'APR3':
                near_x = 11.0
                near_y = None
            if key == 'ENG':
                near_x = 11.2
                near_y = None
            if key == 'ALF2':
                near_y = 1.0
                rotation_offset = 0.0
            if key == 'MPA1':
                near_x = 12.0
                near_y = None
            if key == 'MS1b':
                near_y = 0.4
            if key == 'MS1':
                near_y = 2.3

            print l

	    np.save('m_mB_rad_data/'+key,[mass,rad,rhoc,massBaryon,lambda_dimensionless])	    


            label_line(l, key, near_y=near_y, near_x=near_x, rotation_offset=rotation_offset)


            i += 1


    #plot approximate central pressure isocurves
    #twice nuclear saturation density
    if False:
        x = [11.0, 15.8]
        y = [0.4, 2.5]
        ax.plot(x, y, "r--")

        txt = ax.text(15.5, 2.35, r'$2 \rho_n$', rotation=32, ha='center', va='center', size=8)
        txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=3))




if __name__ == "__main__":
    main(sys.argv)
    plt.subplots_adjust(left=0.15, bottom=0.16, right=0.98, top=0.95, wspace=0.1, hspace=0.1)
    plt.savefig('mr.pdf')



