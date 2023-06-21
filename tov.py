#!/usr/bin/env python
# coding: utf-8
#Mass-Radius calculation of neutron star for general EOS
#rahul.kashyap/21.Jul.2020: imported from kilonovastandardization project
__author__ = "Rahul Kashyap"
__maintainer__ = "Rahul Kashyap"
__email__ = "rahulkashyap411@gmail.com"

"""
Imported from natj <nattila.joonas@gmail.com>; https://github.com/natj/tov
Additions: 
(1) surface finding algorithm 
(2) EOS table format with interpolation
(3) calculation of tidal deformability 

#rahul/21.Jan.2019: imported from kilonovastandardization project
"""

import time
from tovlib import *
from scipy.integrate import solve_ivp
## Root finding by reducing the radial interval and by using scipy.rootfind functions. 
## Works and matches with the MATLAB code.. although there are problems of (1) Speed and (2) accuracy 
from scipy import optimize

#eos_file = './eos_tables/eosDD2F-SF.lorene'
#eos_file = './eos_tables/eosDD2.lorene'
#eos_file = './eos_tables/eos0.lorene'  #eosBHB
#eos_file = './eos_tables/BLH_new_14-Apr-2020.lorene'
eos_file = './eos_tables/eosSLy.lorene'

is_sorted = lambda x: (np.diff(x)>=0).all()

def eos_from_pres(pres_in,eos_file='./eos_tables/BLH_new_14-Apr-2020.lorene'):
    
    #if(pres_in.any()<0):
    #    print('pressure input is negative; stop integration and go to previous radial point')
    lg_pres_in = np.log10(pres_in)
    
    fm=1.e-13 #1femotometer in cm
    dens_conversion=const.CGS_AMU/(fm**3)
    edens_conversion=const.CGS_C**2
    
    ds=np.loadtxt(eos_file,comments='#',skiprows=9)
    #rho=ds[:,1]; edens=ds[:,2]; pres=ds[:,3]
    #print('all density > 0? ',all(i >= 0 for i in rho),'; all edens > 0? ',all(i >= 0 for i in edens),'; all pres > 0? ',all(i >= 0 for i in pres))
    #print('is_sorted(rho)? ',is_sorted(rho),'; is_sorted(edens)? ',is_sorted(edens),'; is_sorted(pres) ',is_sorted(pres))

    lg_rho, lg_edens, lg_pres = (np.log10(ds[:,1]),np.log10(ds[:,2]),np.log10(ds[:,3]))

    ind=lg_rho.argsort()                 # sorting mass and other array; necessary for interpolation routine
    sort_ind = ind[::1]
    lg_rho=lg_rho[sort_ind]
    lg_edens=lg_edens[sort_ind]
    lg_pres=lg_pres[sort_ind]
    
    
    #print(is_sorted(rho))
    #dp_dedens_arr=diff(pres)/diff(edens)/edens_conversion
    #dp_dedens_arr = np.insert(dp_dedens_arr, 0, dp_dedens_arr[0], axis=0)

    dlg_p_dlg_edens_arr=np.gradient(lg_pres,lg_edens) #https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
    gamma_arr = dlg_p_dlg_edens_arr
    Gamma_arr=((np.power(10,lg_edens)*edens_conversion+np.power(10,lg_pres))/(np.power(10,lg_edens)*edens_conversion)) * gamma_arr 
    
    rho_out=10**(splev(lg_pres_in,splrep(lg_pres,lg_rho,k=3,s=0)) )
    edens_out=10**(splev(lg_pres_in,splrep(lg_pres,lg_edens,k=3,s=0))) 
    Gamma_out=splev(lg_pres_in,splrep(lg_pres,gamma_arr,k=3,s=0))
    cs_out=np.sqrt((pres_in/edens_out)*(splev(lg_pres_in,splrep(lg_pres,gamma_arr,k=3,s=0))))  #sound speed: sqrt(dpres/dedens) = 
    #gamma_out=5./3.
    return rho_out*dens_conversion, edens_out, Gamma_out, cs_out


def eos_from_dens(rho_in,eos_file='./eos_tables/BLH_new_14-Apr-2020.lorene'):
    
    lg_rho_in = np.log10(rho_in)
    #print(lg_rho_in)
    fm=1.e-13 #1femotometer in cm
    dens_conversion=const.CGS_AMU/(fm**3)
    edens_conversion=const.CGS_C**2
    
    ds=np.loadtxt(eos_file,comments='#',skiprows=9)
    lg_rho, lg_edens, lg_pres = (np.log10(ds[:,1]*dens_conversion),np.log10(ds[:,2]),np.log10(ds[:,3]))
    #print(lg_rho)
    pres_out=10**(splev(lg_rho_in,splrep(lg_rho,lg_pres,k=3,s=0)) )
    edens_out=10**(splev(lg_rho_in,splrep(lg_rho,lg_edens,k=3,s=0))) 
    
    return pres_out, edens_out





###################################################################################3

fm=1.e-13 #1femotometer in cm
dens_conversion=const.CGS_AMU/(fm**3)
edens_conversion=const.CGS_C**2

"""
eos_dir = './eos_tables/'
eos_list=['eosDD2','eosDD2F-SF','eos0','eosSLy']  #['eosSLy'] #['eosSLy','eosDD2','eos0']  #['eosSLy'] #,'eosDD2']
c_list=['r','b','c','k','g']
p_new=np.logspace(np.log10(6.e17),np.log10(1.e38),300)
fig,((ax1,ax2),(ax3,ax4))=plt.subplots(nrows=2,ncols=2,figsize=(12,12))

for eoskey,color in zip(eos_list,c_list):
    print(f'############ {eoskey} ###########')
    eos_file = eos_dir+eoskey+'.lorene'
    eos_tabdata = np.loadtxt(eos_file,comments='#',skiprows=9)
    rho, edens, pres = (eos_tabdata[:,1],eos_tabdata[:,2],eos_tabdata[:,3])
    #print(eoskey,'\n####\n',rho,'\n#####\n',edens)
    #for p_in in p_new:
    rho_new, edens_new, gamma_new, cs_new = eos_from_pres(p_new,eos_file)
    
    print('is_sorted(edens_new)',is_sorted(edens_new),'is_sorted(rho_new)',is_sorted(rho_new))
    print(f'rho limits: min:{min(rho_new):1.3e}, max={max(rho_new):1.3e}')
    ax1.plot(rho*dens_conversion,pres,c=color,label=eoskey)
    ax1.plot(rho_new,p_new,c=color,marker='x',label='interpolated '+eoskey)
    #ax1.set_xlim([1.e14,1.e16])
    #ax1.set_ylim([1.e32,1.e36])
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set(xlabel=r'mass density ($g/cm^{3}$)',ylabel=r'pressure ($dyne/cm^2$)')
    ax1.legend(loc=4)

    ax2.plot(rho*dens_conversion,edens*edens_conversion,c=color,label=eoskey)
    ax2.plot(rho_new,edens_new*edens_conversion,c=color,marker='x',label='interpolated '+eoskey)
    #ax2.set_xlim([1.e14,1.e16])
    #ax2.set_ylim([1.e34,1.e37])
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set(xlabel=r'mass density ($g/cm^{3}$)',ylabel=r'edens ($g/cm^{3}$)')
    ax2.legend(loc=4)

    ax3.plot(rho_new,gamma_new,label=eoskey)
    ax3.set_xscale('log')
    ax3.set(xlabel=r'mass density ($g/cm^{3}$)',ylabel=r'$\Gamma$')
    ax3.legend(loc=4)
    
    ax4.plot(rho_new,cs_new/const.CGS_C,label=eoskey)
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.set(xlabel=r'mass density ($g/cm^{3}$)',ylabel=r'sound speed (in $c$)')
    
    
ax4.axhline(1/np.sqrt(3),ls='--',c='k',label='conformal limit cs: $c/\sqrt{3}$')
ax4.legend(loc=4)
plt.tight_layout()
plt.savefig('eos_plots_phase.png',dpi=150)
#plt.show()

"""
## saving eos data in a file in a specific format -- 
"""key='BLQ'
rho_new,edens_new,gamma_new=eos_from_pres(p_new,'./eos_tables/BLQ_gibbs_180_10-Mar-2020.lorene') #returns in cgs
index=np.array([int(i) for i in range(len(p_new))])
np.savetxt(f'./eos_tables/{key}_data.out',np.c_[index, rho_new/dens_conversion, edens_new/edens_conversion, p_new, gamma_new],header='\n#i rho edens pres gamma',comments=f'#density, pressure and energy density in cgs units for eos={key}.\n')

key='BLH'
rho_new,edens_new,gamma_new=eos_from_pres(p_new,'./eos_tables/BLH_new_14-Apr-2020.lorene') #returns in cgs
index=np.array([int(i) for i in range(len(p_new))])
np.savetxt(f'./eos_tables/{key}_data.out',np.c_[index, rho_new/dens_conversion, edens_new/edens_conversion, p_new, gamma_new],header='\n#i rho edens pres gamma',comments=f'#density, pressure and energy density in cgs units for eos={key}.\n')
"""

def calc_tidal_deformability(C, Y):
    # """ Compute the dimensionless tidal deformability parameter Lambda from the compactness C and 
    # the Lindblom y-potential at the surface of a polytropic star"""
    # Eq.(C1,C2) of Lindblom & Indik 2014
    zeta = 4. * C**3 * (13. - 11.*Y + C*(3.*Y-2.) + 2.*(C**2)*(1.+Y)) + 3. * ((1.-2.*C)**2) * (2. - Y + 2.*C*(Y-1.)) *np.log(1.-2.*C) + 2. * C * (6. - 3.*Y + 3.*C*(5.*Y-8.))
    Lambda_dimensionless = (16./(15.*zeta)) * ((1.-2.*C)**2) * (2. + 2.*C*(Y-1.) - Y)  #dimensionless tidal deformability
    #lambda_dimensional =   Lambda_dimensionless/const.CGS_G *(const.CGS_G*m*const.CGS_MSUN/const.CGS_C**2)**5
    return Lambda_dimensionless

def tov(r,y):
    #print('inside tov',eos_file)
    [P, m, m_baryon, yp] = y
    #if(P<0):
    #    sys.exit()
    #if eos_file != None:
    #P,dummy = eos_from_dens(rhoc,eos_file)
    rho,eden,eos_gamma, cs = eos_from_pres(P,eos_file)

    G=const.CGS_G; c=const.CGS_C 

    dPdr = -G*(eden + P/c**2)*(m + 4.0*np.pi*r**3*P/c**2)
    dPdr = dPdr/(r*(r - 2.0*G*m/c**2))
    dmdr = 4.0*np.pi*r**2*eden
    dm_baryondr = dmdr/np.sqrt(1-2*G*m/(r*c**2))
    #G=cgs.G; c=cgs.c

    rho=eden*const.CGS_C**2
    dypdr= -yp**2/r -(r + (G/c**4)*4*np.pi*r**3*(P-rho))*yp/(r*(r-2*G*m/c**2)) + (G**2/c**4)*(4*(m+4*np.pi*r**3*P/c**2)**2)/(r*(r-2*G*m/c**2)**2) + 6/(r-2*const.CGS_G*m/const.CGS_C**2) - 4*np.pi*(r**2)*(5*rho+9*P+(rho+P)**2/(P*eos_gamma))*G/(c**4 * (r-2*G*m/c**2))

    return [dPdr, dmdr, dm_baryondr,dypdr]

def tovsolve(rhoc,r_arr,eos_file='./eos_tables/eosSLy.lorene'):
    #print('inside tovsolve',eos_file)
    print('rhoc=',rhoc)
    P,dummy = eos_from_dens(rhoc,eos_file)
    rho,eden,Gamma,cs = eos_from_pres(P,eos_file)
    
    rad_low = r_arr[0]; rad_high = r_arr[-1]
    print(rad_low,rad_high)
    
    rmin = r_arr[0]
    r3=rmin**3
    m = 4./3.*np.pi*r3*eden
    m_baryon = 4./3.*np.pi*r3*eden*(1-2*const.CGS_G*m/(rmin*const.CGS_C**2))**(-0.5)
    yp=2.
    #psol = odeint(tov, [P, m, m_baryon, yp], r_arr, rtol=1.0e-6, atol=1.0e-4,tfirst=True)

    psol = solve_ivp(tov, [rad_low, rad_high] ,[P, m, m_baryon, yp,eos_file], method='RK45',t_eval=r_arr)
    #print m, m_baryon, rhoc
    #return r_arr, psol[:,0], psol[:,1], psol[:,2], psol[:,3]
    return psol.t, psol.y[0], psol.y[1], psol.y[2], psol.y[3] 


def find_surface(pmin, rhoc, rad_high,eos_file='./eos_tables/eosSLy.lorene'):
    print('inside find_surface',eos_file)
    print('rhoc=',rhoc)
    int_pts=2000
    rad_low=1.e-3
    r_arr = np.linspace(rad_low, rad_high,int_pts)
    #r = np.logspace(-4,6.3,N)

    #star = tovsolve(rhoc,r_arr)

    Pc,dummy = eos_from_dens(rhoc,eos_file)
    rhoc,eden_c,Gamma_c,cs_c = eos_from_pres(Pc,eos_file)
    #print(f'for pressure in: {Pc:1.2e}, density is {rhoc:1.2e}, edens: {eden_c:1.2e}, Gamma:{Gamma_c} and sound speed: {cs_c:1.2e}')
    rmin = r_arr[0]
    #print(f'rmin:{rmin:1.2e}')
    r3=rmin**3
    m = 4./3.*np.pi*r3*eden_c
    #print('g_00',np.sqrt(1-2*const.CGS_G*m/(rmin*const.CGS_C**2)))
    m_baryon = m/np.sqrt(1-2*const.CGS_G*m/(rmin*const.CGS_C**2))
    yp=2.
    
    var_vec=[Pc, m, m_baryon, yp]
    
    #psol = odeint(tov, var_vec, r_arr, rtol=1.0e-6, atol=1.0e-4,tfirst=True)
    #end_pres = psol[:,0][-1]
    #return end_pres-pmin

    psol = solve_ivp(tov, [rad_low, rad_high] ,var_vec, method='RK45')
    return psol.y[0][-1]-pmin
    ############################

    
    
    
###
time1=time.time()
print('eos_file=',eos_file)
eos_key=eos_file.split('./eos_tables/')[1].split('.lorene')[0]
f = open(f"tov_{eos_key}.txt", "w+")
f.write(f"#for eos = {eos_file} \n#Grav_Mass (solar mass) Radius (km) Lambda (dimensionless) rhoc (gm/cm^3) Compactness (dimensionless) Baryon_Mass (solar mass) \n")
f.close()
#######################################################3
pmin=1.e-10
len_seq=5
rhoc_arr=np.logspace(np.log10(6.e14),np.log10(5.e15),len_seq)

tov_data=[]
for rhoc in rhoc_arr:
    #rhoc=4.e14 #cgs
    time5= time.time()
    rstar = optimize.brentq(lambda rad_high: find_surface(pmin,rhoc,rad_high,eos_file=eos_file), 6.e5, 3.e6,rtol=1.e-4)
    #rstar = optimize.bisect(lambda rad_high: find_surface(pmin,rhoc,rad_high), 6.e5, 3.e6,rtol=1.e-5)
    
    #print(rstar/1.e5)

    time4= time.time()
    print('time elapsed in root finding:',time4-time5)

    rad_low = 1.e-3; rad_high = rstar
    int_pts=2000
    r_arr = np.linspace(rad_low, rad_high,int_pts)
    #r = np.logspace(-4,6.3,N)

    #star = tovsolve(rhoc,r_arr)

    Pc,dummy = eos_from_dens(rhoc,eos_file)
    rhoc,eden_c,Gamma_c,cs_c = eos_from_pres(Pc,eos_file)
    print(f'for pressure in: {Pc:1.2e}, density is {rhoc:1.2e}, edens: {eden_c:1.2e}, Gamma:{Gamma_c} and sound speed: {cs_c:1.2e}')
    rmin = r_arr[0]
    #print(f'rmin:{rmin:1.2e}')
    r3=rmin**3
    m = 4./3.*np.pi*r3*eden_c
    #print('g_00',np.sqrt(1-2*const.CGS_G*m/(rmin*const.CGS_C**2)))
    m_baryon = m/np.sqrt(1-2*const.CGS_G*m/(rmin*const.CGS_C**2))
    yp=2.
    psol = solve_ivp(tov, [rad_low, rad_high] ,[Pc, m, m_baryon, yp], method='RK45') #,t_eval=r_arr)
    [P, m, m_baryon, yp] = [psol.y[0][-1], psol.y[1][-1], psol.y[2][-1], psol.y[3][-1]]

    print(f'FINAL: rstar: {rstar/1.e5:1.2f}, grav. mass: {m/const.CGS_MSUN:1.2f}, bary. mass: {m_baryon/const.CGS_MSUN:1.2f}, yp: {yp}')
    C=(const.CGS_G/const.CGS_C**2)*m/rstar
    lambda_dimensionless = calc_tidal_deformability(C,yp)
    
    f = open(f"tov_{eos_key}.txt", "a")
    f.write(f"{m/const.CGS_MSUN:1.4f}  {rstar/1.e5:2.4f}  {lambda_dimensionless:2.4f}  {rhoc:1.4e}  {C:2.4f}  {m_baryon/const.CGS_MSUN:1.4f}\n")
    f.close()

    print(f'FINAL: rstar: {rstar/1.e5:1.2f}, grav. mass: {m/const.CGS_MSUN:1.2f}, bary. mass: {m_baryon/const.CGS_MSUN:1.2f}, yp: {yp}, lambda: {lambda_dimensionless:2.4f},   rhoc: {rhoc:1.4e}, compactness:  {C:2.4f}')
##
time2 = time.time()
print('total time elapsed:',time2-time1)

