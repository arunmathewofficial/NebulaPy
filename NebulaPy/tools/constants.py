"""
A set of physical constants, and constants strings.
"""

import numpy as np
from scipy import special

# known wave-length dictonary
wvl_dict = {}
wvl_dict['OIV25'] = 2.589332e+05

#
ev2Ang = 12.39841875e+3







################################################################
################################################################

planck = 6.6260693e-27   #erg s
planckEv = 4.13566743e-15  # ev s
hbar = planck/(2.*np.pi)
light = 29979245800.  # cm/s
q = 4.80320425e-10  # the units of charge

kev2Ang = 1.e-3*12.39841875e+3
ev2Erg = 1.602176487e-12
pi = 3.1415926535897931
boltzmann = 1.3806504e-16  # cgs
stefanBoltzmann = 5.670373e-5  # cgs - ergs cm^-2 K^-4 s^-1
boltzmannEv = 8.617343e-5
invCm2Ev = 1./8.06554465e+3
invCm2ryd = 1./109737.32
rydberg = 109737.31568 # cm^-1
rydbergErg = 2.1798723611035e-11 # erg
rydbergEv = 13.6056923  # eV
ryd2Ev = 13.6056923   # 1 rydberg = 13.604 eV
ryd2erg = 2.17987197e-11  #erg
fine = 7.2973525376e-3  # fine structure constant ~ 1./137  = e^2/(h_bar*c) h_bar = h/(2*pi)
emass = 9.10938215e-28  #  electron mass in gram
bohr = 0.52917720859e-8  # bohr radius in cm
hartree = 4.35974434e-11 #  erg
hartreeEv = 27.21138505
#
# derived constants
hc = planck*light
#
# area of bohr orbit
bohrCross = pi*bohr**2
#
std2fwhm = 2.*np.sqrt(2.*np.log(2.))
#
invCm2Erg = planck*light
#
boltzmannRyd = boltzmannEv/ryd2Ev
#
# collision produces the 8.63e-6 factor
collision = planck**2/((2.*pi*emass)**1.5*np.sqrt(boltzmann))
#
#
freeFree = 1.e+8*(light/(3.*emass))*(fine*planck/pi)**3*np.sqrt((2.*pi)/(3.*emass*boltzmann))
#
sutherland = (2./(3.*np.sqrt(3.)))*np.sqrt(pi/(2.*boltzmann*emass**3))*(planck*fine/pi)**3
#
freeBound = 1.e+8*(8.*fine*(planck**3))/(3.*np.sqrt(3.)*np.pi*(emass**4)*light)*(emass/(2.*np.pi*boltzmann))**1.5
#
freeBounde = 2./(4.*np.pi*planck*boltzmann*light**3*emass*np.sqrt(2.*np.pi*boltzmann*emass))
#
verner = (1.e-8/(planck*light**3*emass**3))*(emass/(2.*pi*boltzmann))**1.5

karzas = (2.**4)*planck*q**2/(3.*np.sqrt(3.)*emass*light)
#
freeFreeLoss = (8./3.)*np.sqrt(pi*boltzmann/(6.*emass**3))*(planck/pi)**2*fine**3
#
freeBoundLoss = ((16.*fine*(planck**2))/(3.*pi*np.sqrt(3.)*(emass**3)*(light**2)))*np.sqrt(emass/(2.*pi*boltzmann))
#
# astronomical
luminositySun = 3.86e+33 # ergs/s
radiusSun = 6.955e+10  # mean radius of Sun in cm
parsec = 3.08568025e+18  # cm
#
El = ['h','he','li','be','b','c','n','o','f','ne','na', \
    'mg','al','si','p','s','cl','ar','k','ca','sc','ti', \
    'v','cr','mn','fe','co','ni','cu','zn',\
    'ga','ge','as','se','br','kr']
Ionstage = ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII', \
    'XIV','XV','XVI','XVII','XVIII','XIX','XX','XXI','XXII','XXIII','XXIV', \
    'XXV','XXVI','XXVII','XXVIII','XXIX','XXX','XXXI','XXXII','XXXIII','XXXIV', \
    'XXXV','XXXVI','XXXVII']
Spd = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'O', 'Q', 'R', 'T', 'U', 'V', 'W', \
       'X','Y', 'Z', 'A','B', 'C', 'S1', 'P1', 'D1', 'E1', 'F1', 'G1', 'H1', 'I1', 'K1', 'L1', \
       'M1', 'N1', 'O1', 'Q1', 'R1', 'T1', 'U1', 'V1', 'W1','X1','Y1', 'Z1', 'A1','B1', 'C1']
#
#  data for Gauss-Laguerre integration
#
ngl = 12
zgl = special.roots_laguerre(ngl)
xgl = zgl[0]
wgl = zgl[1]

