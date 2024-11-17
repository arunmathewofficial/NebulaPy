"""
A set of physical constants, and constants strings.
"""

import numpy as np
from scipy import special

########################################################################################
# Physical constants
pi = 3.1415926535897931
h  = 6.6260693e-27              #erg s
c  = 29979245800                # cm/s
kB = 1.3806504e-16              # cgs
alpha = 7.2973525376e-3  # fine structure constant ~ 1./137  = e^2/(h_bar*c) h_bar = h/(2*pi)
emass = 9.10938215e-28   #  electron mass in gram
q = 4.80320425e-10  # the units of charge

stefanBoltzmann = 5.670373e-5   # cgs - ergs cm^-2 K^-4 s^-1

########################################################################################
# Convertions factors
ev2Erg = 1.602176487e-12
ev2Ang = 12.39841875e+3
Ang2cm = 1.0e-8
kev2Ang = 12.39841875
invCm2Ev = 1./8.06554465e+3

verner = (1.e-8/(h*c**3*emass**3))*(emass/(2.*pi*kB))**1.5
# Astronomical constants
radiusSun = 6.955e+10    # mean radius of Sun in cm
luminositySun = 3.86e+33 # ergs/s
parsec = 3.08568025e+18  # cm


########################################################################################
# NebulaPy PION parameters
nebula_elements = {
    'H': "Hydrogen",
    'He': "Helium",
    'C': "Carbon",
    'N': "Nitrogen",
    'O': "Oxygen",
    'Ne': "Neon",
    'Si': "Silicon",
    'S': "Sulfur",
    'Fe': "Iron"}
mass = {
    'H': 1.6738e-24,
    'He': 6.6464768e-24,
    'C': 1.994374e-23,
    'N': 2.325892e-23,
    'O': 2.6567628e-23,
    'Ne': 3.3509177e-23,
    'Si': 4.6637066e-23,
    'S': 5.3245181e-23,
    'Fe': 9.2732796e-23}
atomic_number = {
    'H': 1,
    'He': 2,
    'C': 6,
    'N': 7,
    'O': 8,
    'Ne': 10,
    'Si': 14,
    'S': 16,
    'Fe': 26}
top_level_ions = {'H1+', 'He2+', 'C6+', 'N7+', 'O8+', 'Ne10+', 'Si14+', 'S16+', 'Fe26+'}
coordinate_system = {
    3: "spherical",
    2: "cylindrical",
    1: "cartesian"
}

########################################################################################
# NebulaPy Line Dictionary
wvl_dict = {}
wvl_dict['OIV25'] = 2.589332e+05



########################################################################################
# NebulaPy SED parameters

SEDModels = ['powr', 'atlas', 'blackbody']

# POWR Model Parameters
PoWRMetallicity = ['mw', 'lmc', 'smc', 'z0.07']
PoWRComposition = ['wne', 'wnl-h20', 'wnl-h40', 'wnl-h50', 'wnl-h60', 'wc']
PoWRMdotMin = -6.184
PoWRMdotMax = -3.784
ClumpFactor = {'mw-wne': 4, 'mw-wnl-h20': 4, 'mw-wnl-h50': 4, 'mw-wc': 10, 'lmc-wne': 10,
               'lmc-wnl-h20': 10, 'lmc-wnl-h40': 10, 'lmc-wc': 10, 'smc-wne': 4, 'smc-wnl-h20': 4,
               'smc-wnl-h40': 4, 'smc-wnl-h60': 4, 'smc-wc': 10, 'z007-wne': 10,
               'z007-wnl-h20': 10, 'z007-wnl-h40': 10, 'z007-wnl-h60': 10, 'z007-wc': 10, 'z086-wo': 0.4,
               'mw-ob-i': 10, 'lmc-ob-i': 10, 'smc-ob-vd3': 10, 'smc-ob-i': 10, 'smc-ob-ii': 10,
               'smc-ob-iii': 10}
data = {
    "Galactic Metallicity": {
        "MW WNE": {"log L": 5.3, "vfinal": 1600, "Dmax": 4, "XH": None, "XHe": None,
                   "XC": 0.98, "XN": 1.0E-4, "XO": 0.015, "XNe": None, "XFe": 1.4E-3},
        "MW WNL-H20": {
            "log L": 5.3,
            "vfinal": 1000,
            "Dmax": 4,
            "XH": 0.2,
            "XHe": 0.78,
            "XC": 1.0E-4,
            "XN": 0.015,
            "XO": None,
            "XNe": None,
            "XFe": 1.4E-3
        },
        "MW WNL-H50": {
            "log L": 5.3,
            "vfinal": 1000,
            "Dmax": 4,
            "XH": 0.5,
            "XHe": 0.48,
            "XC": 1.0E-4,
            "XN": 0.015,
            "XO": None,
            "XNe": None,
            "XFe": 1.4E-3
        },
        "MW WC": {
            "log L": 5.3,
            "vfinal": 2000,
            "Dmax": 10,
            "XH": None,
            "XHe": 0.55,
            "XC": 0.4,
            "XN": None,
            "XO": 0.05,
            "XNe": None,
            "XFe": 1.6E-3
        }
    },
    "LMC Metallicity": {
        "LMC WNE": {
            "log L": 5.3,
            "vfinal": 1600,
            "Dmax": 10,
            "XH": None,
            "XHe": 0.995,
            "XC": 7.0E-5,
            "XN": 4.0E-3,
            "XO": None,
            "XNe": None,
            "XFe": 7.0E-4
        },
        "LMC WNL-H20": {
            "log L": 5.3,
            "vfinal": 1000,
            "Dmax": 10,
            "XH": 0.2,
            "XHe": 0.795,
            "XC": 7.0E-5,
            "XN": 4.0E-3,
            "XO": None,
            "XNe": None,
            "XFe": 7.0E-4
        },
        "LMC WNL-H40": {
            "log L": 5.3,
            "vfinal": 1000,
            "Dmax": 10,
            "XH": 0.4,
            "XHe": 0.595,
            "XC": 7.0E-5,
            "XN": 4.0E-3,
            "XO": None,
            "XNe": None,
            "XFe": 7.0E-4
        },
        "LMC WC": {
            "log L": 5.3,
            "vfinal": 2000,
            "Dmax": 10,
            "XH": None,
            "XHe": 0.55,
            "XC": 0.4,
            "XN": None,
            "XO": 0.05,
            "XNe": 1.0E-3,
            "XFe": 7.0E-4
        }
    },
    "SMC Metallicity": {
        "SMC WNE": {
            "log L": 5.3,
            "vfinal": "1600‡",
            "Dmax": 4,
            "XH": None,
            "XHe": 0.998,
            "XC": 2.5E-5,
            "XN": 1.5E-3,
            "XO": None,
            "XNe": None,
            "XFe": 3.0E-4
        },
        "SMC WNL-H20": {
            "log L": 5.3,
            "vfinal": "1600‡",
            "Dmax": 4,
            "XH": 0.2,
            "XHe": 0.798,
            "XC": 2.5E-5,
            "XN": 1.5E-3,
            "XO": None,
            "XNe": None,
            "XFe": 3.0E-4
        },
        "SMC WNL-H40": {
            "log L": 5.3,
            "vfinal": "1600‡",
            "Dmax": 4,
            "XH": 0.4,
            "XHe": 0.598,
            "XC": 2.5E-5,
            "XN": 1.5E-3,
            "XO": None,
            "XNe": None,
            "XFe": 3.0E-4
        },
        "SMC WNL-H60": {
            "log L": 5.3,
            "vfinal": "1600‡",
            "Dmax": 4,
            "XH": 0.6,
            "XHe": 0.398,
            "XC": 2.5E-5,
            "XN": 1.5E-3,
            "XO": None,
            "XNe": None,
            "XFe": 3.0E-4
        },
        "SMC WC": {
            "log L": 5.3,
            "vfinal": 2000,
            "Dmax": 10,
            "XH": None,
            "XHe": 0.547,
            "XC": 0.4,
            "XN": None,
            "XO": 0.05,
            "XNe": 2.4E-3,
            "XFe": 3.0E-4
        }
    },
    "sub-SMC Metallicity (0.07 solar)": {
        "Z0.07 WNE": {
            "log L": 5.3,
            "vfinal": 1600,
            "Dmax": 10,
            "XH": None,
            "XHe": 0.999,
            "XC": 1.0E-5,
            "XN": 6.1E-4,
            "XO": 1.0E-5,
            "XNe": None,
            "XFe": 9.2E-5
        },
        "Z0.07 WNL-H20": {
            "log L": 5.3,
            "vfinal": 1600,
            "Dmax": 10,
            "XH": 0.2,
            "XHe": 0.799,
            "XC": 1.0E-5,
            "XN": 6.1E-4,
            "XO": 1.0E-5,
            "XNe": None,
            "XFe": 9.2E-5
        },
        "Z0.07 WNL-H40": {
            "log L": 5.3,
            "vfinal": 1600,
            "Dmax": 10,
            "XH": 0.4,
            "XHe": 0.599,
            "XC": 1.0E-5,
            "XN": 6.1E-4,
            "XO": 1.0E-5,
            "XNe": None,
            "XFe": 9.2E-5
        },
        "Z0.07 WC": {
            "log L": 5.3,
            "vfinal": 2000,
            "Dmax": 10,
            "XH": None,
            "XHe": 0.549,
            "XC": 0.4,
            "XN": None,
            "XO": 0.05,
            "XNe": 8.3E-4,
            "XFe": 9.2E-5
        }
    }
}

# Atlas Model Parameter
AtlasMetallicity = ['-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-0.0', '0.0',
                    '+0.0', '+0.2', '0.2', '+0.5', '0.5']

AtlasGravity = ['0.0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0',
                '4.5', '5.0', '+0.0', '+0.5', '+1.0', '+1.5', '+2.0', '+2.5',
                '+3.0', '+3.5', '+4.0', '+4.5', '+5.0']

# Blackbody Model Parameters
blackbody_temp_table = [3.500000e+03, 3.750000e+03, 4.000000e+03, 4.250000e+03,
                        4.500000e+03, 4.750000e+03, 5.000000e+03, 5.250000e+03,
                        5.500000e+03, 5.750000e+03, 6.000000e+03, 6.250000e+03,
                        6.500000e+03, 6.750000e+03, 7.000000e+03, 7.250000e+03,
                        7.500000e+03, 7.750000e+03, 8.000000e+03, 8.250000e+03,
                        8.500000e+03, 8.750000e+03, 9.000000e+03, 9.250000e+03,
                        9.500000e+03, 9.750000e+03, 1.000000e+04, 1.025000e+04,
                        1.050000e+04, 1.075000e+04, 1.100000e+04, 1.125000e+04,
                        1.150000e+04, 1.175000e+04, 1.200000e+04, 1.225000e+04,
                        1.250000e+04, 1.275000e+04, 1.300000e+04, 1.400000e+04,
                        1.500000e+04, 1.600000e+04, 1.700000e+04, 1.800000e+04,
                        1.900000e+04, 2.000000e+04, 2.100000e+04, 2.200000e+04,
                        2.300000e+04, 2.400000e+04, 2.500000e+04, 2.600000e+04,
                        2.700000e+04, 2.800000e+04, 2.900000e+04, 3.000000e+04,
                        3.100000e+04, 3.200000e+04, 3.300000e+04, 3.400000e+04,
                        3.500000e+04, 3.600000e+04, 3.700000e+04, 3.800000e+04,
                        3.900000e+04, 4.000000e+04, 4.100000e+04, 4.200000e+04,
                        4.300000e+04, 4.400000e+04, 4.500000e+04, 4.600000e+04,
                        4.700000e+04, 4.800000e+04, 4.900000e+04]


periodic_tab = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si',
                'p', 's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co',
                'ni', 'cu', 'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr']

########################################################################################
# Ion class
class Ion:
    def __init__(self, element: str, charge: int, mass: float = None, IP: float = None):
        self.element = element
        self.charge = charge
        self.mass = mass
        self.ip = IP

# Dictionary for Ion
Ions = {
    'H': Ion('H', 0, mass=0.0, IP=0.0),
    'H1+': Ion('H', 1),
    'He': Ion('He', 0),
    'He1+': Ion('He', 1,),
    'He2+': Ion('He', 2, ),
    'C': Ion('C', 0),
    'C1+': Ion('C', 1, ),
    'C2+': Ion('C', 2, ),
    'C3+': Ion('C', 3, ),
    'C4+': Ion('C', 4, ),
    'C5+': Ion('C', 5, ),
    'C6+': Ion('C', 6, ),
    'N': Ion('N', 0),
    'N1+': Ion('N', 1, ),
    'N2+': Ion('N', 2, ),
    'N3+': Ion('N', 3, ),
    'N4+': Ion('N', 4, ),
    'N5+': Ion('N', 5, ),
    'N6+': Ion('N', 6, ),
    'N7+': Ion('N', 7, ),
    'O': Ion('O', 0),
    'O1+': Ion('O', 1, ),
    'O2+': Ion('O', 2, ),
    'O3+': Ion('O', 3, ),
    'O4+': Ion('O', 4, ),
    'O5+': Ion('O', 5, ),
    'O6+': Ion('O', 6, ),
    'O7+': Ion('O', 7, ),
    'O8+': Ion('O', 8, ),
    'Ne': Ion('Ne', 0),
    'Ne1+': Ion('Ne', 1, ),
    'Ne2+': Ion('Ne', 2, ),
    'Ne3+': Ion('Ne', 3, ),
    'Ne4+': Ion('Ne', 4, ),
    'Ne5+': Ion('Ne', 5, ),
    'Ne6+': Ion('Ne', 6, ),
    'Ne7+': Ion('Ne', 7, ),
    'Ne8+': Ion('Ne', 8, ),
    'Ne9+': Ion('Ne', 9, ),
    'Ne10+': Ion('Ne', 10, ),
    'Si': Ion('Si', 0),
    'Si1+': Ion('Si', 1, ),
    'Si2+': Ion('Si', 2, ),
    'Si3+': Ion('Si', 3, ),
    'Si4+': Ion('Si', 4, ),
    'Si5+': Ion('Si', 5, ),
    'Si6+': Ion('Si', 6, ),
    'Si7+': Ion('Si', 7, ),
    'Si8+': Ion('Si', 8, ),
    'Si9+': Ion('Si', 9, ),
    'Si10+': Ion('Si', 10, ),
    'Si11+': Ion('Si', 11, ),
    'Si12+': Ion('Si', 12, ),
    'Si13+': Ion('Si', 13, ),
    'Si14+': Ion('Si', 14, ),
    'S': Ion('S', 0),
    'S1+': Ion('S', 1, ),
    'S2+': Ion('S', 2, ),
    'S3+': Ion('S', 3, ),
    'S4+': Ion('S', 4, ),
    'S5+': Ion('S', 5, ),
    'S6+': Ion('S', 6, ),
    'S7+': Ion('S', 7, ),
    'S8+': Ion('S', 8, ),
    'S9+': Ion('S', 9, ),
    'S10+': Ion('S', 10, ),
    'S11+': Ion('S', 11, ),
    'S12+': Ion('S', 12, ),
    'S13+': Ion('S', 13, ),
    'S14+': Ion('S', 14, ),
    'S15+': Ion('S', 15, ),
    'S16+': Ion('S', 16, ),
    'Fe': Ion('Fe', 0),
    'Fe1+': Ion('Fe', 1, ),
    'Fe2+': Ion('Fe', 2, ),
    'Fe3+': Ion('Fe', 3, ),
    'Fe4+': Ion('Fe', 4, ),
    'Fe5+': Ion('Fe', 5, ),
    'Fe6+': Ion('Fe', 6, ),
    'Fe7+': Ion('Fe', 7, ),
    'Fe8+': Ion('Fe', 8, ),
    'Fe9+': Ion('Fe', 9, ),
    'Fe10+': Ion('Fe', 10, ),
    'Fe11+': Ion('Fe', 11, ),
    'Fe12+': Ion('Fe', 12, ),
    'Fe13+': Ion('Fe', 13, ),
    'Fe14+': Ion('Fe', 14, ),
    'Fe15+': Ion('Fe', 15, ),
    'Fe16+': Ion('Fe', 16, ),
    'Fe17+': Ion('Fe', 17, ),
    'Fe18+': Ion('Fe', 18, ),
    'Fe19+': Ion('Fe', 19, ),
    'Fe20+': Ion('Fe', 20, ),
    'Fe21+': Ion('Fe', 21, ),
    'Fe22+': Ion('Fe', 22, ),
    'Fe23+': Ion('Fe', 23, ),
    'Fe24+': Ion('Fe', 24, ),
    'Fe25+': Ion('Fe', 25, ),
    'Fe26+': Ion('Fe', 26, ),
}

##########################################################################################################







planck = 6.6260693e-27   #erg s
planckEv = 4.13566743e-15  # ev s
hbar = planck/(2.*np.pi)
light = 29979245800.  # cm/s
q = 4.80320425e-10  # the units of charge

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

