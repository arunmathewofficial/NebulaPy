#import numpy as np
import NebulaPy.src as nebula
from NebulaPy.tools import util
import time
from pypion.ReadData import ReadData
import matplotlib.pyplot as plt
import numpy as np
#import astropy.units as unit
import os
#from NebulaPy.tools import constants as const
#import pandas as pd
import ChiantiPy.tools.filters as chfilters


# info: code to test emission measure calculation.
import ChiantiPy.core as ch
import ChiantiPy.tools.data as chdata


# constants
cm2au = 6.68459e-14  # cm to au conversion factor

# Macbook
OutputDir = '/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/XraySpectrum'  # Output image directory
# Razer Blade
OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'



#Razer Blade -> Set up paths and filenames
OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'  # Output image directory
SiloDir = '/home/tony/Desktop/CWBs-2026/Silo-n128'  # Directory containing silo files
Filebase = 'wr140_NEMO_d07e13_d2l6n128'  # Base name of the silo files
start_time = 1.24e6  # in sec
finish_time = None
time_unit = 'sec'
out_frequency = None


# Batch the silo files according to the time instant
batched_silos = util.batch_silos(
    SiloDir,
    Filebase,
    start_time=start_time,
    finish_time=finish_time,
    time_unit=time_unit,
    out_frequency=out_frequency
)

key = input(" Press 'y' to continue, anything else to exit: ").strip().lower()

if key == "y":
    print(" Continuing execution...")
else:
    util.nebula_info("Resetting parameters before the next run")
    exit(0)



# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# load chemistry
pion.load_chemistry()

# Calculates and stores geometric grid parameters.
# For example, in a spherical geometry, it extracts radius and shell volumes
# from the first silo file in the batch and saves them into a geometry container.
pion.load_geometry(scale='cm')
N_grid_level = pion.geometry_container['Nlevel']
#mesh_edges_min = pion.geometry_container['edges_min'] * cm2au
#mesh_edges_max = pion.geometry_container['edges_max'] * cm2au
mesh_edges_min = pion.geometry_container['edges_min']
mesh_edges_max = pion.geometry_container['edges_max']
N_grid = pion.geometry_container['Ngrid']
cell_volume = pion.get_2D_cell_volumes()



EM = nebula.emission_measure(Tmin=100, Tmax=1.e9, Nbins=100)

runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()

    print(f" ---------------------------")
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # Extract temperature and electron number density
    temperature = pion.get_parameter('Temperature', silo_instant)
    density = pion.get_parameter('Density', silo_instant)
    ne = pion.get_ne(silo_instant)
    grid_mask = pion.geometry_container['mask']


    EM = nebula.emission_measure(temperature, Tmin=100, Tmax=1.e9, Nbins=100)

    #em = EM.DEM2D(density=density, temperature=temperature,
    #              ne=ne,
    #              mask=grid_mask,
    #              ngrid=N_grid,
    #              mesh_edges_min=mesh_edges_min,
    #              mesh_edges_max=mesh_edges_max,
    #              volume = cell_volume,
    #              temp_bin=temperature_bin,
    #              hw=0.05)

    #dem = em['dem_bin']

    # Plot the differential emission measure for the current time instant
    plt.figure()
    plt.title(f"Differential Emission Measure")
    #plt.plot(temperature_bin, np.log10(dem))
    Filename = f"{Filebase}_DEM_{sim_time.value:.2e}kyr.png"
    OutImageFile = os.path.join(OutputDir, Filename)
    plt.savefig(OutImageFile, bbox_inches="tight", dpi=300)
    plt.close()

    print(f" time: {sim_time:.6e}, Saved snapshot {step} to {Filename}")

    dt = time.time() - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")



# CIE TEST ####################################################################
'''
# fixing the following issues
# info: CIE calculation match for Nebulapy.
database = nebula.database(verbose=True)
database.load_cie()
Temperature = np.logspace(np.log10(200), np.log10(1e9), 100)
plt.figure()
plt.title("Fe Ionisation Fractions (CIE)")
# Loop over all Fe ions (Fe I to Fe XXVI → 1 to 26)
for i in range(1, 28):
    chianti_ion = f'fe_{i}'

    CIE_data = database.get_cie_fraction(chianti_ion, Temperature)

    plt.plot(np.log10(Temperature), np.log10(CIE_data),
             linewidth=2,
             label=chianti_ion.upper())
plt.ylim(-2, 0)  # Adjust y-axis limits for better visibility
plt.xlabel("log10(Temperature [K])")
plt.ylabel("Ion Fraction")
plt.legend(ncol=2, fontsize=8)  # better layout for many ions

outfile = OutputDir + "/cie_fe_all.png"
plt.savefig(outfile, dpi=300)
plt.close()
'''

'''
# initial conditions ##############################################################

ionlist = ['fe_25']
print(f" Testing for ion {ionlist[0]}")
print(f" -------------------------------\n")

temperature = np.array([10**7.0, 10**7.0])
ne = np.array([1.e+9, 1.e+9])
density = np.array([1.e+9, 1.e+9])
species_density = np.array([1.0, 1.0])
shell_volume = np.array([1.0, 1.0])

print("NEBULAPY #########################################################################")


# info: nebulapy calculation for the same conditions.

elements = ["H", "He", "Fe"]
elemental_abundances = {"H": 0.70, "He": 0.29879, "Fe": 1.21e-03}

spectrum = nebula.spectrum(
    min_wavelength=1.0,  # Minimum wavelength in Angstroms
    max_wavelength=13,  # Maximum wavelength in Angstroms
    min_photon_energy=None,  # Minimum photon energy in keV
    max_photon_energy=None,  # Maximum photon energy in keV
    elements=elements,
    elemental_abundances=elemental_abundances,
    CIE=True,
    bremsstrahlung=True,
    freebound=True,
    line=True,
    twophoton=True,
    filtername=None,
    filterfactor=None,
    user_grid=True,
    grid_size=4000,
    allLines=True,
    verbose=True
)

print("\n")
print(spectrum.spectrum_container.keys())
wavelength = spectrum.spectrum_container['wavelength_grid']

# Calculate number density
# info: number denisty calculation is just for testing the module
mass_fraction = elemental_abundances['Fe']
# Fraction of number denisty of Fe to H, from chianti database
# sun_photospheric_2015_scott abundance.
A_Fe = 2.96e-5
print(f" Fe Abundance: {A_Fe:.2e} cm^-3")
# get the CIE ion fraction for the given
# ion and temperature
cie = nebula.cieMode(verbose=True)
cie.load_cie()
ionfrac = cie.get_cie_fraction(ionlist[0], temperature)
print(f"CIE {ionlist[0]} ion fraction: {ionfrac}")
A_ion = A_Fe * ionfrac
print(f"Ion Abuandance {ionlist[0]}: {A_ion}")

# nebulapy calculation
spectrum.compute_spectrum_1D(
    temperature=temperature,
    ne=ne,
    species_density=species_density,
    shell_volume=shell_volume
)

total_emission = spectrum.spectrum_container['spectrum']



# scale by the number density of the ion and emission measure
# EM is set to unity here
em = 1.0
total_emission = total_emission * A_ion[:, np.newaxis]


energy = 12.39841984 / wavelength
# sort in increasing energy
idx = np.argsort(energy)
energy = energy[idx]

total_emission = total_emission[:, idx]

plt.figure()
plt.title(f" NebulaPy {ionlist[0]} spectrum at T = {temperature[0]:.2e} K")
plt.plot(energy, total_emission[0], linewidth=1, color='black', label='NebulaPy Spectrum')
plt.legend()
xy = plt.axis()
outfile = OutputDir + f"/spectrum_nebulapy.png"
plt.savefig(outfile)

print("\n")
print("CHIANTIPY#########################################################################")


# info: original chianti calculation using spectrum method.
#################################################################################
# this calculation use sun_photospheric_2015_scott abundance


min_abund = 2.e-5
spec = ch.spectrum(
    temperature,
    density[0],
    wavelength,
    ionList=ionlist,
    doLines=True,
    doContinuum=True,
    #minAbund=min_abund,
    em=None,# None will set EM =1.0 for every temperature element.
    verbose=True,
)

energy = 12.39841984 / wavelength
# sort in increasing energy
idx = np.argsort(energy)
energy = energy[idx]
xray_spectrum = spec.Spectrum['intensity'][0][idx]
plt.figure()
plt.title(f" ChiantiPy {ionlist[0]} Spectrum at T = {temperature[0]:.2e} K")
plt.plot(energy, xray_spectrum, linewidth=1, color='black', label='ChiantiPy Spectrum')
plt.legend()
xy = plt.axis()
outfile = OutputDir + f"/spectrum_chiantipy.png"
plt.savefig(outfile)
#################################################################################

# info: So far what is acheived
# info: 1)  the cie ionisation fraction PION obtain is slightly different
#  from what chinati has. This indroduce a small difference in the spectrum.
# info: 2) checked both free-free and free-bound emission, they are in good
#  agreement between NebulaPy and ChiantiPy.
# info: 3) while calculating spectrum for fe_25, I exclude fe_27 which is calculated
#  while calculating the spectrum for fe_25. However, including fe_27 does not
#  change the spectrum.
# todo: 4) the line emission is different between NebulaPy and ChiantiPy,
#  which need to be rectified next.

'''

















'''
#info: Using differential emission measures (DEM)
#################################################################################
#flare.dem → strong high-ionization X-ray/EUV lines (Fe XXIII–XXV)
#coronal_hole.dem → cooler low-density plasma
#active_region.dem → enhanced Fe XII–Fe XVI emission
#prominence.dem → cooler chromospheric/transition-region plasma

import ChiantiPy.tools.io as chio
demDir = os.path.join(os.environ['XUVTOP'], 'dem')
demList = os.listdir(demDir)

for idx, demFile in enumerate(demList):
    print('%i  %s'%(idx, demFile))

demFile_no = 0
flDict = chio.demRead(demList[demFile_no])
flTemp = flDict['temperature'][15:]
flDens = flDict['density'][15:]
flDem = flDict['dem'][15:]
flEm = flDict['em'][15:]
plt.figure()
plt.title(f" DEM {demList[demFile_no]}")
plt.plot(flTemp, flDem, linewidth=1, color='black', label='DEM')
plt.legend()
xy = plt.axis()
outfile = OutputDir + f"/dem.png"
plt.savefig(outfile)

'''


# sun_photospheric_2015_scottabundances from chianti database
''' 
 1  12.00  H
 2  10.93  He
 3   1.05  Li
 4   1.38  Be
 5   2.70  B
 6   8.43  C
 7   7.83  N
 8   8.69  O
 9   4.40  F
10   7.93  Ne
11   6.21  Na
12   7.59  Mg
13   6.43  Al
14   7.51  Si
15   5.41  P
16   7.12  S
17   5.50  Cl
18   6.40  Ar
19   5.04  K
20   6.32  Ca
21   3.16  Sc
22   4.93  Ti
23   3.89  V
24   5.62  Cr
25   5.42  Mn
26   7.47  Fe
27   4.93  Co
28   6.20  Ni
29   4.18  Cu
30   4.56  Zn


Convert log abundance to number abundance
Fe: 7.47
A(Fe)=10^{7.47−12}=10^{−4.53} = 2.95×10^{-5}
'''