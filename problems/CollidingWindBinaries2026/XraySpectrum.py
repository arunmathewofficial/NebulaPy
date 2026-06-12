#import numpy as np
import NebulaPy.src as nebula
from NebulaPy.tools import util
import time
#from pypion.ReadData import ReadData
import matplotlib.pyplot as plt
import numpy as np
#import astropy.units as unit
import os
#from NebulaPy.tools import constants as const
#import pandas as pd
import ChiantiPy.tools.filters as chfilters
import matplotlib.pyplot as plt  # Plotting
from mpl_toolkits.axes_grid1 import make_axes_locatable  # For attaching colorbars to axes
from matplotlib.ticker import MultipleLocator, ScalarFormatter  # For controlling tick formatting
from NebulaPy.src.CIE import cieMode

# info: code to test emission measure calculation.
import ChiantiPy.core as ch
import ChiantiPy.tools.data as chdata
from NebulaPy.src.EmissionMeasure import emissionMeasure

import NebulaPy.src.Chianti as nebula_chainti

# constants
cm2au = 6.68459e-14  # cm to au conversion factor

# Macbook
OutputDir = '/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/XraySpectrum'  # Output image directory

# Colliding wind binaries
#Razer Blade -> Set up paths and filenames
OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'  # Output image directory
SiloDir = '/home/tony/Desktop/CWBs-2026/Silo-n128'  # Directory containing silo files
Filebase = 'wr140_NEMO_d07e13_d2l6n128'  # Base name of the silo files
start_time = 1.24e6  # in sec
finish_time = None
time_unit = 'sec'
out_frequency = None
SimulationName = "CWB"


# Bowshock
#Razer Blade -> Set up paths and filenames
#OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'  # Output image directory
#SiloDir = '/home/tony/Desktop/multi-ion-bowshock/sim-output/silo'  # Directory containing silo files
#Filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'  # Base name of the silo files
#start_time = 161  # in kyr
#finish_time = 161.5
#time_unit = 'kyr'
#out_frequency = None
#SimulationName = "Bowshock"

# Batch the silo files according to the time instant
batched_silos = util.batch_silos(
    SiloDir,
    Filebase,
    start_time=start_time,
    finish_time=finish_time,
    time_unit=time_unit,
    out_frequency=out_frequency
)

'''
key = input(" Press 'y' to continue, anything else to exit: ").strip().lower()

if key == "y":
    print(" Continuing execution...")
else:
    util.nebula_info("Resetting parameters before the next run")
    exit(0)
'''

# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# loading geometry attributes from the first silo file in the batch
# and saves them into a geometry container.
pion.load_geometry(scale='cm')
N_grid_level = pion.geometry_container['Nlevel']
mesh_edges_min = pion.geometry_container['edges_min']
mesh_edges_max = pion.geometry_container['edges_max']
N_grid = pion.geometry_container['Ngrid']
grid_volume = pion.get_grid_volumes_2D()
grid_mask = pion.geometry_container['mask']

# loading chemistry container for pion simulation data
pion.load_chemistry()
elements = pion.get_elements()

# initializing spectrum class
NebulaSpectrum = nebula.spectrum(
    min_wavelength=1.0,  # Minimum wavelength in Angstroms
    max_wavelength=20,  # Maximum wavelength in Angstroms
    min_photon_energy=None,  # Minimum photon energy in keV
    max_photon_energy=None,  # Maximum photon energy in keV
    elements=elements,
    doBremsstrahlung=True,
    doFreebound=True,
    doLine=True,
    doTwophoton=True,
    CIE=True,
    filtername=None,
    filterfactor=None,
    userGrid=True,
    gridSize=4000,
    allLines=True,
    verbose=True
)

wavelength = NebulaSpectrum.WavelengthGrid

runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()

    print(" [ SIMULATION SNAPSHOT ] " + "─" * 70)
    sim_time = pion.get_simulation_time(silo_instant, time_unit='sec')
    print(f" Step: {step}  |  Simulation time: {sim_time:.6e} s")

    # Extract temperature and electron number density
    temperature = pion.get_parameter('Temperature', silo_instant)

    temperature = np.asarray(temperature, dtype=np.float64)

    #density = pion.get_parameter('Density', silo_instant)
    ne = pion.get_ne(silo_instant)
    species_densities = pion.get_species_number_densities(silo_instant)

    NebulaSpectrum.generateSpectrum(temperature=temperature, ne=ne,
                               species_densities=species_densities,
                               grid_volume=grid_volume, grid_mask=grid_mask)


    NebulaSpectrum = NebulaSpectrum.Spectrum

    energy = 12.39841984 / wavelength
    # sort in increasing energy
    idx = np.argsort(energy)
    energy = energy[idx]

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.set_xlabel(r"Wavelength [$\AA$]", fontsize=12)
    ax.set_ylabel(r"$L_\lambda$ [erg s$^{-1}$ $\AA^{-1}$]", fontsize=12)

    ax.plot(
        wavelength,
        NebulaSpectrum,
        color="blue",
        linewidth=1.4,
        label="NEQ Spectrum"
    )

    # Logarithmic luminosity axis
    ax.set_yscale("log")

    # Publication-style ticks
    ax.minorticks_on()

    ax.tick_params(
        axis='both',
        which='major',
        direction='in',
        top=True,
        right=True,
        length=6,
        width=1.2,
        labelsize=11
    )

    ax.tick_params(
        axis='both',
        which='minor',
        direction='in',
        top=True,
        right=True,
        length=3,
        width=1.0
    )

    # Slightly thicker frame
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    ax.legend(
        loc='best',
        frameon=False,
        fontsize=10
    )

    fig.tight_layout()

    outfile = f"{OutputDir}/spectrum_final.png"
    fig.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(f"time: {sim_time:.6e}, saved snapshot {step} to {outfile}")

    '''
    #####################################################################################
    print(" ChiantiPy Calculation ############################")
    temperature = np.array([10 ** 7, 10 ** 7.0])
    ne = np.array([1.e+9, 1.e+9])
    density = np.array([1.e+9, 1.e+9])
    species_density = np.array([1.0, 1.0])
    shell_volume = np.array([1.0, 1.0])
    min_abund = 2.e-5
    ionList = ['h_2']
    spec = ch.spectrum(
        temperature,
        density[0],
        wavelength,
        ionList=ionList,
        doLines=False,
        doContinuum=True,
        # minAbund=min_abund,
        em=None,  # None will set EM =1.0 for every temperature element.
        verbose=True,
    )
    energy = 12.39841984 / wavelength
    # sort in increasing energy
    idx = np.argsort(energy)
    energy = energy[idx]
    xray_spectrum = spec.Spectrum['intensity'][0][idx]

    plt.figure()
    plt.title(f" ChiantiPy {ionList[0]} Spectrum at T = {temperature[0]:.2e} K")
    plt.plot(energy, xray_spectrum, linewidth=1, color='black', label='ChiantiPy Spectrum')
    plt.legend()
    xy = plt.axis()
    outfile = OutputDir + f"/SpectrumChiantipy{SimulationName}{ionList[0]}.png"
    plt.savefig(outfile)
    ##############################################################################################
    '''













# initial conditions ##############################################################
'''
ionlist = ['c_5']
print(f" Testing for ion {ionlist[0]}")
print(f" -------------------------------\n")

temperature = np.array([10**7, 10**7.0])
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



