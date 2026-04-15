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


# info: code to test emission measure calculation.
import ChiantiPy.core as ch
import ChiantiPy.tools.data as chdata


# constants
cm2au = 6.68459e-14  # cm to au conversion factor

# Macbook
OutputDir = '/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/XraySpectrum'  # Output image directory


'''
#Razer Blade -> Set up paths and filenames
OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'  # Output image directory
SiloDir = '/home/tony/Desktop/CWBs-2026/Silo-n128'  # Directory containing silo files
Filebase = 'wr140_NEMO_d07e13_d2l6n128'  # Base name of the silo files
start_time = None #1.24e6  # in sec
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


EM = nebula.emission_measure()


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


    temperature_bin = np.arange(2.0, 8.5, 0.02)

    em = EM.DEM2D(density=density, temperature=temperature,
                  ne=ne,
                  mask=grid_mask,
                  ngrid=N_grid,
                  mesh_edges_min=mesh_edges_min,
                  mesh_edges_max=mesh_edges_max,
                  volume = cell_volume,
                  temp_bin=temperature_bin,
                  hw=0.05)

    dem = em['dem_bin']

    # Plot the differential emission measure for the current time instant
    plt.figure()
    plt.title(f"Differential Emission Measure")
    plt.plot(temperature_bin, np.log10(dem))
    Filename = f"{Filebase}_DEM_{sim_time.value:.2e}kyr.png"
    OutImageFile = os.path.join(OutputDir, Filename)
    plt.savefig(OutImageFile, bbox_inches="tight", dpi=300)
    plt.close()

    print(f" time: {sim_time:.6e}, Saved snapshot {step} to {Filename}")

    dt = time.time() - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")

'''


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




temperature = [1.e+7, 1.e+7]
density = [1.e+9]
ne = [1.e+9, 1.e+9]
chianti_ion = 'fe_25'
em = 1.0


# info: chianti calculation using sub modules in ChiantiPy.
'''
i = ch.ion(chianti_ion, temperature, density, em=em)
c = ch.continuum(chianti_ion, temperature=temperature, em=em)
wvl = np.arange(2.6e3, 2.9e3, 0.1)
c.freeFree(wvl, includeAbund=False, includeIoneq=False)
c.freeBound(wvl, includeAbund=False, includeIoneq=False)
c.freeFree(wvl)
c.freeBound(wvl)
i.twoPhoton(wvl)
i.spectrum(wvl, filter=(chfilters.gaussian, 0.1), allLines=True)
ff_emission_chianti = c.FreeFree['intensity']
fb_emission_chianti = c.FreeBound['intensity']
#two_photon_emission_chianti = i.TwoPhoton['intensity']
line_emission_chianti = i.Spectrum['intensity']
continuum_emission_chianti = ff_emission_chianti + fb_emission_chianti
#total_emission_chianti = fb_emission_chianti + ff_emission_chianti + two_photon_emission_chianti + line_emission_chianti
total_emission_chianti = fb_emission_chianti + ff_emission_chianti + line_emission_chianti

plt.figure()
plt.title(f" {chianti_ion} Emission Spectrum at T = {temperature[0]:.2e} K")

# chianti plots
plt.plot(wvl, total_emission_chianti[0], linewidth=1, color='black', label='ChiantiPy Total')
#plt.plot(wvl, two_photon_emission_chianti[0], linewidth=1, color='green', linestyle='dashed', label='ChiantiPy 2-Photon')
plt.plot(wvl, ff_emission_chianti[0], linewidth=1, color='red', linestyle='dashed', label='ChiantiPy Free-Free')
plt.plot(wvl, fb_emission_chianti[0], linewidth=1, color='darkorange', linestyle='dotted', label='ChiantiPy Free-Bound')
plt.plot(wvl, line_emission_chianti[0], label=f'ChiantiPy Line')
energy = 12.39841984 / wvl
#print(energy)
# sort in increasing energy
idx = np.argsort(energy)
energy_sorted = energy[idx]
line_sorted = line_emission_chianti[0][idx]
#plt.plot(energy_sorted, total_emission_chianti[0][idx], label=f'ChiantiPy Total')
#plt.plot(energy_sorted, line_sorted, label=f'ChiantiPy Line')
#plt.plot(energy_sorted, two_photon_emission_chianti[0][idx], linewidth=1, color='green', linestyle='dashed', label='ChiantiPy 2-Photon')
#plt.plot(energy_sorted, ff_emission_chianti[0][idx], linewidth=1, color='red', linestyle='dashed', label='ChiantiPy Free-Free')
#plt.plot(energy_sorted, fb_emission_chianti[0][idx], linewidth=1, color='darkorange', linestyle='dotted', label='ChiantiPy Free-Bound')

#plt.ylim(1.e-4, 5)
#plt.yscale('log')
plt.legend()
xy = plt.axis()
outfile = OutputDir + "/xray_spectrum.png"
plt.savefig(outfile)
'''



spectrum = nebula.spectrum(
    min_wavelength=1.0,  # Minimum wavelength in Angstroms
    max_wavelength=13,  # Maximum wavelength in Angstroms
    min_photon_energy=None,  # Minimum photon energy in keV
    max_photon_energy=None,  # Maximum photon energy in keV
    N_point=2000,
    bremsstrahlung=True,
    freebound=True,
    line=True,
    twophoton=False,
    filtername=None,
    filterfactor=None,
    allLines=True,
    verbose=True
)


elements = ["H", "He", "Fe"]
# we use fe abundance =   2.95e-05 todo: is this mass fraction?
elemental_abundances = {"H": 0.70, "He": 0.2999705, "Fe": 2.95e-05}
spectrum.build_species_attributes(elements=elements, elemental_abundances=elemental_abundances)

# nebulapy calculation
wvl = spectrum.wavelength
spectrum.compute_emission(temperature=temperature, density=density, ne=ne)
ff_emission = spectrum.intensity['bremsstrahlung']
fb_emission = spectrum.intensity['freebound']
# todo: line emission rate do not match, need to fix it.
line_emission = spectrum.intensity['line']

continuum_emission = ff_emission + fb_emission
total_emission = continuum_emission + line_emission

abuandance = elemental_abundances['Fe']
# get the CIE ion fraction for the given ion and temperature
database = nebula.database(verbose=True)
database.load_cie()
ionfrac = database.get_cie_fraction(chianti_ion, temperature)
print(f"CIE {chianti_ion} ion fraction: {ionfrac}")


# line emission for the same conditions
#resolving_power = len(wvl)/1000
#line = ch.ion(chianti_ion, temperature=temperature, eDensity=ne, em=1.0, verbose=True)
#line.spectrum(wvl, filter=(chfilters.gaussian, resolving_power), allLines=True)
# my calculation for line emission
#my_line_emission = spectrum.intensity['line']


energy = 12.39841984 / wvl
# sort in increasing energy
idx = np.argsort(energy)
energy = energy[idx]

plt.figure()
plt.title(f" {chianti_ion} x-ray spectrum at T = {temperature[0]:.2e} K")
plt.plot(energy, total_emission[0][idx], linewidth=1, color='black', label='ChiantiPy Spectrum')
plt.legend()
xy = plt.axis()
outfile = OutputDir + "/xray_nebulapy.png"
plt.savefig(outfile)


# info: original chianti calculation using spectrum method.
'''
ionlist= ['fe_25']
temperature = np.array([1e+7, 1e+7])
density = 1e9
wavelength = np.linspace(1, 13, 2000)
min_abund = 2.e-5
spec = ch.spectrum(
    temperature,
    density,
    wavelength,
    ionList=ionlist,
    doLines=True,
    doContinuum=True,
    #minAbund=min_abund,
    em=1.e27,
    verbose=True,
)


energy = 12.39841984 / wavelength
# sort in increasing energy
idx = np.argsort(energy)
energy = energy[idx]
xray_spectrum = spec.Spectrum['intensity'][0][idx]
plt.figure()
plt.title(f" {ionlist[0]} x-ray spectrum at T = {temperature[0]:.2e} K")
plt.plot(energy, xray_spectrum, linewidth=1, color='black', label='ChiantiPy Spectrum')
plt.legend()
xy = plt.axis()
outfile = OutputDir + "/xray_chainti_test.png"
plt.savefig(outfile)
'''