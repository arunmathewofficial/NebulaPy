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




database = nebula.database(verbose=True)
database.load_cie()


# fixing the following issues
# info: 1. free-free calculation match for Nebulapy and ChiantiPy.
# info: 2. free-bound calculation match for Nebulapy and ChiantiPy.
# info: 3. bug fixed in free-free, now returns 2D array of shape (N_temperature, N_wavelength).
# todo: 4. line emission rate do not match, need to fix it.
# todo: 5. add two-photon emission calculation and comparison.
# todo: 6. make sure all the emission calculations return same physical quantity and units.

'''

spectrum = nebula.spectrum(
    min_wavelength=200,  # Minimum wavelength in Angstroms
    max_wavelength=300,  # Maximum wavelength in Angstroms
    min_photon_energy=None,  # Minimum photon energy in keV
    max_photon_energy=None,  # Maximum photon energy in keV
    N_point=400,
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
elemental_abundances = {"H": 0.70, "He": 0.28, "Fe": 0.02}
spectrum.build_species_attributes(elements=elements, elemental_abundances=elemental_abundances)

temperature = [2.e+6, 3.e+7]
density = [1.e+9]
ne = [1.e+9, 1.e+9]
chianti_ion = 'fe_14'

# my calculation
wvl = spectrum.wavelength
spectrum.compute_emission(temperature=temperature, density=density, ne=ne)
my_ff_emission = spectrum.intensity['bremsstrahlung']
my_fb_emission = spectrum.intensity['freebound']

# chianti calculation
c = ch.continuum(chianti_ion, temperature=temperature, em=None)
c.freeFree(wvl, includeAbund=False, includeIoneq=False)
c.freeBound(wvl, includeAbund=False, includeIoneq=False)
ff_emission_chianti = c.FreeFree['intensity']
fb_emission_chianti = c.FreeBound['intensity']

# line emission for the same conditions
resolving_power = len(wvl)/1000
line = ch.ion(chianti_ion, temperature=temperature, eDensity=ne, em=1.0, verbose=True)
line.spectrum(wvl, filter=(chfilters.gaussian, resolving_power), allLines=True)
# my calculation for line emission
my_line_emission = spectrum.intensity['line']

plt.figure()
plt.title(f" {chianti_ion} Emission Spectrum at T = {temperature[0]:.2e} K")


plt.plot(wvl, my_ff_emission[0],linewidth=2, color='black', label='NebulaPy Free-Free')
plt.plot(wvl, my_fb_emission[0], linewidth=2, color='blue', label='NebulaPy Free-Free')
plt.plot(wvl, ff_emission_chianti[0], linewidth=2, color='red', linestyle='dashed', label='ChiantiPy Free-Free')
plt.plot(wvl, fb_emission_chianti[0], linewidth=2, color='darkorange', linestyle='dotted', label='ChiantiPy Free-Bound')

#plt.plot(wvl, line.Spectrum['intensity'][0], label=f'ChiantiPy Line {}')
#plt.plot(wvl, my_line_emission[0], label='NebulaPy Line')
#plt.ylim(0, 400)
plt.legend()
xy = plt.axis()
outfile = OutputDir + "/xray_spectrum.png"
plt.savefig(outfile)


'''