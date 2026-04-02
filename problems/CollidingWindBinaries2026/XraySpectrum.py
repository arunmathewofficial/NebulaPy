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



# info these need to be removed once the code is tested
import ChiantiPy.core as ch
import ChiantiPy.tools.data as chdata


# constants
cm2au = 6.68459e-14  # cm to au conversion factor


#Razer Blade -> Set up paths and filenames
OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'  # Output image directory
SiloDir = '//home/tony/Desktop/CWBs-2026/Silo-n128'  # Directory containing silo files
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


    temperature_bin = np.arange(3.0, 8.5, 0.02)

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
    Filename = f"{Filebase}_FeNumDen_{sim_time.value:.2f}kyr.png"
    OutImageFile = os.path.join(OutputDir, Filename)
    plt.savefig(OutImageFile, bbox_inches="tight", dpi=300)
    plt.close()

    print(f" time: {sim_time:.6e}, Saved snapshot {step} to {Filename}")

    dt = time.time() - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")









# info: code to test the spectrum calculation for a single temperature
#  and density point, which can be used to compare with ChiantiPy's
#  free-free emission calculation for the same conditions.


'''
spectrum = nebula.spectrum(
    min_wavelength=1.0,  # Minimum wavelength in Angstroms
    max_wavelength=100.0,  # Maximum wavelength in Angstroms
    min_photon_energy=None,  # Minimum photon energy in keV
    max_photon_energy=None,  # Maximum photon energy in keV
    N_point=4000,
    bremsstrahlung=True,
    freebound=False,
    lines=True,
    twophoton=False,
    filtername=None,
    filterfactor=None,
    allLines=True,
    verbose=True
)

elements = ["H", "He", "O"]
elemental_abundances = {"H": 0.70, "He": 0.28, "O": 0.02}
spectrum.build_species_attributes(elements=elements,
                                  elemental_abundances=elemental_abundances)
                                  
temperature = [3.e+7]
density = [1.e+9]
ne = [1.e+9, 1.e+9, 1.e+9]

# my calculation
wvl = spectrum.wavelength
spectrum.compute_emission(temperature=temperature, density=density, ne=ne)
emission_my = spectrum.intensity

# chianti calculation
c = ch.continuum('o_8', temperature=temperature, em=1.e+27,)
c.freeFree(wvl, includeAbund=False, includeIoneq=False)
emission_chianti = c.FreeFree['intensity']

plt.figure()
plt.title(f"OVIII Free–Free Emission Spectrum at T = {temperature[0]:.2e} K")
plt.plot(wvl, emission_my, linewidth=2, color='black', label='NebulaPy Free-Free')
plt.plot(wvl, emission_chianti, linewidth=2, color='red', label='ChiantiPy Free-Free', linestyle='dotted')
xy = plt.axis()
outfile = OutDir + "/xray_spectrum.png"
plt.savefig(outfile)

'''







