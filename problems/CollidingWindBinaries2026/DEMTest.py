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


#Razer Blade -> Set up paths and filenames
OutputDir = '/home/tony/Desktop/CWBs-2026/Postprocessing/X-raySpectrum'  # Output image directory
SiloDir = '/home/tony/Desktop/CWBs-2026/Silo-n128'  # Directory containing silo files
SiloDir = "/home/tony/Desktop/multi-ion-bowshock/high-res-silo-200kyr"
Filebase = 'wr140_NEMO_d07e13_d2l6n128'  # Base name of the silo files
Filebase = 'Ostar_mhd-nemo-dep_d2n0384l3'
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

