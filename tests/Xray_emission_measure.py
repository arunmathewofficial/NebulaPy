import NebulaPy.src as nebula
from NebulaPy.tools import util
import time
from pypion.ReadData import ReadData
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as unit
import os


# Set up paths and filenames
silo_dir = '/home/tony/Desktop/Equi_NonEqui/nonequi-equi/'  # Directory containing silo files
filebase = 'e7_WRwind_d1l5n256_v1500'  # Base name of the silo files
output_path = '/home/tony/Desktop/Equi_NonEqui/nonequi-equi/'  # Output file for results

# Batch the silo files according to the time instant
batched_silos = util.batch_silos(silo_dir, filebase)

# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
pion.load_chemistry()

# Initialize spherical grid parameters (e.g., radius, shell volumes)
# This sets up the grid using the first silo file in the batch
pion.spherical_grid(batched_silos[0])

# Retrieve the radius and shell volumes from the geometry container
radius = pion.geometry_container['radius']
shell_volume = pion.geometry_container['shell_volumes']
elements = pion.get_elements()

xray_emission = nebula.xray(
    min_photon_energy=0.5,  # Minimum photon energy in keV
    max_photon_energy=10.0,  # Maximum photon energy in keV
    energy_point_count=1000,
    elements=elements,
    Tmin=1e+5, Tmax=1e+9,
    bremsstrahlung=True,
    freebound=False,
    lines=False,
    twophoton=False,
    multiprocessing=True,
    ncores=12,
    verbose=True
)

runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time

    # Read the data from the current silo file
    dataio = ReadData(silo_instant)
    basic = dataio.get_1Darray('Density')  # Retrieve basic simulation data, such as density
    sim_time = (basic['sim_time'] * unit.s).to(unit.kyr)  # Convert simulation time to kiloyears
    dataio.close()  # Close the data file

    # Print the current time instant
    print(f" ---------------------------")
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    temperature = pion.get_parameter('Temperature', silo_instant)
    density = pion.get_parameter('Density', silo_instant)
    ne = pion.get_ne(silo_instant)
    elemental_mass_fraction = pion.get_elemental_mass_frac(silo_instant)
    ion_fractions = pion.get_tracer_values(silo_instant)

    spectrum = xray_emission.xray_intensity(
        temperature=temperature,
        density=density, ne=ne,
        elemental_abundances=elemental_mass_fraction,
        ion_fractions=ion_fractions,
        shell_volume=shell_volume
    )

    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")


    # making plots
    generated_wvl_array = xray_emission.xray_containter['wvl_array']
    filtered_temperature = xray_emission.xray_containter['temperature']
    for i in range(len(spectrum)):
        plt.figure(figsize=(8, 6))  # Set the figure size
        # Format and replace '.' with 'p' to avoid issues in filenames
        temp_str = f"{filtered_temperature[i]:.2e}"
        out_filename = filebase + f"_time{sim_time.value:.2f}_T{temp_str}.png"
        out_file = os.path.join(output_path, out_filename)

        # Plot the spectrum with the corresponding temperature
        plt.plot(generated_wvl_array, spectrum[i], linestyle='-', color='b', label=f'T = {filtered_temperature[i]:.2e} (K)')

        plt.xlabel(r'$\lambda \, (\AA)$', fontsize=14)
        plt.ylabel('Spectrum', fontsize=14)
        plt.legend(fontsize=14, frameon=False)

        plt.savefig(out_file)  # Save as a PNG file
        plt.close()  # Close the plot to free memory






    '''
    # differential emission measure
    DEM = pion.DEM(
        temperature=temperature,
        ne=ne,
        shellvolume=shell_volume,
        Tmin=1e+5,
        Tmax=1e+9,
        Nbins=100
    )
    '''








