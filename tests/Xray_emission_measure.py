import NebulaPy.src as nebula
from NebulaPy.tools import util
import time
import matplotlib.pyplot as plt
import numpy as np


'''
# Set up paths and filenames
silo_dir = '/home/tony/Desktop/NebulaPy/tests/wind-wind-jm'  # Directory containing silo files
filebase = 'e7_WRwind_d1l5n256_v0750'  # Base name of the silo files
output_file = '/home/tony/Desktop/NebulaPy/tests/xray.txt'  # Output file for results

# Batch the silo files according to the time instant
batched_silos = util.batch_silos(silo_dir, filebase)

# Initialize the Pion class from NebulaPy, which handles the simulation data
nebula_pion = nebula.pion(batched_silos, verbose=True)

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
nebula_pion.get_chemistry()

elements = nebula_pion.chemistry_container['tracer_elements']
'''

elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si', 'S', 'Fe']
solar_abundance = []
ionisation_equilibrium = []

xray_emission = nebula.xray(
    min_photon_energy=0.04,  # Minimum photon energy in keV
    max_photon_energy=0.062,  # Maximum photon energy in keV
    energy_point_count=1000,
    elements=elements,
    verbose=True
)

temperature = [2e+6]
ne = [1e+9]
em = [1e+27]


#mp.freeze_support()  # Ensures proper support for frozen scripts
# Assuming that the method is being called within an object method
# You need to call xray_intensity on an instance of the object, e.g.,
#if __name__ == '__main__':
#    mp.freeze_support()


runtime = 0.0
start_time = time.time()
spectrum = xray_emission.xray_intensity(
    temperature=temperature,
    ne=ne,
    elementalabunds=solar_abundance,
    ionfractions=ionisation_equilibrium,
    emissionmeasure=em,
    freefree=True, freebound=True, lines=True, twophoton=False,
    multiprocessing=True,
    ncores=12)
finish_time = time.time()  # Record the finish time
# Calculate the time spent on the current step
dt = finish_time - start_time
# Update the runtime with the time spent on the current step
#runtime += dt
print(f" runtime: {runtime:.4e} | dt: {dt:.4e} s")

generated_wvl_array = xray_emission.xray_containter['wvl_array']

# Create a plot
plt.figure(figsize=(8, 6))  # Set the figure size
plt.plot(generated_wvl_array, spectrum[0], linestyle='-', color='b', label=f'T = {temperature[0]:.2e} (K)')
plt.xlabel(r'$\lambda \, (\AA)$', fontsize=14)
plt.ylabel('Spectrum', fontsize=14)
plt.legend(fontsize=14, frameon=False)
plt.show()


