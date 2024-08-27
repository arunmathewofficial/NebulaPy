import NebulaPy.src as nebula
from NebulaPy.tools import util
import time
import matplotlib.pyplot as plt
import numpy as np



# Set up paths and filenames
silo_dir = '/home/tony/Desktop/OIV_Luminosity/wind-wind-equi/'  # Directory containing silo files
filebase = 'e7_WRwind_d1l5n256_v1500'  # Base name of the silo files
output_file = '/home/tony/Desktop/OIV_Luminosity/wind-wind-equi/'  # Output file for results

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

'''
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time

    # Read the data from the current silo file
    #dataio = ReadData(silo_instant)
    #basic = dataio.get_1Darray('Density')  # Retrieve basic simulation data, such as density
    #sim_time = (basic['sim_time'] * unit.s).to(unit.kyr)  # Convert simulation time to kiloyears
    #dataio.close()  # Close the data file

    temperature = pion.get_parameter('Temperature', silo_instant)
    ne = pion.get_ne(silo_instant)

    dem = pion.differential_emission_measure(
        temperature=temperature,
        ne=ne,
        shellvolume=shell_volume,
        Tmin=1.e+6,
        Tmax=1e+9,
        Nbins=100
    )

    # Create a plot
    plt.figure(figsize=(8, 6))  # Set the figure size
    plt.plot(np.log10(dem['Tb']), np.log10(dem['DEM']), linestyle='-', color='b')
    plt.ylabel(r'$\rm log(DEM)$ (cm$^{-3}$)', fontsize=14)
    plt.xlabel(r'$\rm log(T_b) (K)$', fontsize=14)
    plt.show()



    elemental_massfrac = pion.get_elemental_mass_frac(silo_instant)
    tracer_values = pion.get_tracer_values(silo_instant)





'''
#elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si', 'S', 'Fe']
elements = ["H", "He"]
solar_abundance = 1.0
ionisation_equilibrium = 1.0

xray_emission = nebula.xray(
    min_photon_energy=0.5,  # Minimum photon energy in keV
    max_photon_energy=10.0,  # Maximum photon energy in keV
    energy_point_count=1000,
    elements=elements,
    verbose=True
)

temperature = [2e+7]
ne = [1e+2]
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
    elemental_abundances=solar_abundance,
    ion_fractions=ionisation_equilibrium,
    emission_measure=em,
    bremsstrahlung=True, freebound=True,
    lines=True, twophoton=False,
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

