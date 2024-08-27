import NebulaPy.src as nebula
from NebulaPy.tools import util
from pypion.ReadData import ReadData
import astropy.units as unit
import numpy as np
import time

# MIMIR: Set up paths and filenames
#silo_dir = '/mnt/massive-stars/data/mpv10-wind-test'  # Directory containing silo files
#output_file = '/net/maedoc.ap.dias.ie/maedoc/home_cr/arun/Desktop/plots/line_luminosity_OIV25.txt' # Output file for results

# RAZERBLADE: Set up paths and filenames
silo_dir = '/home/tony/Desktop/OIV_Luminosity/mpv10-wind-test'  #Directory containing silo files
output_file = '/home/tony/Desktop/OIV_Luminosity/line_luminosity_OIV25.txt'  # Output file for results

# filebase
filebase = 'e7_WRwind_d1l5n256_v1500'  # Base name of the silo files

# Set up the ion and line emission parameters
ion_pion_name = 'O3+'  # The ion of interest (Oxygen IV)
mass_O = 2.6567628e-23  # Mass of an Oxygen atom (in grams)
line = 2.589332e+05  # Emission line of interest
print(rf" calculating line luminosity of {ion_pion_name} {line} Angstrom")
# Batch the silo files according to the time instant
start_time = 0.0
finish_time = 85.0
batched_silos = util.batch_silos(
    silo_dir,
    filebase,
    start_time=None,
    finish_time=finish_time,
    time_unit='kyr',
    out_frequency=None
)

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

line_emission_obj = nebula.line_emission(ion_pion_name, verbose=True)  # Initialize the emission line calculation

# Open the output file and write the header
with open(output_file, "w") as file:
    file.write("# " + util.nebula_version()+'\n')
    file.write("# flow velocity: 1500 km/s \n")
    file.write(f"# time(kyr)"
               f" L[{line}](erg/s)\n")

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

    # Extract necessary physical parameters for the current time instant
    density = pion.get_parameter('Density', silo_instant)  # Retrieve density
    temperature = pion.get_parameter('Temperature', silo_instant)  # Retrieve temperature
    ion_massfrac = pion.get_ion_values(ion_pion_name, silo_instant)  # Retrieve ion mass fraction

    # Calculate ion number density (number of ions per unit volume)
    ion_num_density = ion_massfrac * density / mass_O

    # Retrieve the electron number density
    ne = pion.get_ne(silo_instant)

    # Calculate the line luminosity for the specific emission line
    line_emission_obj.lineluminosity_spherical(
        line=line,
        temperature=temperature,
        ne=ne,
        ns=ion_num_density,
        dV=shell_volume
    )


    # Print the calculated luminosity for the current time instant
    print(f" L_{line} = {line_emission_obj.line_emission_container['luminosity']} erg/s")

    # Append the calculated luminosity to the output file
    with open(output_file, "a") as file:
        file.write(f"{sim_time.value:.6e}"
                   f"\t{line_emission_obj.line_emission_container['luminosity']:.6e}\n"
                   )

    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")

