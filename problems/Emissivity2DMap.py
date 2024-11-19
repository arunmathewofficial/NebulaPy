import NebulaPy.src as nebula
from NebulaPy.tools import util
from pypion.ReadData import ReadData
import astropy.units as unit
import numpy as np
import time
import warnings
# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")


output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO_images'  # Output image directory
# Set up paths and filenames
silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO'  # Directory containing silo files
filebase = 'BN_grad_d2l4n128'  # Base name of the silo files

# Set up the ion and line emission parameters
pion_ion = 'O2+'  # The ion of interest (Oxygen IV)
lines = [4960.295, 5008.240]
print(rf" calculating line luminosity of {pion_ion} lines: {lines} Angstrom")

# Batch the silo files according to the time instant
batched_silos = util.batch_silos(
    silo_dir,
    filebase,
    start_time=None,
    finish_time=None,
    time_unit=None,
    out_frequency=None
)

# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
pion.load_chemistry()

# load geometry associated with the simulation
pion.load_geometry(batched_silos[0])


runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time


    print(f" ---------------------------")
    # Print the current simulation time instant
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # Extract necessary physical parameters
    temperature = pion.get_parameter('Temperature', silo_instant)  # Retrieve temperature
    ne = pion.get_ne(silo_instant)

    # Ensure the full array is printed
    np.set_printoptions(threshold=np.inf)
    print(ne[0][0])

    temperature = [[5000, 6000], [6000, 5000]]
    ne = [[1, 1], [1, 1]]

    chianti = nebula.chianti(temperature=temperature,
                             ne=ne,
                             pion_ion=pion_ion,
                             verbose=True)

    chianti.get_emissivity()

    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")

