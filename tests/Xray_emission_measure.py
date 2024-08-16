import NebulaPy.src as nebula
from NebulaPy.tools import util
import multiprocessing as mp

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

# elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si']
elements = ['H', 'He']

xray_emission = nebula.xray(
    min_photon_energy=0.3,  # Minimum photon energy in keV
    max_photon_energy=7.0,  # Maximum photon energy in keV
    elements=elements,
    ncores=1,
    verbose=True
)

temperature = [1000, 10000]
ne = [1.5, 10]

#mp.freeze_support()  # Ensures proper support for frozen scripts
# Assuming that the method is being called within an object method
# You need to call xray_intensity on an instance of the object, e.g.,

xray_emission.xray_intensity(temperature=temperature, ne=ne)



