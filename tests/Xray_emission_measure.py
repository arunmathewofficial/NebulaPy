import NebulaPy.src as nebula
from NebulaPy.tools import util

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

elements = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Si']

xray_emission = nebula.xray(0.3, 7.0, elements, verbose=True)