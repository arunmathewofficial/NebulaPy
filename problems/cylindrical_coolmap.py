from pypion.ReadData import ReadData
import astropy.units as unit
import numpy as np
import time

output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO_images'  # Output image directory

from NebulaPy.tools import util

# Set up paths and filenames
silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO'  # Directory containing silo files
filebase = 'BN_grad_d2l4n128'  # Base name of the silo files

# Batch the silo files according to the time instant
start_time = 0.0
finish_time = 85.0
batched_silos = util.batch_silos(
    silo_dir,
    filebase,
    start_time=None,
    finish_time=None,
    time_unit=None,
    out_frequency=None
)

import NebulaPy.src as nebula


# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# Calculates and stores geometric grid parameters.
# For example, in a spherical geometry, it extracts radius and shell volumes
# from the first silo file in the batch and saves them into a geometry container.
pion.load_geometry(batched_silos[0])

print(pion.geometry_container)

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
pion.load_chemistry()
print(pion.chemistry_container)




cooling = nebula.cooling(
    database='/home/tony/Desktop/NebulaPy/NebulaPy-DB',
    pion_ion='O3+',
    verbose=True
)

temperature = [[5000, 6000], [200, 100]]
ne = [[100, 100], [100, 100]]

#temperature = [5000, 6000]
#ne = [100, 100]

# temperature and ne array can be 1D, 2D, 3D
cooling_map = cooling.generate_cooling_map(
    temperature=temperature,
    ne=ne
)

print(cooling_map)

del cooling