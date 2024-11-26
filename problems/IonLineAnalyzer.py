"""
Ion Line Analyzer
Description: Identifies the important spectral lines for the ion in the silo file.
Author: Arun Mathew
Date: 22 Nov 2024
"""

import numpy as np
import os
import time
from NebulaPy.tools import util
from NebulaPy.tools import constants as const
import NebulaPy.src as nebula
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import pandas as pd
import warnings
# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")

# MIMIR
#output_dir = '/mnt/massive-stars/data/arun_simulations/simple-bowshock-coolmap'  # Output image directory
# Set up paths and filenames
#silo_dir = '/mnt/massive-stars/data/nemo/simple-bowshock'  # Directory containing silo files
#filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'  # Base name of the silo files
#database_path = '/net/maedoc.ap.dias.ie/maedoc/home_cr/arun/Desktop/NebulaPy/NebulaPy-DB'

#Razer Blade
output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO_images'  # Output image directory
# Set up paths and filenames
silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO'  # Directory containing silo files
filebase = 'BN_grad_d2l4n128'  # Base name of the silo files
#database_path = '/home/tony/Desktop/NebulaPy/NebulaPy-DB'

ion_list = ['H1+', 'He1+', 'C2+', 'N1+', 'N2+', 'O1+', 'O2+', 'Ne1+', 'Ne2+', 'S1+', 'S2+', 'S3+']

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

# Initialize the Pion class from NebulaPy to handle simulation data and processing
pion = nebula.pion(batched_silos, verbose=True)
# Load chemistry data from the first time instant's silo file
pion.load_chemistry()
# Load the geometry of the simulation grid from the first silo file
pion.load_geometry(batched_silos[0])
geometry = pion.geometry_container
N_grid_level = geometry['Nlevels']  # Number of grid levels

print(" ---------------------------")
print(" task: finding dominant spectral lines for ions")
print(f" ion check:")
# check if the ion in the current silo file and create a new list with ions that pass the check
ion_list = [ion for ion in ion_list if pion.ion_check(ion, top_ion_check=True, terminate=False)]
print(f" ion list: {ion_list}")

filename = filebase + '_dominantLines.txt'
outfile = os.path.join(output_dir, filename)

with open(outfile, "w") as file:
    file.write(f"Generated by {util.nebula_version()}\n")
    file.write("\n")
    file.write("Task: Finding dominant spectral lines for ions\n")
    file.write(f"Ion List: {ion_list}\n")
    file.write("\n\n")

runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time

    print(f" ---------------------------")
    # Print the current simulation time instant
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")
    with open(outfile, "a") as file:
        file.write(f"Step: {step} | Simulation Time: {sim_time:.6e}\n\n")

    # Extract temperature for the current time instant
    temperature = pion.get_parameter('Temperature', silo_instant)  # Retrieve temperature
    # Retrieve the electron number density
    ne = pion.get_ne(silo_instant)

    # ion ************************************************************************
    for ion in ion_list:

        print(f" computing the dominant lines of the ion {ion} for each grid level(s)")
        ion_lines = nebula.line_emission(ion=ion, verbose=True)
        dominant_lines = ion_lines.get_dominant_lines(temperature=temperature, ne=ne, Nlines=10,
                                                      geometry_container=geometry)

        # Convert data into pandas DataFrame for better alignment
        df_wvls = pd.DataFrame(dominant_lines['wvls'])
        df_emiss = pd.DataFrame(dominant_lines['emiss'])

        with open(outfile, "a") as file:
            # Write ion information
            file.write(f"Ion: {ion}\n\n")
            # Write header
            file.write("Grid Level | " + " | ".join([f"Level {i}" for i in range(len(df_wvls))]) + "\n")
            file.write("-" * (12 + (len(df_wvls.columns) * 10)) + "\n")

            # Write wavelengths
            file.write("Wavelength  | ")
            for i in range(len(df_wvls)):
                file.write("  ".join(f"{wvl:.6e}" for wvl in df_wvls.iloc[i].values) + " | ")
            file.write("\n")

            # Write emissivity
            file.write("Emissivity  | ")
            for i in range(len(df_emiss)):
                file.write("  ".join(f"{emis:.6e}" for emis in df_emiss.iloc[i].values) + " | ")
            file.write("\n\n")
    # end of each ion ********************************************************************

    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")
    with open(outfile, "a") as file:
        file.write(f"Runtime: {runtime:.4e} s | dt: {dt:.4e} s\n\n")