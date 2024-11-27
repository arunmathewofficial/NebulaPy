"""
Ion Line Analyzer
Description:
    This script identifies the dominant spectral lines for a list of ions across simulation snapshots stored in silo files.
    It processes the simulation data to determine the most significant emission lines for each ion, based on the temperature
    and electron density at different grid levels. The output is saved to a text file with the computed results for analysis.

Features:
    - Loads and processes chemistry and grid geometry data from simulation snapshots.
    - Filters a user-defined list of ions to identify those present in the simulation.
    - Computes the top dominant emission lines for each ion at different grid levels.
    - Outputs detailed results, including wavelengths and emissivities, in a structured format.

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
# Restrict OpenMP threads
os.environ["OMP_NUM_THREADS"] = "4"  # Limit OpenMP threads for better performance control

# Output directory and filebase configuration for MIMIR
output_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/low-res'
silo_dir = '/mnt/massive-stars/data/nemo/simple-bowshock'
filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'

# OutPut directory and filebase configurationRazer Blade
#output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO_images'
#silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO'
#filebase = 'BN_grad_d2l4n128'  # Base name of the silo files

# List of ions to analyze
ion_list = ['H1+', 'He1+', 'C2+', 'N1+', 'N2+', 'O1+', 'O2+', 'Ne1+', 'Ne2+', 'S1+', 'S2+', 'S3+']

# Batch the silo files for analysis within the specified time range
start_time = 100.0
finish_time = None
batched_silos = util.batch_silos(
    silo_dir,
    filebase,
    start_time=start_time,
    finish_time=None,
    time_unit='kyr',
    out_frequency=20
)

# Initialize the Pion class to handle simulation data and load chemistry and geometry
pion = nebula.pion(batched_silos, verbose=True)
pion.load_chemistry()
pion.load_geometry(batched_silos[0])
geometry = pion.geometry_container
N_grid_level = geometry['Nlevels']

print(f" ---------------------------")
# Filter the ion list to include only those present in the simulation
print(" task: finding dominant spectral lines for ions")
print(" ion check: ")
ion_list = [ion for ion in ion_list if pion.ion_check(ion, top_ion_check=True, terminate=False)]
print(f" filtered ion list: {ion_list}")

# Prepare output file for results
filename = filebase + '_dominantLines.txt'
outfile = os.path.join(output_dir, filename)

# Write the initial header to the output file
with open(outfile, "w") as file:
    file.write(f"Generated by {util.nebula_version()}\n\n")
    file.write("Task: Finding dominant spectral lines for ions\n")
    file.write(f"Filtered ion list: {ion_list}\n\n")

runtime = 0.0

# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()

    print(f" ---------------------------")
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")
    with open(outfile, "a") as file:
        file.write(f"Step: {step} | Simulation Time: {sim_time:.6e}\n\n")

    # Extract temperature and electron number density for the current time instant
    temperature = pion.get_parameter('Temperature', silo_instant)
    ne = pion.get_ne(silo_instant)

    # Analyze each ion
    for ion in ion_list:
        print(f" computing dominant lines for ion: {ion}")
        ion_lines = nebula.line_emission(ion=ion, verbose=True)
        dominant_lines = ion_lines.get_dominant_lines(
            temperature=temperature, ne=ne, Nlines=10, geometry_container=geometry
        )

        # Format the results as a DataFrame
        df_wvls = pd.DataFrame(dominant_lines['wvls'])
        df_emiss = pd.DataFrame(dominant_lines['emiss'])

        # Write results for the current ion to the output file
        with open(outfile, "a") as file:
            file.write(f"Ion: {ion}\n\n")
            file.write("Grid Level | " + " | ".join([f"Level {i}" for i in range(len(df_wvls))]) + "\n")
            file.write("-" * (12 + (len(df_wvls.columns) * 4)) + "\n")
            file.write("Lines (Amstrong) | " + " | ".join(df_wvls.applymap("{:.6f}".format).values.flatten()) + "\n")
            file.write("Emissivity (ergs s^-1 str^-1) | " + " | ".join(df_emiss.applymap("{:.6e}".format).values.flatten()) + "\n\n")

    # Update and log the runtime for the current step
    silo_instant_finish_time = time.time()
    dt = silo_instant_finish_time - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | step runtime: {dt:.4e} s")
    with open(outfile, "a") as file:
        file.write(f"Runtime: {runtime:.4e} s | Step Runtime: {dt:.4e} s\n\n")
