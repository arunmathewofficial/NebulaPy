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

#import numpy as np
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

# OutPut directory and filebase configuration for Razer Blade
output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/out'
silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/silo'
filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'  # Base name of the silo files

# Batch the silo files for analysis within the specified time range
batched_silos = util.batch_silos(
    silo_dir,
    filebase,
    start_time=150,
    finish_time=None,
    time_unit= 'kyr',
    out_frequency=5
)

# Initialize the PION class to handle simulation data, load geometry, and chemistry
pion = nebula.pion(batched_silos, verbose=True)
# following are routine steps when processing PION silo files.
# generate chemistry container
pion.load_chemistry()
# generate geometry container
pion.load_geometry(scale='cm')

print(f" ---------------------------")
# The task begins here: identifying the dominant spectral lines for the list of ions.
# List of ions to analyze
ion_list = ['H1+', 'He1+', 'C2+', 'N1+', 'N2+', 'O1+', 'O2+', 'Ne1+', 'Ne2+', 'S1+', 'S2+', 'S3+']
# Filter the ion list to include only those present in the simulation
print(" task: determining the dominant spectral lines for the given ions")
# Check the batch of ions listed in chemistry container
ion_list = pion.ion_batch_check(ion_list=ion_list, top_ion_check=True, terminate=False)

# Prepare output file for results
filename = filebase + '_dominantLines.txt'
outfile = os.path.join(output_dir, filename)

# get geometry container
geometry = pion.geometry_container
N_grid_level = geometry['Nlevel']
grid_mask = geometry['mask']
cell_volume = pion.get_cylindrical_cell_volume().value


# Write the initial header to the output file
with open(outfile, "w") as file:
    file.write(f"Generated by {util.nebula_version()}\n\n")
    file.write("Task: Determining the dominant spectral lines for the following ions\n")
    # Write column headers
    file.write("   ".join(map(str, range(len(ion_list)))) + "\n")
    # Write ion species
    file.write("   ".join(ion_list) + "\n\n")
    file.write(f"FileBase: {filebase}\n\n")
    file.write("Spectral Line are given in Armstrong\n")
    file.write("Line Luminosity in erg s^-1\n\n")

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
        species_line_emission = nebula.line_emission(ion=ion, verbose=True)
        # Retrieve the ion number density
        n_species = pion.get_ion_number_density(ion, silo_instant)

        dominant_lines = species_line_emission.get_species_dominant_lines(
            temperature=temperature, ne=ne, species_density=n_species,
            cell_volume=cell_volume, grid_mask=grid_mask, Nlines=16)

        if 'lines' not in dominant_lines:
            # Write results for the current ion to the output file
            with open(outfile, "a") as file:
                file.write(f"No free-bound emission associated with {dominant_lines['spectroscopic']}\n\n")
                print(" skipping ...")

        else:
            # Convert the single string to a list so it can be combined with the array of numbers
            spectral_lines = [dominant_lines['spectroscopic'] + " " + str(line) for line in dominant_lines['lines']]
            # Create a DataFrame for spectroscopic lines and their corresponding luminosities
            dataframe_lines = pd.DataFrame({'Spectral Line': spectral_lines})
            dataframe_luminosity = pd.DataFrame({'Luminosity (erg/s)': dominant_lines['luminosity']})

            # Combine the two DataFrames
            final_dataframe = pd.concat([dataframe_lines, dataframe_luminosity], axis=1)

            # Write results for the current ion to the output file
            with open(outfile, "a") as file:
                final_dataframe.to_csv(file, index=False, header=False, sep='\t', float_format="%.5e")
                file.write("\n")

    # Update and log the runtime for the current step
    silo_instant_finish_time = time.time()
    dt = silo_instant_finish_time - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | step runtime: {dt:.4e} s")
    with open(outfile, "a") as file:
        file.write(f"Runtime: {runtime:.4e} s | Step Runtime: {dt:.4e} s\n\n")
