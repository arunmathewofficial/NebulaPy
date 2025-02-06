"""
Line Luminosity

Description:
This script calculates the temporal evolution of luminosity for selected spectral lines
in an astrophysical simulation. The lines correspond to different ions,
including He1+, C2+, N1+, N2+, O1+, O2+, Ne1+, Ne2+, S1+, S2+, and S3+, and their
corresponding emission lines. The luminosity for each line is computed using the
simulation data and stored for later analysis.

Features:
- Reads silo files generated from the simulation to extract temperature, electron density, and ion number densities.
- Computes the line luminosity for specific emission lines over a time range.
- Outputs the results to a text file, with each row representing a time step and the luminosities of different lines.

Author: Arun Mathew
Date: 01 Feb 2025
"""

import os
import time
import warnings
from NebulaPy.tools import util
import NebulaPy.src as nebula

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")

# Input-output file configuration for a low-resolution simulation on MIMIR
output_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/low-res/time-lines-luminosity'  # Change as needed
silo_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/low-res/silo'
filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'  # Base name of the silo files
filename = filebase + '_lines_luminosity_LowRes.txt'
outfile = os.path.join(output_dir, filename)

# Input-output file configuration for a medium-resolution simulation on MIMIR
output_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/med-res/time-lines-luminosity'  # Change as needed
silo_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/med-res/silo'
filebase = 'Ostar_mhd-nemo-dep_d2n0192l3'  # Base name of the silo files
filename = filebase + '_lines_luminosity_MedRes.txt'
outfile = os.path.join(output_dir, filename)

# Input-output file configuration for a low-resolution simulation on Razer Blade machine.
#output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/out'  # Change as needed
#silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/silo'
#filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'  # Base name of the silo files
# Prepare output file for results
#filename = filebase + '_lines_luminosity_LowRes.txt'
#outfile = os.path.join(output_dir, filename)

# Batch the silo files for analysis within the specified time range
batched_silos = util.batch_silos(
    silo_dir,
    filebase,
    start_time=None,
    finish_time=None,
    time_unit='kyr',
    out_frequency=2
)

# Initialize the PION class to handle simulation data
pion = nebula.pion(batched_silos, verbose=True)
# Load chemistry and geometry data
pion.load_chemistry()
pion.load_geometry(scale='cm')

print(f" ---------------------------")
print(" task: calculating temporal evolution of luminosity for given lines")

# Set up the ion and line emission parameters
# Define the ions and their respective emission lines
He1P_pion_ion = 'He1+'
He1P_lines = [303.78, 303.786, 256.317, 256.318, 243.026, 243.027]  # Emission line(s) of interest
print(rf" {He1P_pion_ion:<4} lines: {', '.join(map(str, He1P_lines))}  Angstrom")

C2P_pion_ion = 'C2+'
C2P_lines = [1906.683, 1908.734, 977.02]  # Emission line(s) of interest
print(rf" {C2P_pion_ion:<4} lines: {', '.join(map(str, C2P_lines))}  Angstrom")

N1P_pion_ion = 'N1+'
N1P_lines = [6585.273, 6549.861, 1218026.8, 2053388.09]
print(rf" {N1P_pion_ion:<4} lines: {', '.join(map(str, N1P_lines))}  Angstrom")

N2P_pion_ion = 'N2+'
N2P_lines = [573394.5, 989.799, 1752.16, 1749.674, 1753.995, 1748.646]
print(rf" {N2P_pion_ion:<4} lines: {', '.join(map(str, N2P_lines))}  Angstrom")

O1P_pion_ion = 'O1+'  # The ion of interest (Oxygen II)
O1P_lines = [3729.844, 3727.092, 7331.722, 2470.97, 2471.094, 7321.094, 7322.177, 7332.808]
print(rf" {O1P_pion_ion:<4} lines: {', '.join(map(str, O1P_lines))} Angstrom")

O2P_pion_ion = 'O2+'  # The ion of interest (Oxygen III)
O2P_lines = [5008.24, 883564.0, 518145.0, 4960.295, 1666.15, 4364.436, 832.929, 833.715, 1660.809]
print(rf" {O2P_pion_ion:<4} lines: {', '.join(map(str, O2P_lines))} Angstrom")

Ne1P_pion_ion = 'Ne1+'
Ne1P_lines = [128139.42]
print(rf" {Ne1P_pion_ion:<4} lines: {', '.join(map(str, Ne1P_lines))} Angstrom")

Ne2P_pion_ion = 'Ne2+'
Ne2P_lines = [155545.19, 3869.849, 3968.585, 360230.55, 489.495, 491.041]
print(rf" {Ne2P_pion_ion:<4} lines: {', '.join(map(str, Ne2P_lines))} Angstrom")

S1P_pion_ion = 'S1+'
S1P_lines = [6718.295, 6732.674, 4069.749, 10323.317, 4077.5]
print(rf" {S1P_pion_ion:<4} lines: {', '.join(map(str, S1P_lines))} Angstrom")

S2P_pion_ion = 'S2+'
S2P_lines = [335008.38, 9532.252, 187055.74, 9070.048, 6313.649, 3722.454]
print(rf" {S2P_pion_ion:<4} lines: {', '.join(map(str, S2P_lines))} Angstrom")

S3P_pion_ion = 'S3+'
S3P_lines = [105104.95]
print(rf" {S3P_pion_ion:<4} lines: {', '.join(map(str, S3P_lines))} Angstrom")

# Initialize the emission line calculations for each ion
He1P_line_emission = nebula.line_emission(He1P_pion_ion, verbose=True)
C2P_line_emission = nebula.line_emission(C2P_pion_ion, verbose=True)
N1P_line_emission = nebula.line_emission(N1P_pion_ion, verbose=True)
N2P_line_emission = nebula.line_emission(N2P_pion_ion, verbose=True)
O1P_line_emission = nebula.line_emission(O1P_pion_ion, verbose=True)
O2P_line_emission = nebula.line_emission(O2P_pion_ion, verbose=True)
Ne1P_line_emission = nebula.line_emission(Ne1P_pion_ion, verbose=True)
Ne2P_line_emission = nebula.line_emission(Ne2P_pion_ion, verbose=True)
S1P_line_emission = nebula.line_emission(S1P_pion_ion, verbose=True)
S2P_line_emission = nebula.line_emission(S2P_pion_ion, verbose=True)
S3P_line_emission = nebula.line_emission(S3P_pion_ion, verbose=True)

# Check the requested lines in the database for each ion
print(f" ---------------------------")
print(f" checking requested lines in database:")
He1P_line_emission.line_batch_check(He1P_lines)
C2P_line_emission.line_batch_check(C2P_lines)
N1P_line_emission.line_batch_check(N1P_lines)
N2P_line_emission.line_batch_check(N2P_lines)
O1P_line_emission.line_batch_check(O1P_lines)
O2P_line_emission.line_batch_check(O2P_lines)
Ne1P_line_emission.line_batch_check(Ne1P_lines)
Ne2P_line_emission.line_batch_check(Ne2P_lines)
S1P_line_emission.line_batch_check(S1P_lines)
S2P_line_emission.line_batch_check(S2P_lines)
S3P_line_emission.line_batch_check(S3P_lines)

# Get geometry information
geometry = pion.geometry_container
N_grid_level = geometry['Nlevel']
grid_mask = geometry['mask']
cell_volume = pion.get_cylindrical_cell_volume().value

# Write initial header to the output file
with open(outfile, "w") as file:
    file.write(f"#File generated by {util.nebula_version()}\n\n")
    file.write("#Task: Temporal evolution of luminosity for selected spectral lines\n\n")
    file.write(f"#PION Simulation Reference: NEMO Bowshock 2025 {filebase}\n\n")
    file.write("#Dataset Description:\n")
    file.write("#This dataset provides the luminosity evolution of selected spectral lines over time.\n")
    file.write("#Each row represents a different time step, with:\n")
    file.write("# - The first column indicating time (in kyr).\n")
    file.write("# - Subsequent columns representing the luminosity of specific spectral lines (in erg/s).\n\n")

# Loop over each time instant in the batched silo files
runtime = 0.0
write_heading = True
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()

    print(f" ---------------------------")
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # Extract temperature and electron number density
    temperature = pion.get_parameter('Temperature', silo_instant)
    ne = pion.get_ne(silo_instant)

    He1P_num_density = pion.get_ion_number_density(He1P_pion_ion, silo_instant)  # Retrieve species number density
    C2P_num_density = pion.get_ion_number_density(C2P_pion_ion, silo_instant)  # Retrieve species number density
    N1P_num_density = pion.get_ion_number_density(N1P_pion_ion, silo_instant)  # Retrieve species number density
    N2P_num_density = pion.get_ion_number_density(N2P_pion_ion, silo_instant)  # Retrieve species number density
    O1P_num_density = pion.get_ion_number_density(O1P_pion_ion, silo_instant)  # Retrieve species number density
    O2P_num_density = pion.get_ion_number_density(O2P_pion_ion, silo_instant)  # Retrieve species number density
    Ne1P_num_density = pion.get_ion_number_density(Ne1P_pion_ion, silo_instant)  # Retrieve species number density
    Ne2P_num_density = pion.get_ion_number_density(Ne2P_pion_ion, silo_instant)  # Retrieve species number density
    S1P_num_density = pion.get_ion_number_density(S1P_pion_ion, silo_instant)  # Retrieve species number density
    S2P_num_density = pion.get_ion_number_density(S2P_pion_ion, silo_instant)  # Retrieve species number density
    S3P_num_density = pion.get_ion_number_density(S3P_pion_ion, silo_instant)  # Retrieve species number density


    # 1. Calculate the line luminosity for the specific emission line
    He1P_lines_luminosity = He1P_line_emission.line_luminosity_cylindrical(
        lines=He1P_lines,
        temperature=temperature,
        ne=ne,
        species_density=He1P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 2. Calculate the line luminosity for the specific emission line
    C2P_lines_luminosity = C2P_line_emission.line_luminosity_cylindrical(
        lines=C2P_lines,
        temperature=temperature,
        ne=ne,
        species_density=C2P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 3. Calculate the line luminosity for the specific emission line
    N1P_lines_luminosity = N1P_line_emission.line_luminosity_cylindrical(
        lines=N1P_lines,
        temperature=temperature,
        ne=ne,
        species_density=N1P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 4. Calculate the line luminosity for the specific emission line
    N2P_lines_luminosity = N2P_line_emission.line_luminosity_cylindrical(
        lines=N2P_lines,
        temperature=temperature,
        ne=ne,
        species_density=N2P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 5. Calculate the line luminosity for the specific emission line
    O1P_lines_luminosity = O1P_line_emission.line_luminosity_cylindrical(
        lines=O1P_lines,
        temperature=temperature,
        ne=ne,
        species_density=O1P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 6. Calculate the line luminosity for the specific emission line
    O2P_lines_luminosity = O2P_line_emission.line_luminosity_cylindrical(
        lines=O2P_lines,
        temperature=temperature,
        ne=ne,
        species_density=O2P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 7. Calculate the line luminosity for the specific emission line
    Ne1P_lines_luminosity = Ne1P_line_emission.line_luminosity_cylindrical(
        lines=Ne1P_lines,
        temperature=temperature,
        ne=ne,
        species_density=Ne1P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 8. Calculate the line luminosity for the specific emission line
    Ne2P_lines_luminosity = Ne2P_line_emission.line_luminosity_cylindrical(
        lines=Ne2P_lines,
        temperature=temperature,
        ne=ne,
        species_density=Ne2P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 9. Calculate the line luminosity for the specific emission line
    S1P_lines_luminosity = S1P_line_emission.line_luminosity_cylindrical(
        lines=S1P_lines,
        temperature=temperature,
        ne=ne,
        species_density=S1P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 10. Calculate the line luminosity for the specific emission line
    S2P_lines_luminosity = S2P_line_emission.line_luminosity_cylindrical(
        lines=S2P_lines,
        temperature=temperature,
        ne=ne,
        species_density=S2P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    # 11. Calculate the line luminosity for the specific emission line
    S3P_lines_luminosity = S3P_line_emission.line_luminosity_cylindrical(
        lines=S3P_lines,
        temperature=temperature,
        ne=ne,
        species_density=S3P_num_density,
        cell_volume=cell_volume,
        grid_mask=grid_mask)

    #Combine all dictionaries
    Dictionaries = He1P_lines_luminosity
    Dictionaries.update(C2P_lines_luminosity)
    Dictionaries.update(N1P_lines_luminosity)
    Dictionaries.update(N2P_lines_luminosity)
    Dictionaries.update(O1P_lines_luminosity)
    Dictionaries.update(O2P_lines_luminosity)
    Dictionaries.update(Ne1P_lines_luminosity)
    Dictionaries.update(Ne2P_lines_luminosity)
    Dictionaries.update(S1P_lines_luminosity)
    Dictionaries.update(S2P_lines_luminosity)
    Dictionaries.update(S3P_lines_luminosity)

    if write_heading:
        with open(outfile, "a") as file:
            # Write the time (or any other desired variable)
            file.write(f"Time ")
            # Write the keys from Dictionaries
            file.write(" ".join(f"{key}" for key in Dictionaries.keys()))

            file.write("\n")
        write_heading = False

    with open(outfile, "a") as file:
        file.write(f"{sim_time.value:.6e} ")
        file.write(" ".join(f"{v:.6e}" for v in Dictionaries.values()))
        file.write("\n")

    del Dictionaries

    # Update runtime
    silo_instant_finish_time = time.time()
    dt = silo_instant_finish_time - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | step runtime: {dt:.4e} s")
