"""
2D Emissivity Map Generator for Multi-Line Single Ion Emissions
Description: Generate Emissivity Map from 2D PION simulation silo files
Author: Arun Mathew
Date: 21 Nov 2024
"""

import NebulaPy.src as nebula  # Import the NebulaPy package for astrophysical simulation processing
from NebulaPy.tools import util  # Utility functions from NebulaPy
from pypion.ReadData import ReadData  # Import ReadData from pypion for reading simulation data
import astropy.units as unit  # Astropy units for consistent handling of physical units
import matplotlib.pyplot as plt  # Matplotlib for creating plots
from mpl_toolkits.axes_grid1 import make_axes_locatable  # For adding colorbars to plots
from matplotlib.ticker import MultipleLocator, ScalarFormatter  # For formatting plot ticks and colorbar
import numpy as np  # Numpy for numerical operations
import time  # For tracking runtime
import warnings  # For managing warnings in the code
import os  # For file system operations

# Suppress specific warnings (divide by zero encountered in log10)
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")

# Define the output directory for the emissivity maps
output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO_images'  # Output image directory

# Set the directory containing the silo files and the base name of these files
silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO'  # Directory with simulation data (silo files)
filebase = 'BN_grad_d2l4n128'  # Base name of the silo files

# Define the ion of interest and the emission lines to calculate
pion_ion = 'O2+'  # Ion of interest (Oxygen IV)
lines = [4960.295, 5008.240]  # Wavelengths of the emission lines in Angstroms

# Print the line luminosity calculation information
print(rf" calculating line luminosity of {pion_ion} lines: {lines} Angstrom")

# Batch the silo files based on time intervals (start_time, finish_time, etc.)
batched_silos = util.batch_silos(
    silo_dir,  # Directory of silo files
    filebase,  # Base file name
    start_time=None,  # Optional start time
    finish_time=None,  # Optional finish time
    time_unit=None,  # Optional time units
    out_frequency=None  # Optional output frequency
)

# Initialize the Pion class from NebulaPy to handle simulation data and processing
pion = nebula.pion(batched_silos, verbose=True)
# Load chemistry data from the first time instant's silo file
pion.load_chemistry()
# Load the geometry of the simulation grid from the first silo file
pion.load_geometry(batched_silos[0])

# Extract geometry-related parameters
coordinate_sys = pion.geometry_container['coordinate_sys']  # Coordinate system
N_grid_level = pion.geometry_container['Nlevels']  # Number of grid levels
mesh_edges_min = pion.geometry_container['edges_min']  # Minimum mesh edges
mesh_edges_max = pion.geometry_container['edges_max']  # Maximum mesh edges

# Convert line wavelengths to string format for later use
lines_str = [str(line) for line in lines]

# Define a tolerance for electron density (avoid division by zero)
electron_tolerance = 1.E-08
# Define the output directory for each ion's emissivity maps
ion_name = pion_ion.replace('+', 'p')
ion_output_dir = os.path.join(output_dir, ion_name)
os.makedirs(ion_output_dir, exist_ok=True)  # Create directory if it doesn't exist

# Initialize the runtime counter
runtime = 0.0

# Loop over each time instant (step) in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Track start time for each simulation step

    # Print the current simulation time instant
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" ---------------------------")
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # Extract physical parameters for the current simulation step
    temperature = pion.get_parameter('Temperature', silo_instant)  # Temperature data
    ne = pion.get_ne(silo_instant)  # Electron density
    rows = len(temperature[0])  # Number of inner arrays in 2D temperature data

    # Initialize a list to store emissivity maps
    emissivity_map = [np.zeros(shape) for shape in [arr.shape for arr in temperature]]
    emissivity_map_dict = {line: emissivity_map for line in lines_str}  # Map emission lines to maps

    # Loop over each grid level in the simulation
    for level in range(N_grid_level):
        ne[level] = np.array(ne[level])
        ne[level][ne[level] == 0] = electron_tolerance  # Replace zero values with tolerance

        # Loop over each row in the grid
        for row in range(rows):
            print(f" computing emissivity for line(s) at level {level}, row {row}", end='\r')

            # Get temperature and electron density for the current row
            temperature_row = temperature[level][row]
            ne_row = ne[level][row]

            # Create a Chianti object for line emissivity calculation
            chinati = nebula.chianti(temperature_row, ne_row, pion_ion=pion_ion, pion_elements=None, verbose=False)

            # Calculate the line emissivity for each specified line
            lines_emissivity_row = chinati.get_line_emissivity(line_list=lines)

            # Store the computed emissivity for the current line
            for line in lines_str:
                emissivity_map_dict[line][level][row] = lines_emissivity_row[line]

        print(f" completed emissivity calculations for line(s) at level {level}")

    # Loop through each emission line and plot the corresponding emissivity map
    for line in lines_str:
        fig, ax = plt.subplots(figsize=(8, 6))  # Create a new figure for the plot

        # Add text annotation for simulation time
        ax.text(0.05, 0.9, 'time = %5.2f kyr' % sim_time.value, transform=ax.transAxes,
                fontsize=12, color='white')

        # Set plot limits based on mesh edges
        ax.set_xlim(mesh_edges_min[0][0].value, mesh_edges_max[0][0].value)
        ax.set_ylim(mesh_edges_min[0][1].value, mesh_edges_max[0][1].value)

        # Plot the emissivity map using a logarithmic scale
        for level in range(N_grid_level):
            plot_data = np.log10(emissivity_map_dict[line][level])  # Logarithmic emissivity data
            extents = [mesh_edges_min[level][0].value, mesh_edges_max[level][0].value,
                       mesh_edges_min[level][1].value, mesh_edges_max[level][1].value]

            # Create the image plot
            image = ax.imshow(plot_data, interpolation='nearest', cmap='inferno', extent=extents, origin='lower',
                              vmin=-25, vmax=-21)

        # Add colorbar with appropriate formatting
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        colorbar = plt.colorbar(image, cax=cax, ticks=MultipleLocator(1))
        colorbar.ax.yaxis.set_major_formatter(ScalarFormatter())
        colorbar.ax.yaxis.get_major_formatter().set_scientific(False)
        colorbar.ax.yaxis.get_major_formatter().set_useOffset(False)
        colorbar.update_ticks()

        # Set axis labels and text annotations
        ax.set_xlabel('z (pc)', fontsize=12)
        ax.set_ylabel('R (pc)', fontsize=12)
        ax.text(0.65, 0.9, pion_ion, transform=ax.transAxes, fontsize=12, color='white')

        # Customize tick labels
        ax.tick_params(axis='both', which='major', labelsize=13)

        # Save the figure to a file
        filename = f"{filebase}_coolmap_{ion_name}[{line}]_{str(step).zfill(4)}.png"
        filepath = os.path.join(ion_output_dir, filename)
        plt.savefig(filepath, bbox_inches="tight", dpi=300)
        plt.close(fig)

        print(f" saving emissivity map image for ion {pion_ion}, line {line}, to {filename}")

    # Track runtime for each time step
    silo_instant_finish_time = time.time()
    dt = silo_instant_finish_time - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")

