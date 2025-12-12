"""
Number Density Map for Fe
2D PION simulation:
Author: Arun Mathew
Date: 11 Dec 2025
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
import warnings
# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")

# Mimir -> Set up paths and filenames
OutputDir = '/net/maedoc.ap.dias.ie/maedoc/home_cr/arun/Desktop/plots/CWBs-2026'  # Output image directory
SiloDir = '/mnt/massive-stars/data/wr140-nemo/wr140_NEMO_d07e13_d2l6n128'  # Directory containing silo files
Filebase = 'wr140_NEMO_d07e13_d2l6n128'  # Base name of the silo files
start_time = None
finish_time = None
time_unit = None
out_frequency = None
ion_list = ['Fe25+', 'Fe26+']

#Razer Blade
#OutputDir = '/home/tony/Desktop/multi-ion-bowshock/sim-output/coolmap'  # Output image directory
# Set up paths and filenames
#SiloDir = '/home/tony/Desktop/multi-ion-bowshock/high-res-silos-200kyr'  # Directory containing silo files
#Filebase = 'Ostar_mhd-nemo-dep_d2n0384l3'  # Base name of the silo files
#start_time = None
#finish_time = None
#time_unit = None
#out_frequency = None
#ion_list = ['Fe25+', 'Fe26+']

# Batch the silo files according to the time instant
batched_silos = util.batch_silos(
    SiloDir,
    Filebase,
    start_time=start_time,
    finish_time=finish_time,
    time_unit=time_unit,
    out_frequency=out_frequency
)

# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# Calculates and stores geometric grid parameters.
# For example, in a spherical geometry, it extracts radius and shell volumes
# from the first silo file in the batch and saves them into a geometry container.
pion.load_geometry(scale='cm')
N_grid_level = pion.geometry_container['Nlevel']
mesh_edges_min = pion.geometry_container['edges_min']
mesh_edges_max = pion.geometry_container['edges_max']

print(mesh_edges_min)
print(mesh_edges_max)

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
pion.load_chemistry()

runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time

    print(f" ---------------------------")
    # Print the current simulation time instant
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # Create a grid of subplots for visualizing ionisation fractions
    fig, axes = plt.subplots(1, 2)
    axes = axes.flat

    # Plot ionisation fractions for each ion
    for i, (ax, ion) in enumerate(zip(axes, ion_list)):

        n_ion = pion.get_ion_number_density(ion, silo_instant)

        ax.set_xlim(mesh_edges_min[0][0].value, mesh_edges_max[0][0].value)
        ax.set_ylim(mesh_edges_min[0][1].value, mesh_edges_max[0][1].value)

        ax.set_xlim(-7.0e13, 7.0e13)
        ax.set_ylim(0, 7.0e13)

        '''
        if not i in (8, 9):
            ax.axes.get_xaxis().set_visible(False)  # Remove the x-axis.
        if not i in (0, 2, 4, 6, 8):
            ax.axes.get_yaxis().set_visible(False)  # Remove the x-axis.
        if i in (8, 9):
            ax.set_xlabel('z (' + distance_unit_str + ')', fontsize=15)
        if i in (0, 2, 4, 6, 8):
            ax.set_ylabel('R (' + distance_unit_str + ')', fontsize=15)
        '''

        # Plot data for each grid level
        for level in range(N_grid_level):
            plot_data = n_ion[level]
            extents = [mesh_edges_min[level][0].value, mesh_edges_max[level][0].value,
                       mesh_edges_min[level][1].value, mesh_edges_max[level][1].value]
            image = ax.imshow(plot_data, interpolation='nearest', cmap='Reds',
                              extent=extents, origin='lower',
                              #vmin=0, vmax=1
                              )

        ax.text(0.1, 0.8, ion_list[i], transform=ax.transAxes, fontsize=14,
                #color=ion_text_color[i]
                )
        ax.tick_params(axis='both', which='major', labelsize=12)

    # Add a single colorbar for the entire figure
    cbar_ax = fig.add_axes([0.125, 0.91, 0.775, 0.015])  # [left, bottom, width, height]
    fig.colorbar(image, cax=cbar_ax, orientation='horizontal', ticks=MultipleLocator(0.2))
    plt.subplots_adjust(hspace=0.000, wspace=0.00)

    # Save the figure
    Filename = f"{Filebase}_FeNumDen_{sim_time.value:.2f}kyr.png"
    plt.savefig(os.path.join(OutputDir, Filename), bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f" time: {sim_time:.6e}, Saved snapshot {step} to {Filename}")


    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")



