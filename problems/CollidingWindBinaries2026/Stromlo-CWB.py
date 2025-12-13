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
start_time = 1.0e6  # in sec
finish_time = None
time_unit = 'sec'
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
    silo_instant_start_time = time.time()

    print(f" ---------------------------")
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # 2 rows (upper/lower) x 2 cols (Fe25+/Fe26+)
    fig, axes = plt.subplots(2, 2, figsize=(4, 7), sharex=True, sharey='row')
    # axes[row, col]: row 0 = upper, row 1 = lower
    last_image = None

    for j, ion in enumerate(ion_list):  # j=0 Fe25+, j=1 Fe26+
        n_ion = pion.get_ion_number_density(ion, silo_instant)

        axU = axes[0, j]
        axL = axes[1, j]

        # Common x-limits
        for ax in (axU, axL):
            ax.set_xlim(-7.0e13, 2.0e13)
            ax.tick_params(axis='both', which='major', labelsize=12)

        # Upper/lower y-limits
        axU.set_ylim(0.0,  2.0e14)
        axL.set_ylim(-2.0e14, 0.0)

        # Plot data for each AMR/grid level
        for level in range(N_grid_level):
            plot_data = np.log10(n_ion[level])

            x_min = mesh_edges_min[level][0].value
            x_max = mesh_edges_max[level][0].value
            y_min = mesh_edges_min[level][1].value
            y_max = mesh_edges_max[level][1].value

            # --- Upper hemisphere (as-is) ---
            extU = [x_min, x_max, y_min, y_max]
            last_image = axU.imshow(
                plot_data,
                interpolation='nearest',
                cmap='viridis',
                extent=extU,
                origin='lower',
                vmin=0, vmax=4.0
            )

            # --- Lower hemisphere (mirror of upper) ---
            # Flip vertically so it plots correctly with origin='lower'
            plot_data_mirror = np.flipud(plot_data)
            # Map y-range to negative values
            extL = [x_min, x_max, -y_max, -y_min]
            last_image = axL.imshow(
                plot_data_mirror,
                interpolation='nearest',
                cmap='viridis',
                extent=extL,
                origin='lower',
                vmin=0, vmax=4.0
            )

        # Panel labels
        axU.text(0.05, 0.90, f"n({ion})", transform=axU.transAxes, fontsize=13)

    # One shared colorbar (use the last imshow handle)
    cbar_ax = fig.add_axes([0.125, 0.93, 0.775, 0.02])
    fig.colorbar(last_image, cax=cbar_ax, orientation='horizontal')

    plt.subplots_adjust(hspace=0.0, wspace=0.00, top=0.90)

    Filename = f"{Filebase}_FeNumDen_{sim_time.value:.2f}kyr.png"
    plt.savefig(os.path.join(OutputDir, Filename), bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f" time: {sim_time:.6e}, Saved snapshot {step} to {Filename}")

    dt = time.time() - silo_instant_start_time
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")




