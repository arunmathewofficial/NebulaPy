import numpy as np
import os
import time
from NebulaPy.tools import util
import NebulaPy.src as nebula
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, ScalarFormatter
from matplotlib.gridspec import GridSpec

output_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO_images'  # Output image directory
# Set up paths and filenames
silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sims/HHeCNO'  # Directory containing silo files
filebase = 'BN_grad_d2l4n128'  # Base name of the silo files
ion = 'H'

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

# Initialize the Pion class from NebulaPy, which handles the simulation data
pion = nebula.pion(batched_silos, verbose=True)

# Calculates and stores geometric grid parameters.
# For example, in a spherical geometry, it extracts radius and shell volumes
# from the first silo file in the batch and saves them into a geometry container.
pion.load_geometry(batched_silos[0], scale='pc')
print(pion.geometry_container)
N_grid_level = pion.geometry_container['Nlevels']
mesh_edges_min = pion.geometry_container['edges_min']
mesh_edges_max = pion.geometry_container['edges_max']

# Extract all chemistry information from the silo files into a chemistry container
# This uses the first time instant's silo file to initialize
pion.load_chemistry()
print(pion.chemistry_container)

cooling = nebula.cooling(
    database='/home/tony/Desktop/NebulaPy/NebulaPy-DB',
    pion_ion=ion,
    verbose=True
)

ion_name = ion.replace('+', 'p')
ion_output_dir = os.path.join(output_dir, ion_name)
os.makedirs(ion_output_dir, exist_ok=True)
runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time

    print(f" ---------------------------")
    # Print the current simulation time instant
    sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    # Extract necessary physical parameters for the current time instant
    density = pion.get_parameter('Density', silo_instant)  # Retrieve density
    temperature = pion.get_parameter('Temperature', silo_instant)  # Retrieve temperature
    # Identify the shape of each density array to ensure compatibility with other parameters.
    shape_list = [arr.shape for arr in density]

    # Initialize arrays for electron number density (ne) and mass fraction sum,
    # with zeroes matching the shape of the density data.
    cooling_rate_map = [np.zeros(shape) for shape in shape_list]

    # Retrieve the electron number density
    ne = pion.get_ne(silo_instant)

    print(" computing cooling rate map for each grid levels")
    for level in range(N_grid_level):
        level_cooling_rate_map = cooling.generate_cooling_rate_map(
            temperature=temperature[level], ne=ne[level])

        cooling_rate_map[level] = level_cooling_rate_map


    print(" generating cooling rate map")
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.text(0.05, 0.9, 'time = %5.2f kyr' % sim_time.value, transform=ax.transAxes,
            fontsize=12, color='white')
    ax.set_xlim(mesh_edges_min[0][0].value, mesh_edges_max[0][0].value)
    ax.set_ylim(mesh_edges_min[0][1].value, mesh_edges_max[0][1].value)

    # Assuming parameter_data, dims_min_1, dims_max_1, and cmap_value are defined
    # Assuming ax is already defined
    for level in range(N_grid_level):
        plot_data = np.log10(cooling_map[level])
        extents = [mesh_edges_min[level][0].value, mesh_edges_max[level][0].value,
                   mesh_edges_min[level][1].value, mesh_edges_max[level][1].value]

        image = ax.imshow(plot_data, interpolation='nearest', cmap='inferno',
                          extent=extents, origin='lower',
                          #vmin=limits[0], vmax=limits[1]
                          )

    # Create divider for existing axes instance
    divider = make_axes_locatable(ax)
    # Append axes to the right of ax1
    cax = divider.append_axes("right", size="5%", pad=0.05)

    # Add colorbar with appropriate ticks
    colorbar = plt.colorbar(image, cax=cax, ticks=MultipleLocator(1))

    # Set the formatter for the colorbar
    colorbar.ax.yaxis.set_major_formatter(ScalarFormatter())
    colorbar.ax.yaxis.get_major_formatter().set_scientific(False)  # Disable scientific notation
    colorbar.ax.yaxis.get_major_formatter().set_useOffset(False)  # Disable offset

    # Update the ticks for the colorbar
    colorbar.update_ticks()

    ax.set_xlabel('z (pc)', fontsize=12)
    ax.set_ylabel('R (pc)', fontsize=12)

    # Correctly add text to the plot

    ax.text(0.65, 0.9, ion, transform=ax.transAxes, fontsize=12, color='white')  # Ensure a string is provided here
    ax.tick_params(axis='both', which='major', labelsize=13)
    # ax.axes.get_xaxis().set_visible(False)  # Remove the x-axis

    # fig.subplots_adjust(wspace=0, hspace=0)  # Remove the whitespace between the images
    filename = f"{filebase}_coolmap_{ion_name}_{str(step).zfill(4)}.png"
    filepath = os.path.join(ion_output_dir, filename)
    plt.savefig(filepath, bbox_inches="tight", dpi=300)
    plt.close(fig)

    print(f' saving cooling rate map to {filename}')




    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")
