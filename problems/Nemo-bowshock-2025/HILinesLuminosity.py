"""
Line Luminosity

Description:
This script calculates the temporal evolution of luminosity for selected spectral lines
in an astrophysical simulation. The lines correspond to different ions,
including He1+, C2+, N1+, N2+, O1+, O2+, Ne1+, Ne2+, S1+, S2+, and S3+, and their
corresponding emission lines. The luminosity for each line is computed using the
simulation data and stored for later analysis.

Features:
- Reads silo files generated from the simulation to extract temperature, electron
density, and ion number densities.
- Computes the line luminosity for specific emission lines over a time range.
- Outputs the results to a text file, with each row representing a time step and the luminosities of different lines.

Notes:
Calculating line luminosity only from the shocked layer.  For the line luminosities. For this we exclude any gas that
satisfies the following:
(1) wind tracer value <0.5,
(2) density <1.1x the background inflow density from the parameter file,
(3) v_x within a few (5%) percent of the inflow value (-35 km/s?)

Author: Arun Mathew
Date: 01 Feb 2025
"""

import os
import time
import warnings
from NebulaPy.tools import util
import NebulaPy.src as nebula
import multiprocessing as mp
import time
import numpy as np
import queue  # Import queue globally

# remove this after coding shocked ISM mask
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator, ScalarFormatter

#


# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")

# Set OMP_NUM_THREADS=1 to ensure single-threaded execution
os.environ["OMP_NUM_THREADS"] = "1"


def compute_luminosity(workerQ, doneQ):
    """Worker function to compute luminosity for each ion"""
    import queue  # Import inside function to avoid errors

    while True:
        try:
            task = workerQ.get(timeout=timeout)  # Get task from queue
            if task is None:
                # print("[Worker] Received stop signal, exiting.")
                break  # Exit if None is received (stop signal)

            ion_name, line_emission, lines, species_density, temperature, ne, cell_volume, grid_mask = task

            print(f" multiprocessing: computing the luminosity of {ion_name:<4} lines")

            luminosity = line_emission.line_luminosity_2D(
                lines=lines,
                temperature=temperature, ne=ne,
                species_density=species_density,
                cell_volume=cell_volume,
                grid_mask=grid_mask, progress_bar=False)

            print(f" multiprocessing: finished computing the luminosity for {ion_name:<4} lines")
            doneQ.put({ion_name: luminosity})  # Store result
        except queue.Empty:  # Use correct exception for empty queue
            util.nebula_exit_with_error(f"multiprocessing - no task in queue")
            break  # Queue is empty, exit loop
        except Exception as e:
            util.nebula_exit_with_error(f'multiprocessing - working process {e}')
            break


def generate_shocked_ism_mask(pion, silo_instant):
    """
    Generate a mask for the shocked ISM region based on specific criteria:
    (1) wind tracer value < 0.5,
    (2) density < 1.1x the background inflow density from the parameter file,
    (3) v_x within a few (5%) percent of the inflow value (-35 km/s)

    Parameters:
    pion: nebula.pion object
    silo_instant: the current silo file being processed.

    Returns:
    numpy.ndarray
        A boolean mask array where True indicates shocked ISM regions.
    """

    print(f" generating shocked ISM mask")

    # Extract necessary parameters from the simulation
    density = pion.get_parameter('Density', silo_instant)
    grid_mask = pion.get_parameter('NG_Mask', silo_instant)
    shocked_ism_mask = [np.ones(shape) for shape in [arr.shape for arr in density]]
    wind_tracer = pion.get_parameter('Tr000_WIND', silo_instant)
    velocity_Y = pion.get_parameter('VelocityY', silo_instant)

    for level in range(N_grid_level):
        shocked_ism_mask[level][wind_tracer[level] > 0.5] = 0
        shocked_ism_mask[level][density[level] < 1.05 * 1.0E-23] = 0
        shocked_ism_mask[level][velocity_Y[level] <= 12.0E+03] = 0
        shocked_ism_mask[level] = shocked_ism_mask[level] * grid_mask[level]

    return shocked_ism_mask


if __name__ == "__main__":

    # Input-output file configuration for the low-resolution simulation on Razer Blade machine.
    output_dir = '/home/tony/Desktop/multi-ion-bowshock/sim-output/time-lines-luminosity'
    silo_dir = '/home/tony/Desktop/multi-ion-bowshock/high-res-silos-200kyr'
    filebase = 'Ostar_mhd-nemo-dep_d2n0384l3'  # Base name of the silo files
    filename = filebase + '_HI_lines_luminosity_HighRes.txt'
    start_time = None
    finish_time = None
    out_frequency = None
    time_unit = 'kyr'


    # Input-output file configuration for the medium-resolution simulation on MIMIR
    # output_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/med-res/time-lines-luminosity'
    # silo_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/med-res/silo'
    # filebase = 'Ostar_mhd-nemo-dep_d2n0192l3'  # Base name of the silo files
    # filename = filebase + '_lines_luminosity_MedRes_3.txt'
    # start_time = 0
    # finish_time = 47
    # out_frequency = 2

    # Input-output file configuration for the high-resolution simulation on MIMIR
    #output_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/high-res/time-lines-luminosity'
    #silo_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/high-res/silo'
    #filebase = 'Ostar_mhd-nemo-dep_d2n0384l3'  # Base name of the silo files
    #filename = filebase + '_lines_luminosity_HighRes_ShockedISM.txt'
    #start_time = None
    #finish_time = 201
    #out_frequency = 2
    #time_unit = 'kyr'




    # creating output file
    outfile = os.path.join(output_dir, filename)

    # Batch the silo files for analysis within the specified time range
    batched_silos = util.batch_silos(
        silo_dir,
        filebase,
        start_time=start_time,
        finish_time=finish_time,
        time_unit=time_unit,
        out_frequency=out_frequency
    )

    # Total number of time instant
    N_time_instant = len(batched_silos)

    # Initialize the PION class to handle simulation data
    pion = nebula.pion(batched_silos, verbose=True)
    # Load chemistry and geometry data
    pion.load_chemistry()
    pion.load_geometry(scale='cm')

    print(f" ---------------------------")
    print(f" task: multiprocessing for computing the temporal evolution of Hα, Hβ luminosities")
    print(" list of spectral lines to be processed:")

    # Set up the ion and line emission parameters
    # Define the ions and their respective emission lines
    # Chianti H-alpha lines 6564.538, 6564.564, 6564.523, 6564.665, 6564.722 (Å)
    # and the repeated entries near 6564 Å correspond to the fine-structure
    # components of the hydrogen Hα line (Balmer α, 3 → 2 transition).

    # collisional lines
    pion_ion = 'H'
    Chianti_HAlphaLines = [6564.538, 6564.564, 6564.523, 6564.665, 6564.722]  # in Angstroms
    print(f" Hα collisional de-excited chianti lines: {', '.join(map(str, Chianti_HAlphaLines))}  \u212B")
    Chianti_HBetaLines = [4862.733, 4862.72]  # in Angstroms
    print(f" Hβ collisional de-excited chianti lines: {', '.join(map(str, Chianti_HBetaLines))}  \u212B")

    # recombination lines
    pyneb_ion = 'H1+'
    PyNeb_HAlphaLines = [6562.816]
    print(f" Hα recombination pyneb lines:: {', '.join(map(str, PyNeb_HAlphaLines))}  \u212B")
    PyNeb_HBetaLines = [4861.332]
    print(f" Hβ recombination pyneb lines: {', '.join(map(str, PyNeb_HBetaLines))}  \u212B")


    # Initialize the emission line calculations for each ion
    H_line_emission = nebula.line_emission(pion_ion, verbose=True)
    H1P_line_emission = nebula.line_emission(pyneb_ion, verbose=True)

    # Check the requested lines in the database for each ion
    print(f" ---------------------------")
    print(f" checking requested lines in CHIANTI database:")
    H_line_emission.chianti_line_batch_check(Chianti_HAlphaLines)
    H_line_emission.chianti_line_batch_check(Chianti_HBetaLines)

    print(f" ---------------------------")
    print(f" checking requested lines in PyNeb database:")
    H1P_line_emission.pyneb_line_batch_check(PyNeb_HAlphaLines)
    H1P_line_emission.pyneb_line_batch_check(PyNeb_HBetaLines)

    # Prepare output file for results
    outfile = os.path.join(output_dir, filename)

    # Get geometry information
    geometry = pion.geometry_container
    N_grid_level = geometry['Nlevel']
    grid_mask = geometry['mask']
    cell_volume = pion.get_cylindrical_cell_volume().value

    mesh_edges_min = pion.geometry_container['edges_min']
    mesh_edges_max = pion.geometry_container['edges_max']

    # Write initial header to the output file
    with open(outfile, "w") as file:
        file.write(f"#File generated by {util.nebula_version()}\n\n")
        file.write("#Task: Multiprocessing for Computing the Temporal Evolution of Luminosity in Spectral Line\n\n")
        file.write(f"#PION Simulation Reference: NEMO Bowshock 2025\n")
        file.write(f"#PION Simulation FileBase: {filebase}\n\n")
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
        print(f" step: {step}/{N_time_instant - 1} | simulation time: {sim_time:.6e}")

        # get temperature
        temperature = pion.get_parameter('Temperature', silo_instant)
        # calculate electron number density
        ne = pion.get_ne(silo_instant)
        # generate ISM shocked mask
        shocked_ism_mask = generate_shocked_ism_mask(pion, silo_instant)

        # Retrieve species number density
        H_num_density = pion.get_ion_number_density(pion_ion, silo_instant)
        H1P_num_density = pion.get_ion_number_density(pyneb_ion, silo_instant)

        # ions for processing
        ions = {
            "H": (H_line_emission, Chianti_HAlphaLines, H_num_density),
            "H1+": (H1P_line_emission, PyNeb_HAlphaLines, H1P_num_density),
            "H": (H_line_emission, Chianti_HBetaLines, H_num_density),
            "H1+": (H1P_line_emission, PyNeb_HBetaLines, H1P_num_density),
        }

        # get keys of ions
        keys = ions.keys()
        species_lines_luminosity = {}

        # Parameters for multiprocessing
        timeout = 0.1
        Ncores = len(ions)
        cpu_count = mp.cpu_count()
        proc = min(Ncores, cpu_count)
        print(f" multiprocessing: utilizing {proc} cores (one core per ion)")

        with mp.Manager() as manager:
            luminosity_workerQ = manager.Queue()
            luminosity_doneQ = manager.Queue()

            print(" multiprocessing: adding tasks to queue")
            # Populate worker queue
            for ion_name, (line_emission, lines, species_density) in ions.items():
                luminosity_workerQ.put(
                    (ion_name, line_emission, lines, species_density, temperature, ne, cell_volume, shocked_ism_mask))

            print(" multiprocessing: starting worker processes")
            # Start worker processes
            luminosity_processes = []
            for _ in range(proc):
                p = mp.Process(target=compute_luminosity, args=(luminosity_workerQ, luminosity_doneQ))
                p.start()
                luminosity_processes.append(p)

            # Allow some time for processing
            time.sleep(2)

            # Send stop signal (None) to workers
            for _ in range(proc):
                luminosity_workerQ.put(None)

            # Wait for processes to complete
            for p in luminosity_processes:
                p.join()

            # Collect results
            while not luminosity_doneQ.empty():
                species_lines_luminosity.update(luminosity_doneQ.get())

        # sort the dictionary according to the species dictionary
        sorted_species_lines_luminosity = {key: species_lines_luminosity[key] for key in keys}
        # combine all sub dictionary into a single dictionary
        species_lines_luminosity_dict = {key: value for subdict
                                         in sorted_species_lines_luminosity.values() for key, value in subdict.items()}

        print(f" saving data to file: {filename}")
        if write_heading:
            with open(outfile, "a") as file:
                # Write the time (or any other desired variable)
                file.write(f"Time ")
                # Write the keys from Dictionaries
                file.write(" ".join(f"{key}" for key in species_lines_luminosity_dict.keys()))

                file.write("\n")
            write_heading = False

        with open(outfile, "a") as file:
            file.write(f"{sim_time.value:.6e} ")
            file.write(" ".join(f"{value:.6e}" for value in species_lines_luminosity_dict.values()))
            file.write("\n")

        del species_lines_luminosity_dict

        '''
        # generating masked grid image
        fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

        ax.text(0.05, 0.9, 'time = %5.2f kyr' % sim_time.value, transform=ax.transAxes,
                fontsize=12, color='white')
        ax.set_xlim(mesh_edges_min[0][0].value, mesh_edges_max[0][0].value)
        ax.set_ylim(mesh_edges_min[0][1].value, mesh_edges_max[0][1].value)

        # Assuming parameter_data, dims_min_1, dims_max_1, and cmap_value are defined
        # Assuming ax is already defined
        for level in range(N_grid_level):
            plot_data = shocked_ism_mask[level]
            extents = [mesh_edges_min[level][0].value, mesh_edges_max[level][0].value,
                       mesh_edges_min[level][1].value, mesh_edges_max[level][1].value]

            image = ax.imshow(plot_data, interpolation='nearest', cmap='inferno',
                              extent=extents, origin='lower')

        # Create divider for existing axes instance
        divider = make_axes_locatable(ax)
        # Append axes to the right of ax1
        cax = divider.append_axes("right", size="5%", pad=0.05)

        # Add colorbar with appropriate ticks
        colorbar = plt.colorbar(image, cax=cax,
                                #ticks=MultipleLocator(1)
                                )

        # Set the formatter for the colorbar
        colorbar.ax.yaxis.set_major_formatter(ScalarFormatter())
        colorbar.ax.yaxis.get_major_formatter().set_scientific(False)  # Disable scientific notation
        colorbar.ax.yaxis.get_major_formatter().set_useOffset(False)  # Disable offset

        # Update the ticks for the colorbar
        colorbar.update_ticks()
        ax.set_xlabel('z (pc)', fontsize=12)
        ax.set_ylabel('R (pc)', fontsize=12)

        ax.text(0.65, 0.9, 'Shocked ISM', transform=ax.transAxes, fontsize=12, color='white')  # Ensure a string is provided here
        ax.tick_params(axis='both', which='major', labelsize=13)
        # ax.axes.get_xaxis().set_visible(False)  # Remove the x-axis
        plt.show()
        '''

        # Update runtime
        silo_instant_finish_time = time.time()
        dt = silo_instant_finish_time - silo_instant_start_time
        runtime += dt
        print(f" runtime: {runtime:.4e} s | step runtime: {dt:.4e} s")

