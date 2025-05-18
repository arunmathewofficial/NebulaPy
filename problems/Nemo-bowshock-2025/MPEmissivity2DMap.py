"""
Spectral Line Emissivity Mapping

Description:
This script computes spatial emissivity maps for selected spectral emission lines
of various ions in an astrophysical simulation. It reads a sequence of simulation
snapshots (silo files), extracts temperature and electron number density, and
calculates the emissivity of each line using atomic data. The computation is
parallelized using multiprocessing, where each worker handles the emission
calculation for one ion.

Supported ions include:
    He1+, C2+, N1+, N2+, O1+, O2+, Ne1+, Ne2+, S1+, S2+, and S3+
Each ion has a user-defined list of spectral lines (wavelengths in Angstroms).

Features:
- Reads simulation data from DB silo files using the `NebulaPy` interface.
- Loads and processes simulation geometry and chemistry information.
- Extracts physical parameters such as temperature and electron density.
- Computes line emissivity maps for selected lines using a cylindrical projection.
- Parallel processing for emissivity calculations using Python's multiprocessing.
- Organizes results in ion-specific subdirectories under the given output path.

Note:
- The script assumes the NebulaPy module is correctly installed and configured.
- The specific ions to process can be toggled on/off in the `ions` dictionary.
- Output is organized by ion and line, suitable for further post-processing or plotting.

Author: Arun Mathew
Date: 18 May 2025
"""
import os
import warnings
from NebulaPy.tools import util
import NebulaPy.src as nebula
import numpy as np
import multiprocessing as mp
import time
import h5py
import matplotlib.pyplot as plt  # Matplotlib for creating plots
from mpl_toolkits.axes_grid1 import make_axes_locatable  # For adding colorbars to plots
from matplotlib.ticker import MultipleLocator, ScalarFormatter  # For formatting plot ticks and colorbar

# Suppress specific warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="divide by zero encountered in log10")

# Set OMP_NUM_THREADS=1 to ensure single-threaded execution
os.environ["OMP_NUM_THREADS"] = "1"

def compute_emissivity(workerQ, doneQ):
    """Worker function to compute emissivity for each lines of ions"""
    import queue  # Import inside function to avoid errors

    while True:
        try:
            task = workerQ.get(timeout=timeout)  # Get task from queue
            if task is None:
                # print("[Worker] Received stop signal, exiting.")
                break  # Exit if None is received (stop signal)

            ion_name, line_emission, lines, temperature, ne = task

            print(f" multiprocessing: computing line emissivity of {ion_name:<4} lines")

            emissivity = line_emission.line_emissivity_map_cylindrical(
                lines=lines, temperature=temperature, ne=ne, progress_bar=False)

            print(f" multiprocessing: finished computing emissivity for {ion_name:<4} lines")
            doneQ.put({ion_name: emissivity})  # Store result
        except queue.Empty:  # Use correct exception for empty queue
            util.nebula_exit_with_error(f"multiprocessing - no task in queue")
            break  # Queue is empty, exit loop
        except Exception as e:
            util.nebula_exit_with_error(f'multiprocessing - working process {e}')
            break


if __name__ == "__main__":

    # Input-output file configuration for the high-resolution simulation on MIMIR
    # output_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/high-res/emissivity'
    # silo_dir = '/mnt/massive-stars/data/arun_simulations/Nemo_BowShock/high-res/silo'
    # filebase = 'Ostar_mhd-nemo-dep_d2n0128l3'  # Base name of the silo files
    # start_time = None
    # finish_time = None
    # out_frequency = 2

    # Input-output file configuration for the low-resolution simulation on Razer Blade machine.
    output_dir = '/home/tony/Desktop/multi-ion-bowshock/sim-output/emissiviity_map'
    silo_dir = '/home/tony/Desktop/multi-ion-bowshock/sim-output/silo'
    filebase = 'Ostar_mhd-nemo_d2n0128l3'  # Base name of the silo files
    start_time = None
    finish_time = None
    out_frequency = 2

    # Batch the silo files for analysis within the specified time range
    batched_silos = util.batch_silos(
        silo_dir,
        filebase,
        start_time=start_time,
        finish_time=finish_time,
        time_unit='kyr',
        out_frequency=out_frequency
    )

    # Total number of time instant
    N_time_instant = len(batched_silos)

    # Initialize the PION class to handle simulation data
    pion = nebula.pion(batched_silos, verbose=True)

    # Calculates and stores geometric grid parameters.
    # For example, in a spherical geometry, it extracts radius and shell volumes
    # from the first silo file in the batch and saves them into a geometry container.
    pion.load_geometry(scale='pc')
    N_grid_level = pion.geometry_container['Nlevel']
    mesh_edges_min = pion.geometry_container['edges_min']
    mesh_edges_max = pion.geometry_container['edges_max']

    # h5 header and metadata
    heading = f"Generated by {util.nebula_version()}"
    metadata = {"minimum_mesh_edges": mesh_edges_min,
                "maximum_mesh_edges": mesh_edges_max,
                'Nlevels': N_grid_level}

    # Extract all chemistry information from the silo files into a chemistry container
    # This uses the first time instant's silo file to initialize
    pion.load_chemistry()

    print(f" ---------------------------")
    print(f" task: multiprocessing for emissivity map of spectral lines")

    # Set up the ion and line emission parameters
    # Define the ions and their respective emission lines
    He1P_pion_ion = 'He1+'
    He1P_lines = [303.78, 303.786, 256.317, 256.318, 243.026, 243.027]  # Emission line(s) of interest
    print(rf" {He1P_pion_ion:<4} lines: {', '.join(map(str, He1P_lines))}  Angstrom")
    He1P_pion_ion_name = He1P_pion_ion.replace('+', 'p')
    He1P_ion_output_dir = os.path.join(output_dir, He1P_pion_ion_name)
    os.makedirs(He1P_ion_output_dir, exist_ok=True)

    C2P_pion_ion = 'C2+'
    C2P_lines = [1906.683, 1908.734, 977.02]  # Emission line(s) of interest
    print(rf" {C2P_pion_ion:<4} lines: {', '.join(map(str, C2P_lines))}  Angstrom")
    C2P_pion_ion_name = C2P_pion_ion.replace('+', 'p')
    C2P_ion_output_dir = os.path.join(output_dir, C2P_pion_ion_name)
    os.makedirs(C2P_ion_output_dir, exist_ok=True)

    N1P_pion_ion = 'N1+'
    N1P_lines = [6585.273, 6549.861, 1218026.8, 2053388.09]
    print(rf" {N1P_pion_ion:<4} lines: {', '.join(map(str, N1P_lines))}  Angstrom")
    N1P_pion_ion_name = N1P_pion_ion.replace('+', 'p')
    N1P_ion_output_dir = os.path.join(output_dir, N1P_pion_ion_name)
    os.makedirs(N1P_ion_output_dir, exist_ok=True)

    N2P_pion_ion = 'N2+'
    N2P_lines = [573394.5, 989.799, 1752.16, 1749.674, 1753.995, 1748.646]
    print(rf" {N2P_pion_ion:<4} lines: {', '.join(map(str, N2P_lines))}  Angstrom")
    N2P_pion_ion_name = N2P_pion_ion.replace('+', 'p')
    N2P_ion_output_dir = os.path.join(output_dir, N2P_pion_ion_name)
    os.makedirs(N2P_ion_output_dir, exist_ok=True)

    O1P_pion_ion = 'O1+'  # The ion of interest (Oxygen II)
    O1P_lines = [3729.844, 3727.092, 7331.722, 2470.97, 2471.094, 7321.094, 7322.177, 7332.808]
    print(rf" {O1P_pion_ion:<4} lines: {', '.join(map(str, O1P_lines))} Angstrom")
    O1P_pion_ion_name = O1P_pion_ion.replace('+', 'p')
    O1P_ion_output_dir = os.path.join(output_dir, O1P_pion_ion_name)
    os.makedirs(O1P_ion_output_dir, exist_ok=True)

    O2P_pion_ion = 'O2+'  # The ion of interest (Oxygen III)
    O2P_lines = [5008.24, 883564.0, 518145.0, 4960.295, 1666.15, 4364.436, 832.929, 833.715, 1660.809]
    print(rf" {O2P_pion_ion:<4} lines: {', '.join(map(str, O2P_lines))} Angstrom")
    O2P_pion_ion_name = O2P_pion_ion.replace('+', 'p')
    O2P_ion_output_dir = os.path.join(output_dir, O2P_pion_ion_name)
    os.makedirs(O2P_ion_output_dir, exist_ok=True)

    Ne1P_pion_ion = 'Ne1+'
    Ne1P_lines = [128139.42]
    print(rf" {Ne1P_pion_ion:<4} lines: {', '.join(map(str, Ne1P_lines))} Angstrom")
    Ne1P_pion_ion_name = Ne1P_pion_ion.replace('+', 'p')
    Ne1P_ion_output_dir = os.path.join(output_dir, Ne1P_pion_ion_name)
    os.makedirs(Ne1P_ion_output_dir, exist_ok=True)

    Ne2P_pion_ion = 'Ne2+'
    Ne2P_lines = [155545.19, 3869.849, 3968.585, 360230.55, 489.495, 491.041]
    print(rf" {Ne2P_pion_ion:<4} lines: {', '.join(map(str, Ne2P_lines))} Angstrom")
    Ne2P_pion_ion_name = Ne2P_pion_ion.replace('+', 'p')
    Ne2P_ion_output_dir = os.path.join(output_dir, Ne2P_pion_ion_name)
    os.makedirs(Ne2P_ion_output_dir, exist_ok=True)

    S1P_pion_ion = 'S1+'
    S1P_lines = [6718.295, 6732.674, 4069.749, 10323.317, 4077.5]
    print(rf" {S1P_pion_ion:<4} lines: {', '.join(map(str, S1P_lines))} Angstrom")
    S1P_pion_ion_name = S1P_pion_ion.replace('+', 'p')
    S1P_ion_output_dir = os.path.join(output_dir, S1P_pion_ion_name)
    os.makedirs(S1P_ion_output_dir, exist_ok=True)

    S2P_pion_ion = 'S2+'
    S2P_lines = [335008.38, 9532.252, 187055.74, 9070.048, 6313.649, 3722.454]
    print(rf" {S2P_pion_ion:<4} lines: {', '.join(map(str, S2P_lines))} Angstrom")
    S2P_pion_ion_name = S2P_pion_ion.replace('+', 'p')
    S2P_ion_output_dir = os.path.join(output_dir, S2P_pion_ion_name)
    os.makedirs(S2P_ion_output_dir, exist_ok=True)

    S3P_pion_ion = 'S3+'
    S3P_lines = [105104.95]
    print(rf" {S3P_pion_ion:<4} lines: {', '.join(map(str, S3P_lines))} Angstrom")
    S3P_pion_ion_name = S3P_pion_ion.replace('+', 'p')
    S3P_ion_output_dir = os.path.join(output_dir, S3P_pion_ion_name)
    os.makedirs(S3P_ion_output_dir, exist_ok=True)

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

    # Loop over each time instant in the batched silo files
    runtime = 0.0
    write_heading = True
    for step, silo_instant in enumerate(batched_silos):
        silo_instant_start_time = time.time()

        print(f" ---------------------------")
        sim_time = pion.get_simulation_time(silo_instant, time_unit='kyr')
        print(f" step: {step}/{N_time_instant - 1} | simulation time: {sim_time:.6e}")

        # Extract temperature and electron number density
        temperature = pion.get_parameter('Temperature', silo_instant)
        ne = pion.get_ne(silo_instant)

        # ions for processing
        ions = {
            "He1+": (He1P_line_emission, He1P_lines),
            "C2+": (C2P_line_emission, C2P_lines),
            "N1+": (N1P_line_emission, N1P_lines),
            "N2+": (N2P_line_emission, N2P_lines),
            "O1+": (O1P_line_emission, O1P_lines),
            "O2+": (O2P_line_emission, O2P_lines),
            "Ne1+": (Ne1P_line_emission, Ne1P_lines),
            "Ne2+": (Ne2P_line_emission, Ne2P_lines),
            "S1+": (S1P_line_emission, S1P_lines),
            "S2+": (S2P_line_emission, S2P_lines),
            "S3+": (S3P_line_emission, S3P_lines)
        }

        # get keys of ions
        keys = ions.keys()
        # dict to save all outputs of a single time instant
        species_lines_emissivity = {}

        # Parameters for multiprocessing
        timeout = 0.1
        Ncores = len(ions)
        cpu_count = mp.cpu_count()
        proc = min(Ncores, cpu_count)
        print(f" multiprocessing: utilizing {proc} cores (one core per ion)")

        with mp.Manager() as manager:
            emissivity_workerQ = manager.Queue()
            emissivity_doneQ = manager.Queue()

            print(" multiprocessing: adding tasks to queue")
            # Populate worker queue
            for ion_name, (line_emission, lines) in ions.items():
                emissivity_workerQ.put((ion_name, line_emission, lines, temperature, ne))

            print(" multiprocessing: starting worker processes")
            # Start worker processes
            emissivity_processes = []
            for _ in range(proc):
                p = mp.Process(target=compute_emissivity, args=(emissivity_workerQ, emissivity_doneQ))
                p.start()
                emissivity_processes.append(p)

            # Allow some time for processing
            time.sleep(2)

            # Send stop signal (None) to workers
            for _ in range(proc):
                emissivity_workerQ.put(None)

            # Wait for processes to complete
            for p in emissivity_processes:
                p.join()

            # Collect results
            while not emissivity_doneQ.empty():
                species_lines_emissivity.update(emissivity_doneQ.get())

        # plotting emissivity map and saving emissivity data into files
        for ion in species_lines_emissivity:

            ion_name = ion.replace('+', 'p')
            ion_output_dir = os.path.join(output_dir, ion_name)

            # saving data to h5 file ########################################################
            h5_filename = f"{filebase}_emiss_{ion_name}_{str(step).zfill(4)}.h5"
            h5_file = os.path.join(ion_output_dir, h5_filename)
            data_title = f"Bow-Shock emissivity map for {ion}"
            ion_emissivity_map_dict = species_lines_emissivity[ion]

            print(f" saving {ion:>5} line emissivity maps data into {h5_filename}")
            with h5py.File(h5_file, "w") as file:
                # Add heading and title
                file.attrs['head'] = heading
                file.attrs['title'] = data_title

                # Save metadata
                meta = file.create_group("metadata")
                for key, value in metadata.items():
                    if isinstance(value, list) or isinstance(value, np.ndarray):
                        meta.create_dataset(key, data=np.array(value))
                    else:
                        meta.attrs[key] = value

                # Save emissivity function map
                emissivity_map_group = file.create_group("emissivity_map")
                for key, value in ion_emissivity_map_dict.items():
                    emissivity_map_group.create_dataset(str(key), data=np.array(value, dtype=np.float32))
            # end of saving data to h5 file #################################################

            # Generating image
            # Loop through each emission line and plot the corresponding emissivity map
            for line in ion_emissivity_map_dict:
                line_name = line.replace(' ', '')

                fig, ax = plt.subplots(figsize=(8, 6))  # Create a new figure for the plot

                # Add text annotation for simulation time
                ax.text(0.05, 0.9, 'time = %5.2f kyr' % sim_time.value, transform=ax.transAxes,
                        fontsize=12, color='white')

                # Set plot limits based on mesh edges
                ax.set_xlim(mesh_edges_min[0][0].value, mesh_edges_max[0][0].value)
                ax.set_ylim(mesh_edges_min[0][1].value, mesh_edges_max[0][1].value)

                # Plot the emissivity map using a logarithmic scale
                for level in range(N_grid_level):
                    plot_data = np.log10(ion_emissivity_map_dict[line][level])  # Logarithmic emissivity data
                    extents = [mesh_edges_min[level][0].value, mesh_edges_max[level][0].value,
                               mesh_edges_min[level][1].value, mesh_edges_max[level][1].value]

                    # Create the image plot
                    image = ax.imshow(plot_data, interpolation='nearest', cmap='inferno', extent=extents,
                                      origin='lower',
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
                ax.text(0.65, 0.9, line, transform=ax.transAxes, fontsize=12, color='white')

                # Customize tick labels
                ax.tick_params(axis='both', which='major', labelsize=13)

                # Save the figure to a file
                image_filename = f"{filebase}_emiss_{line_name}_{sim_time.value:.2f}kyr.png"
                filepath = os.path.join(ion_output_dir, image_filename)
                print(f" saving {line:>14} emissivity map to {image_filename}")

                plt.savefig(filepath, bbox_inches="tight", dpi=300)
                plt.close(fig)
            del ion_emissivity_map_dict

        del species_lines_emissivity

        # Update runtime
        silo_instant_finish_time = time.time()
        dt = silo_instant_finish_time - silo_instant_start_time
        runtime += dt
        print(f" runtime: {runtime:.4e} s | Δt: {dt:.4e} s")
