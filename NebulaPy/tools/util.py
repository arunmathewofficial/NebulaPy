import sys

import numpy as np

from NebulaPy import version
from pypion.SiloHeader_data import OpenData
from pypion.ReadData import ReadData
import astropy.units as unit
import os
import glob
import re


######################################################################################
# Nebulapy exit with error
######################################################################################
def nebula_exit_with_error(errorMessage):
    """
    Custom exit function display error before exiting.
    :param message: Optional exit message.
    """
    RED = "\033[91m"
    RESET = "\033[0m"
    print(f'{RED} error: {errorMessage}{RESET}')
    print(f' NebulaPy {version.__version__} exiting ...')
    sys.exit()

######################################################################################
# Nebula Warning
######################################################################################
def nebula_warning(warnMessage):
    """
    Custom exit function display error before exiting.
    :param message: Optional exit message.
    """
    YELLOW = "\033[93m"
    RESET = "\033[0m"
    print(f'{YELLOW} warning: {warnMessage}{RESET}')

######################################################################################
# Nebula version
######################################################################################
def nebula_version():
    """
    return NebulaPy version
    :param message: Optional exit message.
    """
    return f'NebulaPy {version.__version__}'

######################################################################################
# batch silo files
######################################################################################
def batch_silos(dir, filebase, start_time=None, finish_time=None, time_unit=None,
                out_frequency=None):
    '''
    Organizes silo files into groups based on their corresponding time instants.
    This function scans a specified directory for silo files with a given base name,
    and if the simulation involves multiple grid levels, it groups the files by time instant.
    For each time instant, it creates a sublist containing silo files from different grid levels.

    The function allows running from a particular start time to a finish time, enabling users to
    process specific intervals of simulation data. Additionally, one can restart the processing
    by adjusting the start time, making it flexible for iterative analysis or simulations.

    :param dir: Directory path where the silo files are located.
    :param filebase: Base name of the silo files to search for.
    :param start_time: The start time for processing, in the given time_unit.
    :param finish_time: The finish time for processing, in the given time_unit.
    :param time_unit: The time unit for start_time and finish_time ('sec', 'yr', or 'kyr').
    :return: A list of lists, where each sublist contains silo files for a specific time instant,
             with each file in the sublist representing a different grid level.
    :raises: RuntimeError if no silo files are found or if expected files for certain levels are missing.
    '''
    conversion_factors = {
        'kyr': 3.154e+10,
        'yr': 3.154e+7,
        'sec': 1
    }

    # Convert start_time and finish_time to seconds based on time_unit
    factor = conversion_factors.get(time_unit, 1)  # Default to 1 if time_unit is None or unrecognized
    start_time_sec = start_time * factor if start_time is not None else None
    finish_time_sec = finish_time * factor if finish_time is not None else None

    print(" ---------------------------")
    print(" batching silo files into time instances:")
    print(f" starting time: {start_time} {time_unit}")
    print(f" finishing time: {finish_time} {time_unit}")
    print(f" output frequency: {out_frequency}")

    # Construct search pattern for silo files
    search_pattern = os.path.join(dir, f"{filebase}_*.silo")
    all_silos = glob.glob(search_pattern)

    if not all_silos:
        nebula_exit_with_error(f"no '{filebase}' silo files found in '{dir}'")

    # Assume the number of grid levels; this may need to be adjusted based on actual usage
    # Open header data
    header_data = OpenData(all_silos)
    # Set the directory to '/header'
    header_data.db.SetDir('/header')
    # Retrieve what coordinate system is used
    coord_sys = header_data.db.GetVar("coord_sys")
    # Retrieve no of nested grid levels
    Nlevels = header_data.db.GetVar("grid_nlevels")
    # close the object
    header_data.close()
    # set number of time instance
    Ninstances = 0

    # for uniform grid *********************************************************************
    if Nlevels == 1:
        print(f" grid: uniform")
        silos = sorted(all_silos)
        batched_silos = [[silo] for silo in silos]
        selected_silos = []
        for i, silo in enumerate(batched_silos):
            data = ReadData(silo)
            basic = data.get_1Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()

            if (start_time_sec is not None and sim_time < start_time_sec) or \
                    (finish_time_sec is not None and sim_time > finish_time_sec):
                continue  # Skip this silo if it doesn't meet the criteria
            selected_silos.append(silo)

        batched_silos = selected_silos
        # Keep files based on the output frequency if specified
        if out_frequency is not None:
            Ninstances = len(batched_silos)
            # indices to keep: multiples of out_freq, plus the first and last index
            #indices_to_keep = sorted(set(range(0, Ninstances, out_frequency)) | {0, Ninstances - 1})
            indices_to_keep = sorted(set(range(0, Ninstances, out_frequency)))
            batched_silos = [batched_silos[i] for i in indices_to_keep]

        Ninstances = len(batched_silos)
        print(f" {Ninstances} time instances between {start_time} {time_unit} and {finish_time} {time_unit}")
        print(" batching completed")
        return batched_silos

    # for nested grid *********************************************************************
    else:
        print(f" grid: nested with {Nlevels} levels")
        # Pattern to match level 00 files
        pattern = re.compile(r'_level00_0000\.\d+\.silo$')
        # Find and sort level 00 files, one per time instant
        level00_instants = [file for file in all_silos if pattern.search(file)]
        batched_silos = [[file] for file in sorted(level00_instants)]

        # check if level 00 batched silos are within the asked time range
        selected_silos = []
        for i, silo in enumerate(batched_silos):
            data = ReadData(silo)
            if coord_sys == 3:
                basic = data.get_1Darray('Density')
            elif coord_sys == 2:
                basic = data.get_2Darray('Density')
            elif coord_sys == 1:
                basic = data.get_3Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()

            if (start_time_sec is not None and sim_time < start_time_sec) or \
                    (finish_time_sec is not None and sim_time > finish_time_sec):
                continue  # Skip this silo if it doesn't meet the criteria
            selected_silos.append(silo)

        batched_silos = selected_silos

        # Keep files based on the output frequency if specified
        if out_frequency is not None:
            Ninstances = len(batched_silos)
            # indices to keep: multiples of out_freq, plus the first and last index
            indices_to_keep = sorted(set(range(0, Ninstances, out_frequency)) | {0, Ninstances - 1})
            batched_silos = [batched_silos[i] for i in indices_to_keep]

        # appending other level silos to corresponding time instant
        for i, instant in enumerate(batched_silos):
            # Extract the time instant from the level 00 file
            file_name = instant[0].split('/')[-1]
            instant_extension = file_name.split('_')[-1].replace('.silo', '')

            # Find and append silo files for the same time instant across other levels
            for level in range(1, Nlevels):
                level_pattern = re.compile(f'level{str(level).zfill(2)}_{instant_extension}.silo')

                # Find the file for the current level
                level_instant_file = next((file for file in all_silos if level_pattern.search(file)), None)

                if level_instant_file is None:
                    nebula_exit_with_error(f"missing silo file for level {level} instant {instant_extension}")
                else:
                    # Append the found file to the corresponding time instant group
                    batched_silos[i].append(level_instant_file)

        Ninstances = len(batched_silos)
        print(f" {Ninstances} time instances between {start_time} {time_unit} and {finish_time} {time_unit}")
        print(" batching completed")
        if not batched_silos:
            nebula_exit_with_error('no silo files found in the specified time range, check your selection criteria')
        return batched_silos

######################################################################################
# get element symbol for a specific ion
######################################################################################
def get_element_symbol(ion):
    """
    Extracts the element symbol from an ion identifier by filtering out non-alphabetic characters.

    Parameters:
    ion (str): The ion identifier, which may include both alphabetic and non-alphabetic characters
               (e.g., 'He1+', 'O2+', 'C5+').

    Returns:
    str: The element symbol extracted from the ion identifier (e.g., 'He' from 'He1+', 'O' from 'O2+').
    """
    return ''.join(filter(str.isalpha, ion))



######################################################################################
# general progress bar
######################################################################################
def progress_bar(step, Nstep, prefix='', suffix='', condition=False, completion_msg=''):
    """
    Displays a progress bar in the terminal.

    Args:
    - step (int): Current step (0-indexed). Represents the current progress point.
    - Nstep (int): Total number of steps. Defines the total length of the process.
    - prefix (str): Optional prefix text for the progress bar. Allows user to add custom text at the start.
    - suffix (str): Optional suffix text for the progress bar. Allows user to add custom text at the end.
    - condition (bool): A flag to indicate whether the completion message should be shown once the progress reaches 100%. If `True`, the `completion_msg` is displayed at the end.
    - completion_msg (str): Message to display upon completion. This message appears when the progress reaches 100%.
    """

    step += 1  # Increment step by 1 since it starts from 0, for clearer progress display.
    length = 20  # Length of the progress bar (20 characters long)
    percent = ("{0:.1f}").format(100 * (step / float(Nstep)))  # Percentage progress
    filled_length = int(length * step // Nstep)  # Number of filled segments in the progress bar
    bar = '█' * filled_length + '-' * (length - filled_length)  # Progress bar representation with filled and unfilled segments

    # Format the progress bar string for terminal output
    progress_line = f' {prefix} |{bar}| {percent:>5}% {suffix} '

    if step == 1:
        # Print initial progress line without overwrite
        sys.stdout.write(progress_line)
        sys.stdout.flush()

    # Handle completion
    elif step == Nstep and condition:
        # When the process reaches the final step and the condition is True,
        # clear the progress bar and display a completion message
        sys.stdout.write('\r' + ' ' * len(progress_line) + ' ' * 5 + '\r')  # Clear the current progress line
        sys.stdout.write(f' {completion_msg}\n')  # Print the completion message
        sys.stdout.flush()

    else:
        # Regular progress update by overwriting the previous line
        sys.stdout.write('\r' + progress_line + '\r')
        sys.stdout.flush()






