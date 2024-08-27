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
    RED = "\033[91m"
    RESET = "\033[0m"
    print(f' error: {RED}{warnMessage}{RESET}')
    sys.exit()

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
# file check
######################################################################################
def file_check(file):
    """
    return NebulaPy version
    :param message: Optional exit message.
    """


######################################################################################
# batch silo files
######################################################################################
def batch_silos(dir, filebase, start_time=None, finish_time=None, time_unit=None, out_freq=None):
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
    print(" batching silo files into time instances")

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
        for i, silo in enumerate(batched_silos):
            data = ReadData(silo)
            basic = data.get_1Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()

            if start_time_sec is not None and sim_time < start_time_sec:
                # Remove the current silo from the list
                batched_silos.pop(i)
            elif finish_time_sec is not None and sim_time > finish_time_sec:
                # Remove the current silo and all remaining silos from the list
                batched_silos = batched_silos[:i]
                break

        # Keep files based on the output frequency if specified
        if out_freq is not None:
            Ninstances = len(batched_silos)
            # indices to keep: multiples of out_freq, plus the first and last index
            indices_to_keep = sorted(set(range(0, Ninstances, out_freq)) | {0, Ninstances - 1})
            batched_silos = [batched_silos[i] for i in indices_to_keep]

        Ninstances = len(batched_silos)
        print(f" {Ninstances} instances found")
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
        for i, silo in enumerate(batched_silos):
            data = ReadData(silo)
            basic = data.get_1Darray('Density')
            sim_time = (basic['sim_time'] * unit.s).value
            data.close()
            if start_time_sec is not None and sim_time < start_time_sec:
                # Remove the current silo from the list
                batched_silos.pop(i)
            elif finish_time_sec is not None and sim_time > finish_time_sec:
                # Remove the current silo and all remaining silos from the list
                batched_silos = batched_silos[:i]
                break

        # Keep files based on the output frequency if specified
        if out_freq is not None:
            Ninstances = len(batched_silos)
            # indices to keep: multiples of out_freq, plus the first and last index
            indices_to_keep = sorted(set(range(0, Ninstances, out_freq)) | {0, Ninstances - 1})
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
                    nebula_exit_with_error(f"missing file for level {level} instant {instant_extension}")
                else:
                    # Append the found file to the corresponding time instant group
                    batched_silos[i].append(level_instant_file)

        Ninstances = len(batched_silos)
        print(f" {Ninstances} instances found")
        print(" batching completed")
        return batched_silos
#*************************************************************************************
