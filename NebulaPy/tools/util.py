import sys
from NebulaPy import version
from pypion.SiloHeader_data import OpenData
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
# batch silo files
######################################################################################
def batch_silos(dir, filebase):
    '''
    Organizes silo files into groups based on their corresponding time instants.
    This function scans a specified directory for silo files with a given base name,
    and if the simulation involves multiple grid levels, it groups the files by time instant.
    For each time instant, it creates a sublist containing silo files from different grid levels.

    :param dir: Directory path where the silo files are located.
    :param filebase: Base name of the silo files to search for.
    :return: A list of lists, where each sublist contains silo files for a specific time instant,
             with each file in the sublist representing a different grid level.
    :raises: RuntimeError if no silo files are found or if expected files for certain levels are missing.
    '''

    print(" ---------------------------")
    print(" batching silo files into time instances")

    # Construct search pattern for silo files
    search_pattern = os.path.join(dir, f"{filebase}_*.silo")
    all_silos = glob.glob(search_pattern)

    if not all_silos:
        raise RuntimeError(f"no '{filebase}' silo files found in '{dir}'")

    # Assume the number of grid levels; this may need to be adjusted based on actual usage
    first_silo = all_silos[0]
    # Open the data for the first silo instant silo
    header_data = OpenData(first_silo)
    # Set the directory to '/header'
    header_data.db.SetDir('/header')
    # Retrieve no of nested grid levels
    Nlevels = header_data.db.GetVar("grid_nlevels")
    # close the object
    header_data.close()

    if Nlevels == 1:
        print(" batching completed")
        return sorted(all_silos)

    else:
        # Pattern to match level 00 files
        pattern = re.compile(r'_level00_0000\.\d+\.silo$')
        # Find and sort level 00 files, one per time instant
        level00_instants = [file for file in all_silos if pattern.search(file)]
        batched_silos = [[file] for file in sorted(level00_instants)]

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
                    nebula_exit_with_error(f"missing file for level {level} and instant {instant_extension}")
                else:
                    # Append the found file to the corresponding time instant group
                    batched_silos[i].append(level_instant_file)

        print(" batching completed")
        return batched_silos
#*************************************************************************************
