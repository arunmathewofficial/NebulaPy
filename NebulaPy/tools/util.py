import sys
from NebulaPy import version
import os
import glob

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
    This method searches for silo files with the same file base in a specified
    directory and organizes them based on their corresponding time instants.
    If the simulation is run on a nested grid with multiple levels, then for each
    time instant, an inner array is created that contains the silo files for
    that instant, with each file in the inner array representing a different
    grid level.
    :param dir: Directory path where files are located
    :param filebase: Base name of the files to search for
    :return: List of batched silo files by levels
    '''
    all_silos = []
    batched_silos = []
    print(f" ---------------------------")
    print(" batching silo files into time instances")

    # Iterate through possible levels to find the files
    for level in range(20):
        search_pattern = os.path.join(dir, f"{filebase}_level{str(level).zfill(2)}_0000.*.silo")
        level_files = sorted(glob.glob(search_pattern))

        if level_files:
            # Ensure all files for the level have the same extension pattern
            file_extensions = [os.path.splitext(file)[1] for file in level_files]
            unique_extensions = set(file_extensions)

            if len(unique_extensions) == 1:
                all_silos.append(level_files)
            else:
                nebula_exit_with_error(f"level {level} file missing for file extensions {unique_extensions}")
                break
        else:
            break

    if not all_silos:
        nebula_exit_with_error(f"No '{filebase}' silo files found in '{dir}'")

    try:
        for i in range(len(all_silos[0])):
            instantaneous_set = [all_silos[level][i] for level in range(len(all_silos))]
            batched_silos.append(instantaneous_set)
    except IndexError:
        print(f" batch processing {level} level silo file completed")
        pass

    del all_silos
    return batched_silos

#*************************************************************************************
