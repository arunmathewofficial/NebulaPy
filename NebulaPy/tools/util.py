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
    Find silo files and put them in order of levels.
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
            all_silos.append(level_files)
        else:
            break

    if not all_silos:
        print(f" error: No '{filebase}' silo files found in '{dir}'.")
        quit()

    try:
        for i in range(len(all_silos[0])):
            instantaneous_set = [all_silos[level][i] for level in range(len(all_silos))]
            batched_silos.append(instantaneous_set)
    except IndexError:
        print(f" batch processing {level} level silo file completed")
        pass

    return batched_silos
#*************************************************************************************
