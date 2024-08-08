import glob
import os
import re
from pypion.ReadData import ReadData
from pypion.SiloHeader_data import OpenData
from NebulaPy.tools import util as util
from NebulaPy.tools import constants as const


class pion_silos():
    '''
    This class is not an alternative to Pypion; rather, it is a
    bundle of methods useful for creating synthetic emission
    maps from the Silo file.
    '''


    def __init__(self, verbose):
        self.verbose = verbose
        self.chemistry_container = {}


    ######################################################################################
    # batch silo files
    ######################################################################################
    def batch_silos(self, dir, filebase):
        '''
        Find silo files and put them in order of levels.
        :param dir: Directory path where files are located
        :param filebase: Base name of the files to search for
        :return: List of batched silo files by levels
        '''
        all_silos = []
        batched_silos = []

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
            pass

        return batched_silos




    ######################################################################################
    # get chemistry from the initial silo file
    ######################################################################################
    def get_chemistry(self, instant_silo_set):
        '''
        This method extracts information related to the chemistry and chemical tracers,
        transforming the chemical tracer names to a format that PyPion can directly
        read from the Silo file. This method can be included in the next version of
        PyPion.

        Parameters
        ----------
        instant_silo_set : The instance for which chemical data is to be extracted

        Returns
        -------
        Generates and stores the following in self.chemistry_container:
        - 'dynamics': Dynamics data retrieved from the Silo file
        - 'chemistry': Chemistry flag indicating if chemistry data is available
        - 'E_update': Energy update information (if chemistry flag is true)
        - 'chemistry_code': The code indicating the type of chemistry (if chemistry flag is true)
        - 'microphysics': List of microphysics processes (if chemistry flag is true)
        - 'Ntracers': Number of chemical tracers
        - 'mpv10_elements': List of elements identified for MPv10 chemistry code
        - 'mpv10_tracers': List of tracers corresponding to each element for MPv10 chemistry code
        '''

        # Read data for the given instant silo set
        data = ReadData(instant_silo_set)
        # Open the data for the given instant silo set
        header_data = OpenData(instant_silo_set)
        # Set the directory to '/header'
        header_data.db.SetDir('/header')
        # Retrieve the value of "EP_chemistry" from the header data
        chemistry_flag = header_data.db.GetVar("EP_chemistry")
        # Retrieve the value of "EP_dynamics" from the header data
        dynamics = header_data.db.GetVar("EP_dynamics")

        # save the dynamics and chemistry_flag values in the chemistry_container dictionary
        self.chemistry_container['dynamics'] = dynamics
        self.chemistry_container['chemistry'] = chemistry_flag

        # Define the list of process variable names
        processes = ['EP_coll_ionisation', 'EP_rad_recombination',
                     'EP_cooling', 'EP_raytracing', 'EP_phot_ionisation',
                     'EP_charge_exchange']

        # Define the list of process names corresponding to the process variable names
        processes_name = ['coll_ionisation', 'rad_recombination', 'cooling',
                          'raytracing', 'phot_ionisation', 'charge_exchange']

        # Check if chemistry_flag is true
        if chemistry_flag:
            # Retrieve the value of "EP_update_erg"
            energy_update = header_data.db.GetVar("EP_update_erg")
            # save the energy_update value in the chemistry_container dictionary
            self.chemistry_container['E_update'] = energy_update
            # Retrieve the value of "chem_code"
            chemistry_code = header_data.db.GetVar("chem_code")[0]
            # save the chemistry_code value in the chemistry_container dictionary
            self.chemistry_container['chemistry_code'] = chemistry_code

            # Initialize an empty list to store microphysics processes
            microphysics = []
            # Check if the chemistry_code is not 'MPv10'
            if not chemistry_code == 'MPv10':
                # Exit with an error if the chemistry_code is not 'MPv10'
                util.nebula_exit_with_error("PION is not running MPv10; NelubaPy functionality is limited.")
            else:
                # If verbose is enabled, print the chemistry code
                if self.verbose:
                    print(f" pion chemistry: {chemistry_code}")

                # Loop through each process
                for index, process in enumerate(processes):
                    # Check if the process variable exists in the header data
                    if header_data.db.GetVar(process):
                        # Append the corresponding process name to the microphysics list
                        microphysics.append(processes_name[index])

                # save the microphysics list in the chemistry_container dictionary
                self.chemistry_container['microphysics'] = microphysics
                # Retrieve the number of tracers
                Ntracers = header_data.db.GetVar('num_tracer')
                # initialize an empty list to store elements for MPv10
                mpv10elements = []
                # initialize a 2D list to store tracers for each element
                mpv10tracers = [[] for _ in range(len(const.nebula_elements))]
                # list of element names from the nebula_elements dictionary keys
                element_list = list(const.nebula_elements.keys())
                # save the number of tracers in the chemistry_container dictionary
                self.chemistry_container['Ntracers'] = Ntracers
                # If verbose is enabled, print the number of chemical tracers
                if self.verbose:
                    print(f" number of chemical tracers: {Ntracers}")

                # Loop through each tracer index
                for i in range(Ntracers):
                    # create a tracer index string with leading zeros
                    tracer_index = f'Tracer{i:03}'
                    # retrieve the tracer value
                    chem_tracer = header_data.db.GetVar(tracer_index)[0]

                    # check if the tracer is an element ('X' denoting elemental mass fraction)
                    if 'X' in chem_tracer and chem_tracer.replace("_", "").replace("X", "") in const.nebula_elements:
                        # extract the element name
                        element = chem_tracer.replace("_", "").replace("X", "")
                        # get the full element name from the nebula_elements dictionary
                        element_name = const.nebula_elements[element]
                        # append the element to the mpv10elements list
                        mpv10elements.append(element)
                        # if verbose is enabled, print the found element name
                        if self.verbose:
                            print(f" found {element_name}")

                        # get the index of the element in the element_list
                        element_index = element_list.index(element)
                        # append the tracer with the corresponding element to the mpv10tracers list
                        if 0 <= element_index < len(mpv10tracers):
                            mpv10tracers[element_index].append(f'Tr{i:03}_' + chem_tracer)

                    # check if the tracer is a corresponding ion
                    if re.sub(r'\d{1,2}\+', '', chem_tracer) in const.nebula_elements:
                        # extract the element name
                        element = re.sub(r'\d{1,2}\+', '', chem_tracer)
                        # get the index of the element in the element_list
                        element_index = element_list.index(element)
                        # gppend the tracer with the corresponding ion to the mpv10tracers list
                        mpv10tracers[element_index].append(f'Tr{i:03}_' + chem_tracer)

                # save the mpv10elements and mpv10tracers lists in the chemistry_container dictionary
                self.chemistry_container['mpv10_elements'] = mpv10elements
                self.chemistry_container['mpv10_tracers'] = mpv10tracers



