from .Chianti import chianti
import NebulaPy.tools.constants as const
import numpy as np
import NebulaPy.tools.util as util


import multiprocessing
from multiprocessing import Pool, Manager


class line_emission():

    ######################################################################################
    # initializing the class emissionline
    ######################################################################################
    def __init__(self, ion, verbose=True):
        """
        only single ion is considered here
        """
        self.ion = ion
        self.verbose = verbose
        self.line_emission_container = {}
        self.line_emission_container['ion'] = self.ion


    ######################################################################################
    # line luminosity in 1D spherical setting
    ######################################################################################
    def line_luminosity_spherical(self, lines, temperature, ne, species_density, shell_volume):
        '''

        Parameters
        ----------
        lines
        temperature
        ne
        ns
        dV

        Returns
        -------

        '''

        ion = chianti(pion_ion=self.ion, temperature=temperature, ne=ne, verbose=self.verbose)
        self.line_emission_container['temperature'] = temperature
        self.line_emission_container['ne'] = ne
        self.line_emission_container['spectroscopic'] = ion.chianti_ion.Spectroscopic

        # if the line (wavelength) is given in string, get the corresponding
        # float value
        #if isinstance(line, str):
        #    line = const.wvl_dict[line]

        all_emissivity_data = ion.get_line_emissivity()
        allLines = all_emissivity_data['wvl']
        self.line_emission_container['allLines'] = allLines
        self.line_emission_container['lines'] = lines

        indices = []
        for line in lines:
            if self.verbose:
                print(f' identifying {line} Å from allLines of {ion.chianti_ion.Spectroscopic}')
            index = (np.abs(allLines - line)).argmin()
            tolerance = 10 ** -4
            if np.abs(allLines[index] - line) <= tolerance:
                if self.verbose:
                    print(f' line {line} Å found at index {index} in allLines')
                indices.append(index)
            else:
                util.nebula_exit_with_error('line not found in allLines')

        self.line_emission_container['line_indices'] = indices

        if self.verbose:
            print(f' retrieving emissivity values for {ion.chianti_ion.Spectroscopic} lines(s): {lines}')

        emissivity = np.asarray(all_emissivity_data['emiss'][indices])
        self.line_emission_container['emiss'] = emissivity

        # Calculating line Luminosity
        if self.verbose:
            print(f' calculating line luminosity for {ion.chianti_ion.Spectroscopic} lines(s): {lines}')
        luminosity = [4.0 * const.pi * np.sum(e * species_density * shell_volume) for e in emissivity]
        self.line_emission_container['luminosity'] = luminosity

    ######################################################################################
    # line emissivity map for a given list of lines in cylindrical coordinate system
    ######################################################################################
    def line_emissivity_map_cylindrical(self, lines, temperature, ne):

        # Define a tolerance for electron density (to avoid division by zero)
        electron_tolerance = 1.E-08
        # Get number of grid levels in the temperature data
        NGlevel = len(temperature)
        # Number of inner arrays in 2D temperature data
        rows = len(temperature[0])

        # Initialize a list to store emissivity maps
        emissivity_map = [np.zeros(shape) for shape in [arr.shape for arr in temperature]]
        # Convert line wavelengths to string format for later use
        lines_str = [str(line) for line in lines]
        emissivity_map_dict = {line: emissivity_map for line in lines_str}  # Map emission lines to maps

        # Loop over each grid level in the simulation
        for level in range(NGlevel):
            ne[level] = np.array(ne[level])
            ne[level][ne[level] == 0] = electron_tolerance  # Replace zero values with tolerance

            # Loop over each row in the grid
            for row in range(rows):
                print(f" computing emissivity for line(s) at level {level}, row {row}", end='\r')

                # Get temperature and electron density for the current row
                temperature_row = temperature[level][row]
                ne_row = ne[level][row]

                # Create a Chianti object for line emissivity calculation
                ion = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
                # Calculate the line emissivity for each specified line
                lines_emissivity_row = ion.get_line_emissivity_for_list(line_list=lines)

            # Store the computed emissivity for the current line
            for line in lines_str:
                emissivity_map_dict[line][level][row] = lines_emissivity_row[line]

        print(f" returning emissivity map for line(s) for all level(s)       ")

        return emissivity_map_dict


    ######################################################################################
    # get dominant list of lines for a simulation snapshot in cylindrical coordinate system
    ######################################################################################
    def get_species_dominant_lines(self, temperature, ne, species_density, cell_volume, grid_mask, Nlines):

        '''
        This is a general method for 1d and 2d dataset
        Parameters
        ----------
        temperature
        ne
        cell_volume
        Nlines

        Returns
        -------

        '''

        # Define a tolerance for electron density (to avoid division by zero)
        electron_tolerance = 1.E-08
        # Get number of grid levels in the temperature data
        NGlevel = len(temperature)


        # for cylindrical coordinate ######################################################
        # Number of inner arrays in 2D temperature data
        rows = len(temperature[0])
        # get all lines of the species
        dummy_temperature_array = [1000]
        dummy_ne_array = [1.0]
        species = chianti(pion_ion=self.ion, temperature=dummy_temperature_array, ne=dummy_ne_array, verbose=False)
        spectroscopic_name = species.chianti_ion.Spectroscopic
        if 'line' not in species.species_attributes_container[species.chianti_ion_name]['keys']:
            util.nebula_warning(f"{spectroscopic_name} has no line emission associate")
            return {'spectroscopic': spectroscopic_name}
        else:
            all_lines = species.get_line_emissivity(allLines=False)['wvl']
        del species

        species_all_line_luminosity = np.zeros_like(all_lines)
        # Looping over grid levels
        for level in range(NGlevel):
            print(f" computing luminosity for all {spectroscopic_name} lines at grid level {level}", end='\r')
            ne[level] = np.array(ne[level])
            ne[level][ne[level] == 0] = electron_tolerance

            # Initialize species_line_luminosity outside the row loop
            species_all_lines_luminosity_level = np.zeros_like(all_lines)

            for row in range(rows):
                temperature_row = temperature[level][row]
                ne_row = ne[level][row]
                species_density_row = species_density[level][row]
                cell_volume_row = cell_volume[level][row]
                grid_mask_row = grid_mask[level][row]

                species = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
                all_lines_emissivity_info_row = species.get_line_emissivity(allLines=False)
                del species

                all_lines_emissivity_row = all_lines_emissivity_info_row['emiss']

                for index in range(len(all_lines)):
                    species_all_lines_luminosity_level[index] += 4.0 * const.pi \
                                                                 * np.sum(all_lines_emissivity_row[index]
                                                                          * species_density_row
                                                                          * cell_volume_row
                                                                          * grid_mask_row)
                del all_lines_emissivity_info_row

            species_all_line_luminosity += species_all_lines_luminosity_level

        # start coding here
        print(f" completed the luminosity computation for all {spectroscopic_name} lines", end='\n')
        # find 10 brighest lines from the list
        print(f" retrieving the {Nlines} most luminous {spectroscopic_name} lines", end='\n')

        indices = np.argsort(species_all_line_luminosity)[-Nlines:]
        brightest_lines_luminosity = species_all_line_luminosity[indices]
        brightest_lines = all_lines[indices]

        return {'spectroscopic': spectroscopic_name,'lines': brightest_lines,
                'luminosity': brightest_lines_luminosity}





        

















    ######################################################################################
    # line luminosity for a given list of lines in cylindrical coordinate system
    ######################################################################################

    def line_luminosity_cylindrical(self, lines, temperature, ne):
        pass



