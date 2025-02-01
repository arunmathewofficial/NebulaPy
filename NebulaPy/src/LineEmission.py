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
        """
        Computes the most luminous emission lines for a given species in a 1D or 2D dataset.

        This method calculates the line luminosity for a given ionized species across multiple grid levels.
        It retrieves the brightest emission lines by computing the line emissivity for each cell and summing
        the contributions across the grid.

        Parameters
        ----------
        temperature : list (1D or 2D array-like)
            Temperature values of the grid cells.
        ne : list (1D or 2D array-like)
            Electron density values corresponding to each grid cell.
        species_density : list (1D or 2D array-like)
            Density of the given ionized species in each grid cell.
        cell_volume : list (1D or 2D array-like)
            Volume of each grid cell.
        grid_mask : list (1D or 2D array-like)
            Mask specifying active grid cells.
        Nlines : int
            The number of most luminous emission lines to retrieve.

        Returns
        -------
        dict
            A dictionary containing:
            - 'spectroscopic' : str, the spectroscopic name of the species.
            - 'lines' : list, the N most luminous emission line wavelengths.
            - 'luminosity' : list, the luminosity values of the brightest lines.

        Notes
        -----
        - A small electron density tolerance (`electron_tolerance = 1.E-08`) is set to avoid division by zero.
        - The method loops over each grid level and computes line luminosities.
        - The most luminous lines are selected using `np.argsort()`.
        """

        # Define a tolerance for electron density (to avoid division by zero)
        electron_tolerance = 1.E-08

        # Get the number of grid levels in the temperature dataset
        NGlevel = len(temperature)

        # Handle cylindrical coordinates (assumes 2D data)
        rows = len(temperature[0])  # Number of rows in the temperature array

        # Retrieve the list of possible emission lines for the species
        dummy_temperature_array = [1000]
        dummy_ne_array = [1.0]
        species = chianti(pion_ion=self.ion, temperature=dummy_temperature_array, ne=dummy_ne_array, verbose=False)
        spectroscopic_name = species.chianti_ion.Spectroscopic

        # Check if the species has emission lines
        if 'line' not in species.species_attributes_container[species.chianti_ion_name]['keys']:
            util.nebula_warning(f"{spectroscopic_name} has no line emission associated")
            return {'spectroscopic': spectroscopic_name}
        else:
            all_lines = species.get_line_emissivity(allLines=False)['wvl']
        del species

        # Initialize an array to store total line luminosities
        species_all_line_luminosity = np.zeros_like(all_lines)

        # Loop over each grid level
        for level in range(NGlevel):
            print(f" computing luminosity for all {spectroscopic_name} lines at grid level {level}", end='\r')

            # Ensure electron density values are nonzero
            ne[level] = np.array(ne[level])
            ne[level][ne[level] == 0] = electron_tolerance

            # Initialize luminosity storage for this level
            species_all_lines_luminosity_level = np.zeros_like(all_lines)

            # Loop over each row (assumes a 2D dataset)
            for row in range(rows):
                temperature_row = temperature[level][row]
                ne_row = ne[level][row]
                species_density_row = species_density[level][row]
                cell_volume_row = cell_volume[level][row]
                grid_mask_row = grid_mask[level][row]

                # Compute emissivity for the species at the given conditions
                species = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
                all_lines_emissivity_info_row = species.get_line_emissivity(allLines=False)
                del species

                all_lines_emissivity_row = all_lines_emissivity_info_row['emiss']

                # Compute total luminosity for each emission line
                for index in range(len(all_lines)):
                    species_all_lines_luminosity_level[index] += (
                            4.0 * const.pi * np.sum(all_lines_emissivity_row[index]
                                                    * species_density_row
                                                    * cell_volume_row
                                                    * grid_mask_row)
                    )
                del all_lines_emissivity_info_row

            # Accumulate luminosity across all grid levels
            species_all_line_luminosity += species_all_lines_luminosity_level

        print(f" completed the luminosity computation for all {spectroscopic_name} lines", end='\n')

        # Retrieve the N most luminous lines
        print(f" retrieving the {Nlines} most luminous {spectroscopic_name} lines", end='\n')

        indices = np.argsort(species_all_line_luminosity)[-Nlines:]  # Get indices of the brightest lines
        brightest_lines_luminosity = species_all_line_luminosity[indices]
        brightest_lines = all_lines[indices]

        return {
            'spectroscopic': spectroscopic_name,
            'lines': brightest_lines,
            'luminosity': brightest_lines_luminosity
        }

    ######################################################################################
    # line luminosity for a given list of lines in cylindrical coordinate system
    ######################################################################################

    def line_luminosity_cylindrical(self, lines, temperature, ne, species_density, cell_volume, grid_mask, Nlines):

        # Define a tolerance for electron density (to avoid division by zero)
        electron_tolerance = 1.E-08

        # Get the number of grid levels in the temperature dataset
        NGlevel = len(temperature)

        # Handle cylindrical coordinates (assumes 2D data)
        rows = len(temperature[0])  # Number of rows in the temperature array

        # identify all lines
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



        # Initialize an array to store total line luminosities
        species_all_line_luminosity = np.zeros_like(all_lines)

        # Loop over each grid level
        for level in range(NGlevel):
            print(f" computing luminosity for all {spectroscopic_name} lines at grid level {level}", end='\r')

            # Ensure electron density values are nonzero
            ne[level] = np.array(ne[level])
            ne[level][ne[level] == 0] = electron_tolerance

            # Initialize luminosity storage for this level
            species_all_lines_luminosity_level = np.zeros_like(all_lines)

            # Loop over each row (assumes a 2D dataset)
            for row in range(rows):
                temperature_row = temperature[level][row]
                ne_row = ne[level][row]
                species_density_row = species_density[level][row]
                cell_volume_row = cell_volume[level][row]
                grid_mask_row = grid_mask[level][row]

                # Compute emissivity for the species at the given conditions
                species = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
                all_lines_emissivity_info_row = species.get_line_emissivity(allLines=False)
                del species

                all_lines_emissivity_row = all_lines_emissivity_info_row['emiss']

                # Compute total luminosity for each emission line
                for index in range(len(all_lines)):
                    species_all_lines_luminosity_level[index] += (
                            4.0 * const.pi * np.sum(all_lines_emissivity_row[index]
                                                    * species_density_row
                                                    * cell_volume_row
                                                    * grid_mask_row)
                    )
                del all_lines_emissivity_info_row

            # Accumulate luminosity across all grid levels
            species_all_line_luminosity += species_all_lines_luminosity_level

        print(f" completed the luminosity computation for all {spectroscopic_name} lines", end='\n')

        # Retrieve the N most luminous lines
        print(f" retrieving the {Nlines} most luminous {spectroscopic_name} lines", end='\n')

        indices = np.argsort(species_all_line_luminosity)[-Nlines:]  # Get indices of the brightest lines
        brightest_lines_luminosity = species_all_line_luminosity[indices]
        brightest_lines = all_lines[indices]




