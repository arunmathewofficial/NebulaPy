from .Chianti import chianti
import NebulaPy.tools.constants as const
import numpy as np
import NebulaPy.tools.util as util

class line_emission():

    ######################################################################################
    # initializing the class emissionline
    ######################################################################################
    def __init__(self, ion, verbose):
        """
        only single ion is considered here
        """
        self.ion = ion
        self.verbose = verbose
        self.line_emission_container = {}
        self.line_emission_container['ion'] = self.ion


    ######################################################################################
    # line luminosity
    ######################################################################################
    def lineluminosity_spherical(self, lines, temperature, ne, species_density, shell_volume):
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

        all_emissivity_data = ion.get_emissivity()
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
    # line luminosity
    ######################################################################################
    def line_emissivity_cylindrical(self, lines, temperature, ne, geometry_container={}):

        coordinate_sys = geometry_container['coordinate_sys']  # Coordinate system

        if not coordinate_sys == 'cylindrical':
            if self.verbose:
                util.nebula_exit_with_error(f" coordinate system is not cylindrical")

        N_grid_level = geometry_container['Nlevels']  # Number of grid levels
        # Convert line wavelengths to string format for later use
        lines_str = [str(line) for line in lines]

        # Define a tolerance for electron density (avoid division by zero)
        electron_tolerance = 1.E-08

        rows = len(temperature[0])  # Number of inner arrays in 2D temperature data

        # Initialize a list to store emissivity maps
        emissivity_map = [np.zeros(shape) for shape in [arr.shape for arr in temperature]]
        emissivity_map_dict = {line: emissivity_map for line in lines_str}  # Map emission lines to maps

        # Loop over each grid level in the simulation
        for level in range(N_grid_level):
            ne[level] = np.array(ne[level])
            ne[level][ne[level] == 0] = electron_tolerance  # Replace zero values with tolerance

            # Loop over each row in the grid
            for row in range(rows):
                print(f" computing emissivity for line(s) at level {level}, row {row}", end='\r')

                # Get temperature and electron density for the current row
                temperature_row = temperature[level][row]
                ne_row = ne[level][row]

                # Create a Chianti object for line emissivity calculation
                chinati = chianti(temperature_row, ne_row, pion_ion=self.ion, verbose=False)

                # Calculate the line emissivity for each specified line
                lines_emissivity_row = chinati.get_line_emissivity(line_list=lines)

                # Store the computed emissivity for the current line
                for line in lines_str:
                    emissivity_map_dict[line][level][row] = lines_emissivity_row[line]

            print(f" returning emissivity map for line(s) at level {level}")

        return emissivity_map_dict



