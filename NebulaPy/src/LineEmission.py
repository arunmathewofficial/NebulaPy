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
    # line emissivity in cylindrical coordinate system
    ######################################################################################
    def line_emissivity_cylindrical(self, lines, temperature, ne, geometry_container):

        # get coordinate system
        coordinate_sys = geometry_container['coordinate_sys']
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
                ion = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
                # Calculate the line emissivity for each specified line
                 lines_emissivity_row = ion.get_line_emissivity(line_list=lines)

                # Store the computed emissivity for the current line
                for line in lines_str:
                    emissivity_map_dict[line][level][row] = lines_emissivity_row[line]

            print(f" returning emissivity map for line(s) at level {level}       ")
        return emissivity_map_dict


    ######################################################################################
    # get dominant lines of a ion for give temperature and electron number density
    ######################################################################################
    def get_dominant_lines(self, temperature, ne, Nlines, geometry_container):

        coordinate_sys = geometry_container['coordinate_sys']
        # Define a tolerance for electron density (avoid division by zero)
        electron_tolerance = 1.E-08

        dominant_lines_container = {
            'ion': self.ion,
            'electron_tolerance': electron_tolerance,
        }

        ###########################################################################################
        # For cylindrical coordinate system
        if coordinate_sys == 'cylindrical':
            N_grid_level = geometry_container['Nlevels']  # Number of grid levels
            rows = len(temperature[0])  # Number of inner arrays in 2D temperature data

            wvls = [np.zeros((Nlines,)) for _ in range(N_grid_level)]
            largest_emissivity = [np.zeros((Nlines,)) for _ in range(N_grid_level)]

            # all emission lines of the ion
            all_wvls = []

            # Loop over each grid level in the simulation
            for level in range(N_grid_level):
                ne[level] = np.array(ne[level])
                ne[level][ne[level] == 0] = electron_tolerance  # Replace zero values with tolerance

                previous_largest_emissivity = np.zeros(Nlines)
                previous_largest_emissivity_index = np.zeros(Nlines, dtype=int)

                # Loop over each row in the grid
                for row in range(rows):
                    print(
                        f" retrieving emissivity of all lines of {self.ion} at level {level}, row {row}  ",
                        end='\r'
                    )
                    # Get temperature and electron density for the current row
                    temperature_row = temperature[level][row]
                    ne_row = ne[level][row]

                    # Create a Chianti object for line emissivity calculation
                    ion = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
                    chianti_ion_name = ion.chianti_ion_name

                    if 'line' in ion.species_attributes_container[chianti_ion_name]['keys']:
                        # Get the emissivity for wavelengths
                        emissivity_data = ion.get_emissivity()
                        emissivity_for_wvl = emissivity_data['emiss']
                        all_wvls = np.array(emissivity_data['wvl'])
                        emissivity_for_wvl_1D = np.max(emissivity_for_wvl, axis=1)
                        del ion.species_attributes_container

                        # Find the indices of the Nlines largest values
                        current_largest_emissivity_index = np.argsort(emissivity_for_wvl_1D)[-Nlines:]
                        current_largest_emissivity = emissivity_for_wvl_1D[current_largest_emissivity_index]

                        # Update previous arrays elements using current array element
                        for i in range(Nlines):
                            if current_largest_emissivity[i] > previous_largest_emissivity[i]:
                                previous_largest_emissivity[i] = current_largest_emissivity[i]
                                previous_largest_emissivity_index[i] = current_largest_emissivity_index[i]

                        wvls[level] = all_wvls[previous_largest_emissivity_index]
                        largest_emissivity[level] = previous_largest_emissivity

                    else:
                        if self.verbose:
                            print(f" grid level {level}: {self.ion} has no line emission, skipping ..." + " " * 10)
                        break

            dominant_lines_container['emiss'] = largest_emissivity
            dominant_lines_container['wvls'] = wvls
            print(f" completed identifying the dominant lines for ion {self.ion} in the current silo instance")
        else:
            if self.verbose:
                util.nebula_exit_with_error(" Coordinate system is not cylindrical")
        # End of cylindrical section ##########################################################

        return dominant_lines_container
    ###########################################################################################













    ######################################################################################
    # multi-processing line emissivity calculation
    ######################################################################################
    # Refactored function to minimize CPU usage
    def compute_row_emissivity(self, args):

        row, level, temperature_row, ne_row, lines, lines_str = args
        ion = chianti(pion_ion=self.ion, temperature=temperature_row, ne=ne_row, verbose=False)
        lines_emissivity_row = ion.get_line_emissivity(line_list=lines)
        return {line: lines_emissivity_row[line] for line in lines_str}

    def MP_line_emissivity_cylindrical(self, lines, temperature, ne, geometry_container, ncores=3):
        cpu_count = multiprocessing.cpu_count()
        if ncores >= cpu_count:
            raise ValueError("Using excessive cores for multiprocessing")

        if geometry_container['coordinate_sys'] != 'cylindrical':
            raise ValueError("Coordinate system is not cylindrical")

        N_grid_level = geometry_container['Nlevels']
        lines_str = [str(line) for line in lines]
        electron_tolerance = 1.E-08
        rows = len(temperature[0])

        emissivity_map = [np.zeros(shape) for shape in [arr.shape for arr in temperature]]
        emissivity_map_dict = {line: emissivity_map for line in lines_str}

        with Manager() as manager:
            shared_dict = manager.dict(emissivity_map_dict)

            args = []
            for level in range(N_grid_level):
                ne[level] = np.array(ne[level])
                ne[level][ne[level] == 0] = electron_tolerance
                args.extend(
                    (row, level, temperature[level][row], ne[level][row], lines, lines_str)
                    for row in range(rows)
                )

            with Pool(processes=ncores) as pool:
                results = pool.map(self.compute_row_emissivity, args)

            # Aggregate results back into the shared dictionary
            for (row, level, _, _, _, _), result in zip(args, results):
                for line, value in result.items():
                    shared_dict[line][level][row] = value

            emissivity_map_dict = dict(shared_dict)

        return emissivity_map_dict

