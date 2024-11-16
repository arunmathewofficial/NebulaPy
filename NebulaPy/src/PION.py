import re
import numpy as np
from pypion.ReadData import ReadData
from pypion.SiloHeader_data import OpenData
from NebulaPy.tools import util as util
from NebulaPy.tools import constants as const
import astropy.units as unit

class pion():
    '''
    This class is not an alternative to Pypion; rather, it is a
    bundle of methods useful for creating synthetic emission
    maps from the Silo file.
    '''

    def __init__(self, silo_set, verbose):
        self.silo_set = silo_set
        self.verbose = verbose
        self.geometry_container = {}
        self.chemistry_container = {}


    def load_geometry(self, silo_instant, scale='cm'):
        '''
        This method will load geometry of the simulation from
        the given silo file.

        Parameters
        ----------
        silo_instant
        scale

        Returns
        -------

        '''
        # If verbose is enabled, print the chemistry code
        if self.verbose:
            print(f" ---------------------------")
            print(f" loading geometry:")

        # Open the data for the first silo instant silo
        header_data = OpenData(silo_instant)
        # Set the directory to '/header'
        header_data.db.SetDir('/header')
        # Retrieve what coordinate system is used
        coord_sys = header_data.db.GetVar("coord_sys")

        if coord_sys == 3:
            if self.verbose:
                print(f" geometry: {const.coordinate_system[coord_sys]}")
                self.spherical_grid(silo_instant, scale=scale)
        elif coord_sys == 2:
            if self.verbose:
                print(f" geometry: {const.coordinate_system[coord_sys]}")
                self.cylindrical_grid(silo_instant, scale=scale)
        elif coord_sys == 1:
            if self.verbose:
                print(f" geometry: {const.coordinate_system[coord_sys]}")
                util.nebula_exit_with_error(f"{const.coordinate_system[coord_sys]} not defined, todo list")


    ######################################################################################
    # spherical grid
    ######################################################################################
    def spherical_grid(self, silo_instant, scale='cm'):

        # Open the data for the first silo instant silo
        header_data = OpenData(silo_instant)
        # Set the directory to '/header'
        header_data.db.SetDir('/header')
        # Retrieve what coordinate system is used
        coord_sys = header_data.db.GetVar("coord_sys")
        if not coord_sys == 3:
            util.nebula_exit_with_error(f"geometry mismatch {const.coordinate_system[coord_sys]}")
        # Retrieve no of nested grid levels
        Nlevels = header_data.db.GetVar("grid_nlevels")
        # close the object
        header_data.close()

        # save the dynamics and chemistry_flag values in the chemistry_container dictionary
        self.geometry_container['coordinate_sys'] = const.coordinate_system[coord_sys]
        self.geometry_container['Nlevels'] = Nlevels
        if self.verbose:
            print(f" N grid levels: {Nlevels}")

        # read silo file
        data = ReadData(silo_instant)
        basic = data.get_1Darray('Density')
        mask = data.get_1Darray("NG_Mask")['data']
        mask = np.array(mask)
        # radial axis
        rmax = (basic['max_extents'] * unit.cm)
        rmin = (basic['min_extents'] * unit.cm)
        Ngrid = data.ngrid()
        # close the object
        data.close()

        # calculating radial points
        if self.verbose:
            print(' calculating radial points')
        radius = []
        for level in range(Nlevels):
            level_min = rmin[level].value
            level_max = rmax[level].value
            level_dr = (level_max[0] - level_min[0]) / Ngrid[0]
            r0 = level_min[0] + 0.5 * level_dr
            rn = level_max[0] - 0.5 * level_dr
            r = np.linspace(r0, rn, Ngrid[0])
            radius.append(r)  # append radius of each level

        if Nlevels > 1:
            # last element of the tracer array tracer[Nlevels - 1]
            fine_level = radius[Nlevels - 1]
            # Loop through the tracer array starting from the second-to-last element down to the first element,
            # goes from Nlevels - 2 (second-to-last element) to 0 (first element)
            for i in range(Nlevels - 2, -1, -1):
                # Use the mask array to selectively delete certain elements from tracer[i]. np.where(mask[i] == 0)
                # finds the indices in mask[i] where the value is 0. np.delete(tracer[i], np.where(mask[i] == 0))
                # removes the elements from tracer[i] at those indices.
                coarse_level = np.delete(radius[i], np.where(mask[i] == 0))
                # append the filtered array coarse_level to the result array to fine_level.
                fine_level = np.append(fine_level, coarse_level)
            radius = np.array(fine_level)

        # if the data is single level (uniform grid)
        if Nlevels == 1:
            radius = radius[0] * mask[0]

        # calculating shell volumes
        if self.verbose:
            print(' calculating shell volumes')
        # Calculating the core volume
        core = 4.0 * const.pi * radius[0] ** 3.0 / 3.0
        # Calculating the shell volumes
        shell_volumes = 4.0 * const.pi * (radius[1:] ** 3 - radius[:-1] ** 3) / 3.0
        # Insert the core volume at the beginning of the shell_volumes array
        shell_volumes = np.insert(shell_volumes, 0, core)

        self.geometry_container['radius'] = radius
        self.geometry_container['shell_volumes'] = shell_volumes

    ######################################################################################
    # cylindrical grid
    ######################################################################################
    def cylindrical_grid(self, silo_instant, scale='cm'):
        # Open the data for the first silo instant silo
        header_data = OpenData(silo_instant)
        # Set the directory to '/header'
        header_data.db.SetDir('/header')
        # Retrieve what coordinate system is used
        coord_sys = header_data.db.GetVar("coord_sys")
        if not coord_sys == 2:
            util.nebula_exit_with_error(f"geometry mismatch {const.coordinate_system[coord_sys]}")
        # Retrieve no of nested grid levels
        Nlevels = header_data.db.GetVar("grid_nlevels")
        # close the object
        header_data.close()

        # save the dynamics and chemistry_flag values in the chemistry_container dictionary
        self.geometry_container['coordinate_sys'] = const.coordinate_system[coord_sys]
        self.geometry_container['Nlevels'] = Nlevels
        if self.verbose:
            print(f" N grid levels: {Nlevels}")

        # Read the data from the current silo file
        dataio = ReadData(silo_instant)
        basic = dataio.get_2Darray('Density')  # Retrieve basic simulation data, such as density
        dataio.close()  # Close the data file

        if scale == 'cm':
            dims_max = (basic['max_extents'] * unit.cm)
            dims_min = (basic['min_extents'] * unit.cm)
            self.geometry_container['edges_min'] = dims_min
            self.geometry_container['edges_max'] = dims_max
        elif scale == 'pc':
            dims_max = (basic['max_extents'] * unit.cm).to(unit.pc)
            dims_min = (basic['min_extents'] * unit.cm).to(unit.pc)
            self.geometry_container['edges_min'] = dims_min
            self.geometry_container['edges_max'] = dims_max


    ######################################################################################
    # get simulation time
    ######################################################################################
    def get_simulation_time(self, silo_instant, time_unit='sec'):

        if self.geometry_container['coordinate_sys'] == 'spherical':
            # Read the data from the current silo file
            dataio = ReadData(silo_instant)
            basic = dataio.get_1Darray('Density')  # Retrieve basic simulation data, such as density
            dataio.close()  # Close the data file

            if time_unit == 'sec':
                return (basic['sim_time'] * unit.s)
            elif time_unit == 'kyr':
                return (basic['sim_time'] * unit.s).to(unit.kyr)

        elif self.geometry_container['coordinate_sys'] == 'cylindrical':
            # Read the data from the current silo file
            dataio = ReadData(silo_instant)
            basic = dataio.get_2Darray('Density')  # Retrieve basic simulation data, such as density
            dataio.close()  # Close the data file
            if time_unit == 'sec':
                return (basic['sim_time'] * unit.s)
            elif time_unit == 'kyr':
                return (basic['sim_time'] * unit.s).to(unit.kyr)

    ######################################################################################
    # get chemistry from the initial silo file
    ######################################################################################
    def load_chemistry(self):
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

        # Open the data for the first silo instant silo
        header_data = OpenData(self.silo_set[0])
        # Set the directory to '/header'
        header_data.db.SetDir('/header')
        #print(header_data.header_info())
        # Retrieve the value of "EP_chemistry" from the header data
        chemistry_flag = header_data.db.GetVar("EP_chemistry")
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
                util.nebula_exit_with_error(" PION is not running NEMO v1.0; NelubaPy functionality is limited.")
            else:
                # If verbose is enabled, print the chemistry code
                if self.verbose:
                    print(f" ---------------------------")
                    print(f" loading chemistry:")
                    print(f" chemistry module: NEMO v1.0")

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
                # elements in the tracer list
                tracer_elements = []
                # mass_fraction
                mass_fractions = {}
                # list of element wise tracer list
                elementWiseTracers = [[] for _ in range(len(const.nebula_elements))]
                # list of element names from the nebula_elements dictionary keys
                element_list = list(const.nebula_elements.keys())
                # save the number of tracers in the chemistry_container dictionary
                self.chemistry_container['Ntracers'] = Ntracers
                # If verbose is enabled, print the number of chemical tracers
                if self.verbose:
                    print(f" N chemical tracers: {Ntracers}")

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
                        tracer_elements.append(element)
                        # get the full element name from the nebula_elements dictionary
                        mass_fractions[element] = f'Tr{i:03}_' + chem_tracer
                        # if verbose is enabled, print the found element name
                        if self.verbose:
                            print(f" found {const.nebula_elements[element]}")
                        # get the index of the element in the element_list
                        element_index = element_list.index(element)
                        # append the tracer with the corresponding element to the mpv10tracers list
                        if 0 <= element_index < len(elementWiseTracers):
                            elementWiseTracers[element_index].append(f'Tr{i:03}_' + chem_tracer)

                    # check if the tracer is a corresponding ion
                    if re.sub(r'\d{1,2}\+', '', chem_tracer) in const.nebula_elements:
                        self.chemistry_container[chem_tracer] = f'Tr{i:03}_' + chem_tracer.replace('+', 'p')
                        # extract the element name
                        element = re.sub(r'\d{1,2}\+', '', chem_tracer)
                        # get the index of the element in the element_list
                        element_index = element_list.index(element)
                        # gppend the tracer with the corresponding ion to the mpv10tracers list
                        elementWiseTracers[element_index].append(f'Tr{i:03}_' + chem_tracer.replace('+', 'p'))

                # save mass fraction to chemistry_container dictionary
                #self.chemistry_container['mass_fractions'] = mass_fractions
                self.element_list = tracer_elements
                self.chemistry_container['mass_fractions'] = mass_fractions
                self.chemistry_container['tracer_elements'] = tracer_elements
                self.element_wise_tracer_list = elementWiseTracers
        header_data.close()

    ######################################################################################
    # get elements
    ######################################################################################
    def get_elements(self):
        return np.array(self.chemistry_container['tracer_elements'])

    ######################################################################################
    # get chemical tracers
    ######################################################################################
    def get_chemical_tracers(self):
        """
        Retrieve the list of chemical tracer strings for each tracer in the chemistry
        container dictionary, processed element by element. Each sublist starts with the
        mass fraction of the element followed by the tracers.

        Returns:
            list of lists: Each sublist contains the mass fraction followed by the values of
            the tracers for a specific element.
        """
        elements = self.get_elements()
        tracers = []

        for element in elements:
            # Retrieve tracers for the element
            element_tracers = [self.chemistry_container[f"{element}{q}+" if q > 0 else element]
                               for q in range(const.atomic_number[element])]

            tracers.append(element_tracers)

        return tracers

    ######################################################################################
    # get elemental mass fraction
    ######################################################################################
    def get_elemental_mass_frac(self, silo_instant):

        elements = self.get_elements()
        elemental_mass_fraction = []
        for element in elements:
            # Retrieve mass fraction
            element_tracer = self.chemistry_container['mass_fractions'][element]
            elemental_mass_fraction.append(self.get_parameter(element_tracer, silo_instant))

        return np.array(elemental_mass_fraction)

    ######################################################################################
    # get tracer values
    ######################################################################################
    def get_tracer_values(self, silo_instant):
        """
        Retrieves the chemical tracer values for the given time instant from the
        simulation silo data.
        Parameters:
        ----------
        silo_instant : silo file(s)

        Returns:
        -------
        tracer_values : list of lists
            A 2D list containing the tracer values for each ion in the tracers array.
        """

        # Retrieve the 2D array of chemical tracers.
        tracers = self.get_chemical_tracers()

        # Initialize tracer_values using list comprehension for better efficiency.
        tracer_values = np.array([
            [self.get_parameter(ion, silo_instant) for ion in element_row]
            for element_row in tracers
        ], dtype=object)

        return tracer_values

    ######################################################################################
    # get parameter //todo: this is not clear
    ######################################################################################
    def get_parameter(self, parameter, silo_instant):
        '''
        Method will return the parameter value for a spherical nested grid

        Parameters
        ----------
        parameter physical
        silo_instant

        Returns
        -------
        physical parameter value for a spherical nested grid
        '''

        # 1 dimensional (spherical) ######################################################
        if self.geometry_container['coordinate_sys'] == 'spherical':
            # get nested grid level
            Nlevels = self.geometry_container['Nlevels']

            # pypion ReadDate object
            data = ReadData(silo_instant)
            # get parameter values
            parameter = data.get_1Darray(parameter)['data']
            # get mask
            mask = data.get_1Darray("NG_Mask")['data']
            data.close()
            # if the data is single level (uniform grid)
            if Nlevels == 1:
                return parameter[0] * mask[0]

            # last element of the parameter array parameter[Nlevels - 1]
            fine_level = parameter[Nlevels - 1]
            # Loop through the parameter array starting from the second-to-last element down to the first element,
            # goes from Nlevels - 2 (second-to-last element) to 0 (first element)
            for i in range(Nlevels - 2, -1, -1):
                # Use the mask array to selectively delete certain elements from parameter[i]. np.where(mask[i] == 0)
                # finds the indices in mask[i] where the value is 0. np.delete(parameter[i], np.where(mask[i] == 0))
                # removes the elements from parameter[i] at those indices.
                coarse_level = np.delete(parameter[i], np.where(mask[i] == 0))
                # append the filtered array coarse_level to the result array to fine_level.
                fine_level = np.append(fine_level, coarse_level)
            return np.array(fine_level)
        # end of 1 dimensional ***********************************************************

        # 2 dimensional (cylindrical) ####################################################
        if self.geometry_container['coordinate_sys'] == 'cylindrical':

            # get nested grid level
            Nlevels = self.geometry_container['Nlevels']
            # pypion ReadDate object
            data = ReadData(silo_instant)
            # get parameter values
            parameter = data.get_2Darray(parameter)['data']
            # get mask
            mask = data.get_2Darray("NG_Mask")['data']
            data.close()
            return parameter
        # end of 2 dimensional ***********************************************************
        # 3 dimensional (cartesian) ######################################################
        # end of 3 dimensional ***********************************************************


    ######################################################################################
    # get ion mass fraction values
    ######################################################################################
    def get_ion_values(self, ion, silo_instant):
        '''
        This methods will return the ion mass fraction value set

        Parameters
        ----------
        ion name
        silo_instant

        Returns
        -------
        ion mass fraction
        '''

        ion_tracer = None
        if ion not in self.chemistry_container:
            util.nebula_exit_with_error(f"ion {ion} not in {self.chemistry_container['chemistry_code']} chemistry")
        else:
            ion_tracer = self.chemistry_container[ion]

        return self.get_parameter(ion_tracer, silo_instant)


    ######################################################################################
    # get electron number density
    ######################################################################################
    def get_ne(self, silo_instant):
        '''

        Parameters
        ----------
        silo_instant

        Returns
        -------
        electron number density in each cell for a specific silo file
        '''
        # 1 dimensional (spherical) ######################################################
        if self.geometry_container['coordinate_sys'] == 'spherical':

            density = self.get_parameter("Density", silo_instant)
            ne = np.zeros(len(density))
            for e, element in enumerate(self.element_wise_tracer_list):
                if not element:  # Check if the element is empty
                    continue
                massfrac_sum = np.zeros(len(density))
                element_name = self.element_list[e]
                atomic_number = len(element) - 1
                top_ion = self.get_parameter(element[0], silo_instant)

                for i, ion in enumerate(element[1:], start=1):
                    charge = i - 1
                    ion_density = self.get_parameter(ion, silo_instant)
                    top_ion -= ion_density
                    massfrac_sum += charge * ion_density

                massfrac_sum += atomic_number * np.maximum(top_ion, 0.0)

                ne += massfrac_sum / const.mass[element_name]

            return ne * density

        # 2 dimensional (cylindrical) ####################################################
        if self.geometry_container['coordinate_sys'] == 'cylindrical':

            # Get the number of nested grid levels in the geometry container.
            Nlevels = self.geometry_container['Nlevels']

            # Determine the number of grid levels for electron density calculation.
            if Nlevels == 1:
                print(" calculating electron number density for a single grid level")
            elif Nlevels > 1:
                print(" calculating electron number density for each grid levels")

            # Retrieve the density data from the input file at the current simulation instant.
            density = self.get_parameter("Density", silo_instant)

            # Identify the shape of each density array to ensure compatibility with other parameters.
            shape_list = [arr.shape for arr in density]

            # Initialize arrays for electron number density (ne) and mass fraction sum,
            # with zeroes matching the shape of the density data.
            ne = [np.zeros(shape) for shape in shape_list]
            massfrac_sum = [np.zeros(shape) for shape in shape_list]

            # Loop through each element in the tracer list to calculate contributions from individual ions.
            for e, element in enumerate(self.element_wise_tracer_list):

                # If the current element has no associated tracers, skip to the next element.
                if not element:
                    continue

                # Get the element name and compute its atomic number (total ions minus one).
                element_name = self.element_list[e]
                atomic_number = len(element) - 1

                # Get the top ion density data (for the element's highest ionization state).
                top_ion = self.get_parameter(element[0], silo_instant)

                # For each subsequent ion (starting from the next lowest ionization state), calculate contributions:
                for i, ion in enumerate(element[1:], start=1):
                    charge = i - 1  # Charge is one less than ionization state index.
                    ion_density = self.get_parameter(ion, silo_instant)

                    # Update the top ion and mass fraction sum for each grid level.
                    for level in range(Nlevels):
                        top_ion[level] -= ion_density[
                            level]  # Adjust top ion density to reflect current ion's contribution.
                        massfrac_sum[level] += charge * ion_density[level]  # Update mass fraction sum by ion charge.

                # Finalize mass fraction sum and update electron density for each grid level.
                for level in range(Nlevels):
                    # Add contribution of the top ion with atomic number, ensuring non-negative values.
                    massfrac_sum[level] += atomic_number * np.maximum(top_ion[level], 0.0)
                    # Calculate electron number density using mass fraction sum and atomic mass.
                    ne[level] += massfrac_sum[level] / const.mass[element_name]

            # Scale the electron number density by the density for each grid level.
            ne = [density[level] * ne[level] for level in range(Nlevels)]

            # Return the calculated electron number density array for each grid level.
            return ne


    ######################################################################################
    # differential emission measure
    ######################################################################################
    import numpy as np

    def DEM(self, dem_indices, ne, shellvolume):
        """
        Calculate the differential emission measure (DEM) across temperature bins.

        Parameters:
        ----------
        dem_indices : list of numpy.ndarray
            List where each element is an array of indices corresponding to a temperature bin.
        ne : numpy.ndarray
            Array of electron densities corresponding to the temperature values.
        shellvolume : numpy.ndarray
            Array of shell volumes corresponding to the temperature values.

        Returns:
        -------
        DEM : numpy.ndarray
            Array of differential emission measure values for each temperature bin.
        """

        # Calculate ne * ne * dV for all elements
        volume_ne_square = ne * ne * shellvolume

        # Initialize an array for DEM with the same length as dem_indices
        DEM = np.zeros(len(dem_indices))

        # Use list comprehension and NumPy's array indexing to sum the values for each bin
        for i, indices in enumerate(dem_indices):
            if indices.size > 0:
                DEM[i] = np.sum(volume_ne_square[indices])

        # Remove zero values from DEM
        DEM = DEM[DEM > 0]
        return DEM

    ######################################################################################
    # generate differential emission measure indices
    ######################################################################################
    def generate_dem_indices(self, temperature, Tmin, Tmax, Nbins):
        # Calculate the logarithmic width of each bin
        bin_width = (np.log10(Tmax) - np.log10(Tmin)) / Nbins
        # Half the width of a bin for adjusting bin edges
        half_bin_width = bin_width / 2

        # Generate the logarithmically spaced temperature bin edges.
        # This will create Nbins+1 edges to define the boundaries of Nbins.
        temperature_edges = np.linspace(np.log10(Tmin), np.log10(Tmax), Nbins + 1)

        # Create temperature bins by pairing adjacent edges.
        # Each bin is represented as [bin_min, bin_max].
        temperature_bins = [[temperature_edges[i], temperature_edges[i + 1]] for i in range(Nbins)]

        # Calculate the midpoints of each bin for potential further use.
        # These midpoints are the average of the logarithmic bin edges.
        Tb = np.array([(bin[0] + bin[1]) / 2 for bin in temperature_bins])

        # Initialize an empty list to store indices of temperatures within each bin.
        dem_indices = []

        # Loop through each bin to identify temperature values that fall within the bin's range.
        for i in range(Nbins):
            # Find the indices of temperature values that fall within the current bin.
            # The condition checks if the logarithm of the temperature is within the bin range,
            # slightly adjusted by half_bin_width to ensure proper capturing of boundary values.
            indices = np.where((np.log10(temperature) >= temperature_bins[i][0] - half_bin_width) &
                               (np.log10(temperature) < temperature_bins[i][1] + half_bin_width))[0]

            # Append the array of indices to the dem_indices list.
            dem_indices.append(indices)

        # Filter Tb values for which dem_indices[i] is not empty
        Tb = [Tb[i] for i in range(Nbins) if len(dem_indices[i]) > 0]

        # Return the list of indices for further processing or analysis.
        return {'indices': dem_indices, 'Tb': Tb}








