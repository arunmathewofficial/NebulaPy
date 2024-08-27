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

    ######################################################################################
    # spherical grid
    ######################################################################################
    def spherical_grid(self, silo_instant):

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
            print(f" ---------------------------")
            print(f" geometry: {const.coordinate_system[coord_sys]}")
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
                util.nebula_exit_with_error(" PION is not running MPv10; NelubaPy functionality is limited.")
            else:
                # If verbose is enabled, print the chemistry code
                if self.verbose:
                    print(f" ---------------------------")
                    print(f" loading chemistry ...")
                    print(f" chemistry module: {chemistry_code}")

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
        # end of 2 dimensional ***********************************************************
        # 3 dimensional (cylindrical) ####################################################
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

    ######################################################################################
    # differential emission measure
    ######################################################################################
    def differential_emission_measure(self, temperature, ne, shellvolume, Tmin, Tmax, Nbins):
        """
        Calculate the differential emission measure (DEM) across temperature bins.

        The temperature range from Tmin to Tmax is divided into Nbins logarithmically spaced bins.
        The DEM is computed for each bin based on the provided electron density, shell volume,
        and temperature values. The midpoint of each temperature bin is also calculated.

        Parameters:
        ----------
        temperature : numpy.ndarray
            Array of temperature values.
        ne : numpy.ndarray
            Array of electron densities corresponding to the temperature values.
        shellvolume : numpy.ndarray
            Array of shell volumes corresponding to the temperature values.
        Tmin : float
            The minimum temperature for binning.
        Tmax : float
            The maximum temperature for binning.
        Nbins : int
            The number of logarithmically spaced bins.

        Returns:
        -------
        dem : numpy.ndarray
            Array of differential emission measure values for each temperature bin.
        temperature_midpoints : numpy.ndarray
            Array of the midpoints of logarithmically spaced temperature bins.
        temperature_bins : numpy.ndarray
            Array of logarithmically spaced temperature bins.
        """

        # Calculate the logarithmic step size.
        dex = np.log10(Tmax / Tmin) / Nbins

        # Generate the temperature bins logarithmically spaced.
        temperature_bins = np.logspace(np.log10(Tmin), np.log10(Tmax), Nbins + 1)

        # Calculate the midpoints of each bin.
        temperature_midpoints = np.sqrt(temperature_bins[:-1] * temperature_bins[1:])

        # Initialize the DEM array with zeros.
        dem = np.zeros(Nbins)

        # Calculate the differential emission measure for each bin.
        for i in range(Nbins):
            # Identify the indices of temperature values that fall within the current bin.
            indices = np.where((temperature >= temperature_bins[i]) & (temperature < temperature_bins[i + 1]))[0]

            # Sum the DEM for the current bin.
            dem[i] = np.sum(ne[indices] ** 2 * shellvolume[indices])

        return {'DEM': dem, 'Tb': temperature_midpoints}











