import os
os.environ['XUVTOP']
import ChiantiPy
import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.io as chio
import ChiantiPy.tools.data as chdata
from ChiantiPy.base import specTrails
import numpy as np
import ChiantiPy.tools.util as chianti_util
from NebulaPy.tools import constants as const
from NebulaPy.tools import util as util
from ChiantiPy.core.Continuum import continuum
import ChiantiPy.tools.io as io
from scipy.interpolate import splev, splrep

# Limit the number of threads used by OpenMP and Intel MKL to 1.
# This prevents oversubscription of CPU cores when running multiprocessing tasks
# and ensures efficient resource utilization in parallel workflows.
os.environ["OMP_NUM_THREADS"] = "1"  # Restrict OpenMP to 1 thread
os.environ["MKL_NUM_THREADS"] = "1"  # Restrict Intel MKL to 1 thread


class chianti:
    """
    The class for calculating emission line spectrum.

    Parameters
    ----------
    Keyword arguments
    -----------------
    Examples
    --------
    >>> temperature = 1.e+9
    >>> ne = 1.e+4
    >>> ion = nebula.chianti('o_4', temperature, ne)
    >>> print(ion.emissivity())
    Notes
    -----
    References
    ----------
    """

    ######################################################################################
    #
    ######################################################################################
    def __init__(self, temperature, ne, chianti_ion=None, pion_ion=None, pion_elements=None, verbose=False):

        self.temperature = temperature
        self.ne = ne
        self.verbose = verbose

        # Count the number of arguments that are not None
        non_none_count = sum(arg is not None for arg in [chianti_ion, pion_ion, pion_elements])

        # If more than one argument is not None, raise a ValueError
        if non_none_count > 1:
            util.nebula_exit_with_error("invalid arguments: set only one of 'chianti_ion', 'pion_ion', or 'pion_elements")

        if pion_ion is not None:
            self.chianti_ion_name = self.get_chianti_symbol(pion_ion, make=False)
            self.chianti_ion = ch.ion(self.chianti_ion_name, temperature=self.temperature, eDensity=self.ne,
                                 pDensity='default', radTemperature=None, rStar=None, abundance=None,
                                 setup=True, em=None, verbose=self.verbose)
            self.get_ion_attributes()

        if chianti_ion is not None:
            self.chianti_ion = ch.ion(chianti_ion, temperature=self.temperature, eDensity=self.ne,
                                 pDensity='default', radTemperature=None, rStar=None, abundance=None,
                                 setup=True, em=None, verbose=self.verbose)
            self.chianti_ion_name = chianti_ion

        if pion_elements is not None:
            chianti_element_list = []
            for element in pion_elements:
                element_symbol = self.get_chianti_symbol(element, make=True)
                chianti_element_list.append(element_symbol)
            self.chianti_element_list = chianti_element_list
            self.get_elements_attributes()

    ######################################################################################
    # generate species chianti symbol
    ######################################################################################
    def get_chianti_symbol(self, species, make=False):
        '''
        Converts a PION species symbol to a corresponding CHIANTI species symbol.

        This function takes a species symbol used in PION (a computational tool) and
        converts it into the format required by CHIANTI, a database for atomic data.
        The function can generate either the elemental symbol or the ionized symbol
        depending on the 'make' parameter.

        :param species: str, the PION species symbol, which can represent both neutral
                        and ionized states.
        :param make: bool, if True, returns the elemental symbol (e.g., 'h' for hydrogen).
                     if False, returns the CHIANTI ion symbol (e.g., 'h_2' for H+).
        :return: str, the corresponding CHIANTI symbol for the species.
        '''

        # Convert the input species symbol to lowercase and remove any '+' characters
        # (denoting ionization)
        species = species.lower().replace('+', '')

        # Extract alphabetic characters to identify the element symbol (e.g., 'h' from 'h1' or 'h+')
        element = ''.join(filter(str.isalpha, species))

        if make:
            # If 'make' is True, return only the element symbol (e.g., 'h')
            return element
        else:
            # Extract numeric characters to determine the ionization level (e.g., '1' from 'h1')
            ion_level = ''.join(filter(str.isdigit, species))

            # If no numeric characters are found, set the CHIANTI ionization level to 1
            chianti_level = int(ion_level) + 1 if ion_level else 1

            # Return the element symbol followed by the ionization level, separated
            # by an underscore (e.g., 'h_2')
            return f"{element}_{chianti_level}"

    ######################################################################################
    # get attributes of all elements in chianti element list
    ######################################################################################
    def get_elements_attributes(self):
        '''
        Generates and appends species attributes for each element in the
        `chianti_element_list` to a dictionary called `species_attributes`.

        This function first retrieves abundance data from the `chdata.Abundance`
        dictionary using the specified abundancoe name ('unity'). It then initializes
        an instance of the `specTrails` class to handle species-related data, setting
        its abundance and temperature properties.

        The function calls the `ionGate` method on the `species` object to process
        the elements in `chianti_element_list`, generating species data based on
        the specified parameters.

        The species data is then looped over, and for each element key in the sorted
        `species.Todo` dictionary:
        - The key is converted using `util.convertName` and stored in `species_attributes`.
        - The relevant data is stored in the `species_attributes` dictionary under
          each element's key.
        - Unnecessary entries like 'filename' and 'experimental' are removed
          from the dictionary before final storage.

        :return: None
        '''
        AbundanceName = 'unity'
        abundAll = chdata.Abundance[AbundanceName]['abundance']

        species = specTrails()  # Create an instance of specTrails to manage species data
        species.AbundAll = abundAll
        species.Temperature = self.temperature  # Set the temperature for the species

        species.ionGate(
            elementList=self.chianti_element_list,
            minAbund=None, doLines=True,
            doContinuum=True, doWvlTest=0,
            doIoneqTest=0, verbose=False
        )

        self.species_attributes_container = {}

        # Loop through the sorted keys in the dictionary of species
        if self.verbose:
            print(f" retrieving species attributes")

        count = 0
        for akey in sorted(species.Todo.keys()):
            self.species_attributes_container[akey] = chianti_util.convertName(akey)  # Convert the key and store it
            # If verbose mode is enabled, print the spectroscopic name
            if self.verbose:
                # Print a comma-separated list of names with up to 10 items per line
                print(f" {self.species_attributes_container[akey]['spectroscopic']}", end='')
                count += 1
                # Print a newline after every 10 items
                if count % 10 == 0:
                    print()  # Move to the next line
                else:
                    print(", ", end='')  # Continue on the same line

            self.species_attributes_container[akey]['keys'] = species.Todo[akey]  # Store relevant data
            # Remove unnecessary data from the dictionary
            del self.species_attributes_container[akey]['filename']
            del self.species_attributes_container[akey]['experimental']


        # Finalize the species attributes dictionary
        # At this point, `self.species_attributes` contains all the relevant
        # attributes for the species in `chianti_element_list`

    ######################################################################################
    # get ion attributes
    ######################################################################################
    def get_ion_attributes(self):

        AbundanceName = 'unity'
        abundAll = chdata.Abundance[AbundanceName]['abundance']

        species = specTrails()  # Create an instance of specTrails to manage species data
        species.AbundAll = abundAll
        species.Temperature = self.temperature  # Set the temperature for the species

        ion_list = [self.chianti_ion_name]
        species.ionGate(ionList=ion_list,
            minAbund=None, doLines=True,
            doContinuum=True, doWvlTest=0,
            doIoneqTest=0, verbose=False
        )

        self.species_attributes_container = {}

        # Loop through the sorted keys in the dictionary of species
        if self.verbose:
            print(f" retrieving species attributes")

        count = 0
        for akey in sorted(species.Todo.keys()):
            self.species_attributes_container[akey] = chianti_util.convertName(akey)  # Convert the key and store it
            # If verbose mode is enabled, print the spectroscopic name
            if self.verbose:
                # Print a comma-separated list of names with up to 10 items per line
                print(f" {self.species_attributes_container[akey]['spectroscopic']}", end='')
                count += 1
                # Print a newline after every 10 items
                if count % 10 == 0:
                    print()  # Move to the next line
                else:
                    print(", ", end='')  # Continue on the same line

            self.species_attributes_container[akey]['keys'] = species.Todo[akey]  # Store relevant data
            # Remove unnecessary data from the dictionary
            del self.species_attributes_container[akey]['filename']
            del self.species_attributes_container[akey]['experimental']

    ######################################################################################
    # get all lines of the ion
    ######################################################################################
    def get_allLines(self):
        """
                   Retrieve all spectral lines associated with a specified ion
                   :return: wave-length array
                   """
        if self.verbose:
            print(' retrieving all spectral lines of ', self.chianti_ion.Spectroscopic)
        wvl = np.asarray(self.chianti_ion.Wgfa['wvl'], np.float64)
        wvl = np.abs(wvl)
        return wvl

    ######################################################################################
    # get all lines and transitions of the ion
    ######################################################################################
    def get_allLinesTransitions(self):
        """
        Retrieve all spectral lines associated with a specified ion
        :return: wave-length array
        """
        if self.verbose:
            print(' retrieving spectral lines and transitions of ', self.chianti_ion.Spectroscopic)

        Ref = self.chianti_ion.Elvlc['ref']
        A_value = np.asarray(self.chianti_ion.Wgfa['avalue'], np.float64)
        Pretty1 = self.chianti_ion.Wgfa['pretty1']
        Pretty2 = self.chianti_ion.Wgfa['pretty2']
        wvl = np.asarray(self.chianti_ion.Wgfa['wvl'], np.float64)

        wvl = np.abs(wvl)
        A_value = np.abs(A_value)
        return {'Reference': Ref, 'wvl': wvl, 'Avalue': A_value, 'From': Pretty1, 'To': Pretty2}

    ######################################################################################
    # Get line emissivity, this is an internal method
    ######################################################################################
    def get_line_emissivity(self, allLines=True):
        """
               Retrieve the emissivity values for all spectral lines associated
               with a specified ion.

               :return: Dict Emiss. Emiss has several quantities, namely, ion,
               # wvl(angstrom), emissivity (ergs s^-1 str^-1), pretty1, pretty2.
               """
        if 'line' not in self.species_attributes_container[self.chianti_ion_name]['keys']:
            util.nebula_warning(f'no line emission associate with {self.chianti_ion.Spectroscopic}')
            return None

        else:
            if self.verbose:
                print(' retrieving emissivity values for all spectral lines '
                      'of', self.chianti_ion.Spectroscopic)
            self.chianti_ion.emiss(allLines=allLines)
            emissivity = self.chianti_ion.Emiss
            return emissivity



    ######################################################################################
    # get line emissivity for a list of lines, this is an internal method
    ######################################################################################
    def get_line_emissivity_for_list(self, line_list):
        """
        Retrieves the emissivity values for a given list of spectral lines.

        Parameters:
        ----------
        line_list : list
            A list of spectral line identifiers for which emissivity values need to be retrieved.

        Returns:
        -------
        line_emissivity : dict
            A dictionary where keys are formatted line identifiers (e.g., "Fe XIV 530.3")
            and values are the corresponding emissivity values.
        """

        # Retrieve the full list of available spectral lines from the current object.
        all_lines = np.array(self.get_allLines())  # Get all available lines as a NumPy array.
        all_lines = all_lines[all_lines != 0]  # Remove any zero entries (which may indicate missing or invalid lines).

        # Print a message if verbose mode is enabled.
        if self.verbose:
            print(" retrieving line index for the given line(s)")

        # Check if every requested line exists in the available list of lines.
        missing_lines = [line for line in line_list if line not in all_lines]

        # If any requested lines are missing, terminate execution with an error message.
        if missing_lines:
            util.nebula_exit_with_error(f" following line(s) are not found: {missing_lines}")

        # Retrieve the indices of the requested lines within the all_lines array.
        # np.where(all_lines == line) returns an array of indices where the condition is met.
        # We take the first occurrence with [0][0] since it's assumed that each line appears only once.
        line_indices = [np.where(all_lines == line)[0][0] for line in line_list if line in all_lines]

        # Retrieve the full emissivity array from the method get_line_emissivity.
        # The `allLines=False` argument ensures that emissivity is returned only for relevant lines.
        emissivity = self.get_line_emissivity(allLines=False)['emiss']

        # Initialize a dictionary to store emissivity values for the requested lines.
        line_emissivity = {}

        # Iterate through the requested lines and their corresponding indices.
        for i, index in enumerate(line_indices):
            # Construct a human-readable identifier for the line.
            # The spectroscopic notation (e.g., "Fe XIV") is combined with the wavelength (or another identifier).
            line_str = self.chianti_ion.Spectroscopic + " " + str(line_list[i])

            # Retrieve the emissivity value corresponding to the current line.
            specific_line_emissivity = emissivity[index]

            # Store the retrieved emissivity value in the dictionary.
            line_emissivity[line_str] = specific_line_emissivity

        # Return the dictionary containing emissivity values for the requested lines.
        return line_emissivity

    ######################################################################################
    # get line spectrum
    ######################################################################################
    def get_line_spectrum(self, wavelength, species_density, shell_volume, allLines=True,
                          filtername=None, filterfactor=None):
        """
        Calculates the intensities for spectral lines of a specified ion, considering elemental
        abundance, ionization fraction, and emission measure.

        The method convolves the intensity results to simulate an observed spectrum. By default,
        it uses a Gaussian filter with a resolving power of 1000 to match the units of the continuum
        and line spectrum.

        Note:
        Emissivity has the unit \( \text{ergs} \, \text{s}^{-1} \, \text{str}^{-1} \).
        Intensity is given by:
        \[
        \text{Intensity} = \text{Ab} \times \text{ion\_frac} \times \text{emissivity} \times \frac{\text{em}}{\text{ne}}
        \]
        where the emission measure is given by:
        \[
        \text{em} = \int N_e \, N_H \, d\ell
        \]
        Intensity has the units \( \text{ergs} \, \text{cm}^{-2} \, \text{s}^{-1} \, \text{str}^{-1} \).

        Parameters
        ----------
        wavelength : array-like
            Array of wavelength values.
        Ab : float
            Elemental abundance.
        ion_frac : float
            Ionization fraction.
        em : array-like
            Emission measure values.
        allLines : bool, optional
            Whether to include all spectral lines (default is True).
        select_filter : str, optional
            Filter type to use for convolution, default is 'default'.
        factor : float, optional
            Factor for filter resolution, default is 1000.

        Returns
        -------
        line_spectrum : ndarray
            Array of convolved line intensities across the wavelength range.
        """
        if self.verbose:
            print(f" retrieving emissivity values for all spectral lines of {self.chianti_ion.Spectroscopic}")

        # Get emissivity of the specified ion (units: ergs s^-1 str^-1)
        self.chianti_ion.emiss(allLines)
        emissivity = self.chianti_ion.Emiss['emiss']

        # Number of temperature points and spectral lines
        N_temp = len(self.temperature)
        emission_measure = np.ones(N_temp)
        lines = self.chianti_ion.Emiss['wvl']
        N_lines = len(lines)

        # Initialize the intensity array
        intensity = np.zeros((N_temp, N_lines), dtype=np.float64)

        # Calculate intensity for each temperature
        if self.verbose:
            print(f" calculating line intensity for {self.chianti_ion.Spectroscopic}")
        for temp_idx in range(N_temp):
            intensity[temp_idx] = emissivity[:, temp_idx] * emission_measure[temp_idx] / self.ne[temp_idx]
        if self.verbose:
            print(f" {self.chianti_ion.Spectroscopic} line calculation completed")
        # Define the wavelength range and number of wavelength points
        wvl_range = [wavelength[0], wavelength[-1]]
        N_wvl = len(wavelength)

        # Select filter and factor
        filter = (chfilters.gaussianR, 1000.)
        useFilter = filter[0]
        useFactor = filter[1]

        # Initialize the line spectrum array
        line_spectrum = np.zeros((N_temp, N_wvl), dtype=np.float64)

        # Get indices of lines within the wavelength range
        selected_idx = chianti_util.between(lines, wvl_range)

        if len(selected_idx) == 0:
            if self.verbose:
                print(f' no lines found for {self.chianti_ion.Spectroscopic} in the wavelength '
                      f'range {wvl_range[0]:.2e} - {wvl_range[1]:.2e} Angstrom')
                print(' skipping ...')
        else:
            # Convolve the intensities with the filter for each temperature
            for temp_idx in range(N_temp):
                for wvl_idx in selected_idx:
                    line = lines[wvl_idx]
                    line_spectrum[temp_idx] += useFilter(wavelength, line, factor=useFactor) \
                                               * intensity[temp_idx, wvl_idx]

        # differential volume emission measure
        DVEM = self.ne * species_density * shell_volume
        # multiplying
        line_spectrum = line_spectrum * DVEM[:, np.newaxis]

        return line_spectrum


    ######################################################################################
    # get free-free emission
    ######################################################################################
    def get_bremsstrahlung_emission(self, wavelength, species_density, shell_volume):
        """
        Calculates the free-free emission (bremsstrahlung) for a single ion using the following formula:
        .. math::
           \\frac{dW}{dtdVd\lambda} = \\frac{c}{3m_e}\\left(\\frac{\\alpha h}{\pi}\\right)^3
           \\left(\\frac{2\pi}{3m_e k_B}\\right)^{1/2}\\frac{Z^2}{\lambda^2 T^{1/2}}
           \exp\\left(-\\frac{hc}{\lambda k_B T}\\right) \\bar{g}_{ff},

        where :math:`\nu = c/\lambda`, :math:`\\alpha` is the fine structure constant, :math:`Z` is the nuclear charge,
        and :math:`\\bar{g}_{ff}` is the velocity-averaged Gaunt factor.

        The free-free emission is calculated in units of
        :math:`\mathrm{erg}\ \mathrm{cm}^3\ \mathrm{s}^{-1}\ \mathrm{\mathring{A}}^{-1}\ \mathrm{str}^{-1}`.
        If the emission measure is provided, the result will be multiplied by
        :math:`\mathrm{cm}^{-5}` (for line-of-sight emission measure) or
        :math:`\mathrm{cm}^{-3}` (for volumetric emission measure).

        Parameters:
        ----------
        wavelength : array-like
            The wavelength(s) at which to calculate the emission, in angstroms.
        elemental_abundance : float
            The abundance of the element in the plasma.
        ion_fraction : float
            The fraction of the element in the ionization state of interest.
        emission_measure : float
            The emission measure, which may be line-of-sight or volumetric.

        Returns:
        -------
        free_free_emission : numpy.ndarray
            The calculated free-free emission for the given wavelength(s).
        """

        # Calculate the ion's nuclear charge (Z)
        Zion = self.chianti_ion.Ion - 1

        # Create a continuum object for calculating Gaunt factors
        continuum_spectrum = continuum(
            self.chianti_ion_name,
            temperature=self.temperature,
            abundance=None,
            em=None,
            verbose=self.verbose)

        # If verbose, print the ion's spectroscopic label
        if self.verbose:
            print(f' calculating bremsstrahlung emission for {self.chianti_ion.Spectroscopic}')

        # Ensure wavelength is treated as an array
        wavelength = np.atleast_1d(wavelength)

        # Define the numerical prefactor for the emission formula
        prefactor = ((const.c * 1e8) / (3. * const.emass) *
                     (const.alpha * const.h / const.pi) ** 3 *
                     np.sqrt(2. * const.pi / (3. * const.emass * const.kB)))

        # Include temperature dependence in the prefactor
        prefactor *= Zion ** 2 / np.sqrt(self.temperature)

        # Apply the elemental abundance and ion fraction to the prefactor
        #prefactor *= elemental_abundance * ion_fraction

        # Include the emission measure in the prefactor
        #prefactor *= emission_measure

        # Calculate the exponential factor based on temperature and wavelength
        exp_factor = np.exp(-const.planck * (1.e8 * const.light) / const.boltzmann /
                            np.outer(self.temperature, wavelength)) / (wavelength ** 2)

        # Calculate the Gaunt factor using the continuum spectrum object
        gf_itoh = continuum_spectrum.itoh_gaunt_factor(wavelength)
        gf_sutherland = continuum_spectrum.sutherland_gaunt_factor(wavelength)
        gf = np.where(np.isnan(gf_itoh), gf_sutherland, gf_itoh)

        # Optionally, apply an energy factor if flux is in photons (commented out by default)
        energy_factor = 1.0
        # if chdata.Defaults['flux'] == 'photon':
        #     energy_factor = const.planck * (1.e8 * const.light) / wavelength

        # Calculate the final free-free emission and ensure the result is properly shaped
        bremsstrahlung_emission = (prefactor[:, np.newaxis] * exp_factor * gf / energy_factor).squeeze()

        # If verbose, indicate completion
        if self.verbose:
            print(f' {self.chianti_ion.Spectroscopic} bremsstrahlung emission calculation completed')

        # differential volume emission measure
        DVEM = self.ne * species_density * shell_volume
        # multiplying
        bremsstrahlung_emission = bremsstrahlung_emission * DVEM[:, np.newaxis]

        return bremsstrahlung_emission

    ######################################################################################
    # get free-bound emission
    ######################################################################################
    def get_freebound_emission(self, wavelength, species_density, shell_volume, verner=True):
        """
        Calculates the free-bound (radiative recombination) continuum emissivity of an ion.

        Parameters
        ----------
        wavelength : numpy.ndarray
            Array of wavelengths in Angstroms where the emissivity is computed.
        elemental_abundance : float
            Abundance of the element in the astrophysical environment.
        ion_fraction : float
            Fraction of the ionized species in the environment.
        emission_measure : numpy.ndarray
            Array of emission measures at each temperature.
        verner : bool, optional
            If True, use the Verner-Yakovlev photoionization cross-sections. Default is True.

        Returns
        -------
        numpy.ndarray
            Array of emissivity in units of ergs cm^(-2) s^(-1) sr^(-1) Angstrom^(-1) for the given ion.

        Notes
        -----
        - Uses the Gaunt factors of CHIANTI V10 for recombination to the excited levels.
        - Uses the photoionization cross-sections to develop the free-bound cross-section.
        - Revised to calculate the free-bound cross-section and Maxwell energy distribution.
        """

        # Create a continuum object for calculating Gaunt factors
        continuum_spectrum = continuum(
            self.chianti_ion_name,
            temperature=self.temperature,
            abundance=None,
            em=None,
            verbose=self.verbose
        )

        Nwvl = wavelength.size
        Ntemp = self.temperature.size
        emission_measure = np.ones(Ntemp)

        # Generate a sequence of indices for temperature
        goodT = np.arange(Ntemp)

        if self.verbose:
            print(f' calculating free-bound emission for {self.chianti_ion.Spectroscopic}')

        # Load free-bound level data (Fblvl)
        if hasattr(continuum_spectrum, 'Fblvl'):
            fblvl = continuum_spectrum.Fblvl
            if 'errorMessage' in fblvl.keys():
                continuum_spectrum.FreeBound = fblvl
                print(' ' + fblvl['errorMessage'])
                return np.zeros((Ntemp, Nwvl), dtype=np.float64)
        elif continuum_spectrum.Z == continuum_spectrum.Stage - 1:
            # Fully ionized stage, assign default Fblvl
            continuum_spectrum.Fblvl = {'mult': [1., 1.]}
            fblvl = continuum_spectrum.Fblvl
        else:
            fblvlname = continuum_spectrum.nameDict['filename'] + '.fblvl'
            if os.path.isfile(fblvlname):
                continuum_spectrum.Fblvl = io.fblvlRead(self.chianti_ion_name)
                fblvl = continuum_spectrum.Fblvl
            else:
                if self.verbose:
                    print(f' no Fblvl file for {self.chianti_ion.Spectroscopic}')
                    print(' skipping...')
                return np.zeros((Ntemp, Nwvl), dtype=np.float64)

        # Load recombined ion data (rFblvl)
        if hasattr(continuum_spectrum, 'rFblvl'):
            rFblvl = continuum_spectrum.rFblvl
        else:
            lower = continuum_spectrum.nameDict['lower']
            lowerDict = chianti_util.convertName(lower)
            rFblvlname = lowerDict['filename'] + '.fblvl'
            if os.path.isfile(rFblvlname):
                continuum_spectrum.rFblvl = io.fblvlRead(lower)
                rFblvl = continuum_spectrum.rFblvl
            else:
                if self.verbose:
                    print(f' no Fblvl file for {self.chianti_ion.Spectroscopic}')
                    print(' skipping...')
                return np.zeros((Ntemp, Nwvl), dtype=np.float64)

        # Extract information for the recombined ion
        nlvls = len(rFblvl['lvl'])
        pqn = np.asarray(rFblvl['pqn'], dtype='int64')
        l = rFblvl['l']
        ecm = rFblvl['ecm']
        multr = rFblvl['mult']
        mult = fblvl['mult']

        # Get revised Gaunt factors
        klgbfn = chdata.Klgbfn

        # Initialize arrays for calculations
        expfun = np.zeros((nlvls, Ntemp, Nwvl), dtype=np.float64)
        fbn = np.zeros((nlvls, Ntemp, Nwvl), dtype=np.float64)
        fbIntensity = np.zeros((nlvls, Ntemp, Nwvl), dtype=np.float64)
        ratg = np.zeros(nlvls, dtype=np.float64)
        mygf = np.zeros((nlvls, Nwvl))
        ratg[0] = float(multr[0]) / float(mult[0])
        iprLvlEv = continuum_spectrum.Ipr - const.invCm2Ev * ecm[0]
        iprLvlErg = const.ev2Erg * iprLvlEv
        edgeLvlAng = []

        hnu = const.h * const.c / (1.e-8 * wavelength)
        hnuEv = const.ev2Ang / wavelength

        # Constants for free-bound emission calculation
        K1 = 2. ** 4 * const.h * const.q ** 2 / (3. * np.sqrt(3.) * const.emass * const.c)
        K2 = 1. / (2. * const.emass * const.c ** 2)
        K3 = (1. / (const.h * const.c)) * (1. / (np.sqrt(2. * const.emass))) * (1. / (const.pi * const.kB) ** 1.5)
        K0 = 1.e-8 * K1 * K2 * K3

        if verner:
            # Use Verner-Yakovlev photoionization cross-sections
            lvl1 = 1
            continuum_spectrum.vernerCross(wavelength)
            ilvl = 0
            iprLvlEv = continuum_spectrum.Ipr - const.invCm2Ev * ecm[ilvl]
            edgeLvlAng.append(const.ev2Ang / iprLvlEv)

            for itemp in goodT:
                atemp = self.temperature[itemp]

                xponent = (iprLvlErg - hnu) / (const.kB * atemp)
                expfun[0, itemp] = np.where(xponent <= 0., np.exp(xponent), 0.)

                c1 = const.verner * ratg[0] * expfun[0, itemp] * continuum_spectrum.VernerCross / atemp ** 1.5

                fbn[0, itemp] = (const.h * const.c / (1.e-8 * wavelength)) ** 5 * c1

                fbIntensity[ilvl, itemp] = emission_measure[itemp] * fbn[ilvl, itemp]
        else:
            lvl1 = 0

        # Calculate Gaunt factors and emissivity for remaining levels
        for ilvl in range(lvl1, nlvls):
            pqnIdx = pqn[ilvl] - 1
            lIdx = l[ilvl]
            klgbf = klgbfn[pqnIdx]
            pe = klgbf['pe']
            gf = klgbf['klgbf'][lIdx]
            iprLvlEv = continuum_spectrum.Ipr - const.invCm2Ev * ecm[ilvl]
            edgeLvlAng.append(const.ev2Ang / iprLvlEv)
            iprLvlErg = const.ev2Erg * iprLvlEv

            # Scaled energy relative to the ionization potential
            scaledE = hnuEv / continuum_spectrum.Ipr

            tck = splrep(np.log(pe), np.log(gf), s=0)
            gflog = splev(np.log(scaledE), tck, der=0, ext=3)
            mygf[ilvl] = np.where(hnuEv >= iprLvlEv, np.exp(gflog), 0.)

            ratg[ilvl] = float(multr[ilvl]) / float(mult[0])  # Ratio of statistical weights

            for itemp in goodT:
                atemp = self.temperature[itemp]

                xponent = (iprLvlErg - hnu) / (const.kB * atemp)
                expfun = np.where(xponent <= 0., np.exp(xponent), 0.)

                fbn[ilvl, itemp] = K0 * hnu ** 2 * expfun * iprLvlErg ** 2 * ratg[ilvl] * mygf[ilvl] \
                                   / (atemp ** 1.5 * float(pqn[ilvl]))

                fbIntensity[ilvl, itemp] = emission_measure[itemp] * fbn[ilvl, itemp]

        # Sum up intensities and apply elemental abundance and ion fraction
        freebound_emission = fbIntensity.sum(axis=0)

        if self.verbose:
            print(f' {self.chianti_ion.Spectroscopic} free-bound emission calculation completed')

        # differential volume emission measure
        DVEM = self.ne * species_density * shell_volume
        # multiplying
        freebound_emission = freebound_emission * DVEM[:, np.newaxis]
        return freebound_emission



