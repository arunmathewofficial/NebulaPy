import os
os.environ['XUVTOP']
import ChiantiPy
import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.io as chio
import ChiantiPy.tools.data as chdata
from ChiantiPy.base import specTrails
import numpy as np
import ChiantiPy.tools.util as util

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
    def __init__(self, temperature, ne, ion=None, element_list=None, verbose=False):

        self.temperature = temperature
        self.ne= ne
        self.verbose = verbose

        if ion is not None:
            self.chianti_ion_name = self.get_chianti_symbol(ion, make=0)
            self.chianti_ion = ch.ion(
                self.chianti_ion_name, temperature=self.temperature,
                eDensity=self.ne, pDensity='default',
                radTemperature=None, rStar=None, abundance=None,
                setup=True, em=None, verbose=self.verbose
            )

        if element_list is not None:
            chianti_element_list = []
            for element in element_list:
                element_symbol = self.get_chianti_symbol(element, make=1)
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
    # Get emissivity
    ######################################################################################
    def get_emissivity(self):
        """
        Retrieve the emissivity values for all spectral lines associated
        with a specified ion.

        :return: Dict Emiss. Emiss has several quantities, namely, ion,
        # wvl(angstrom), emissivity (ergs s^-1 str^-1), pretty1, pretty2.
        """
        if self.verbose:
            print(' retrieving emissivity values for all spectral lines '
                  'of', self.chianti_ion.Spectroscopic)
        self.chianti_ion.emiss(True)
        emissivity = self.chianti_ion.Emiss
        return emissivity

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
        print(f" ---------------------------")
        print(f" species attributes")
        for akey in sorted(species.Todo.keys()):
            self.species_attributes_container[akey] = util.convertName(akey)  # Convert the key and store it
            if self.verbose:
                print(f" retrieving attributes for {self.species_attributes_container[akey]['spectroscopic']}")
            self.species_attributes_container[akey]['keys'] = species.Todo[akey]  # Store relevant data

            # Remove unnecessary data from the dictionary
            del self.species_attributes_container[akey]['filename']
            del self.species_attributes_container[akey]['experimental']

        # Finalize the species attributes dictionary
        # At this point, `self.species_attributes` contains all the relevant
        # attributes for the species in `chianti_element_list`

    ######################################################################################
    # get line spectrum
    ######################################################################################
    def get_line_spectrum(self, wavelength, abun, ionfrac, emission_measure, allLines=True, filtername=None, filterfactor=None):
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
        Intensity has the units \( \text{ergs} \, \text{cm}^{-3} \, \text{s}^{-1} \, \text{str}^{-1} \).


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
        lines = self.chianti_ion.Emiss['wvl']
        N_lines = len(lines)

        # Initialize the intensity array
        intensity = np.zeros((N_temp, N_lines), dtype=np.float64)

        # Calculate intensity for each temperature
        for temp_idx in range(N_temp):
            intensity[temp_idx] = abun * ionfrac * emissivity[:, temp_idx] * emission_measure[temp_idx] / self.ne[temp_idx]

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
        selected_idx = util.between(lines, wvl_range)

        if len(selected_idx) == 0:
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

        return line_spectrum









