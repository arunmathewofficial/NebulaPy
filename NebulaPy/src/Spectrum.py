from NebulaPy.tools import constants as const
#import multiprocessing as mp
#import copy
#from datetime import datetime
#from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np
#from ChiantiPy.core import mspectrum
from NebulaPy.tools import util as util
#import ChiantiPy.tools.mputil as mputil
#import ChiantiPy.tools.util as chianti_util

#from .ChiantiMultiProc import *

from NebulaPy.src.LineEmission import line_emission


from ChiantiPy.core import mspectrum
import ChiantiPy.tools.filters as chfilters

class spectrum:

    ######################################################################################
    #
    ######################################################################################
    def __init__(
            self,
            min_wavelength=None,
            max_wavelength=None,
            min_photon_energy=None,
            max_photon_energy=None,
            N_point=1000,
            bremsstrahlung=False,
            freebound=False,
            lines=False,
            twophoton=False,
            filtername=None,
            filterfactor=None,
            allLines=True,
            verbose=True
    ):

        # flags and parameters
        self.bremsstrahlung = bremsstrahlung
        self.freebound = freebound
        self.lines = lines
        self.twophoton = twophoton
        self.filtername = filtername
        self.filterfactor = filterfactor
        self.allLines = allLines
        self.verbose = verbose

        if not (bremsstrahlung or freebound or lines or twophoton):
            util.nebula_exit_with_error("No emission processes specified")

        # wavelength and photon energy inputs
        if min_wavelength is not None and max_wavelength is not None:
            self.min_wvl, self.max_wvl = min_wavelength, max_wavelength

        elif min_photon_energy is not None and max_photon_energy is not None:
            self.min_wvl = const.kev2Ang / max_photon_energy
            self.max_wvl = const.kev2Ang / min_photon_energy

        else:
            util.nebula_exit_with_error(
                "Provide either wavelength range or photon energy range."
            )

        if self.min_wvl <= 0.0 or self.max_wvl <= 0.0:
            util.nebula_exit_with_error("Wavelength bounds must be positive.")

        if self.min_wvl >= self.max_wvl:
            util.nebula_exit_with_error(
                "Minimum wavelength must be smaller than maximum wavelength."
            )

        self.N_wvl = N_point
        self.wavelength = np.linspace(self.min_wvl, self.max_wvl, self.N_wvl)

        # Verbose output
        if self.verbose:
            print("--- Spectrum Calculation --------------------------------")
            print(f" Bremsstrahlung : {self.bremsstrahlung} | Free-bound : {self.freebound}")
            print(f" Lines          : {self.lines} | Two-photon : {self.twophoton}")
            print(f" Grid points    : {self.N_wvl}")
            print("---------------------------------------------------------")



    ######################################################################################
    # Build Species Attributes
    ######################################################################################
    def build_species_attributes(self, elements, elemental_abundances):
        """
        Build species attributes and validate elemental abundances.

        Parameters
        ----------
        elements : list
            Elements used in the spectral calculation.
        elemental_abundances : dict
            Elemental abundances given as mass fractions. Must sum to unity.
        """

        # Validate abundance input
        if not isinstance(elemental_abundances, dict):
            util.nebula_exit_with_error(
                "Provide elemental abundances (mass fractions)"
            )

        total_abundance = sum(elemental_abundances.values())

        # Check all elements exist
        missing = set(elements) - set(elemental_abundances.keys())
        if missing:
            util.nebula_exit_with_error(
                f"Missing abundances for elements {', '.join(missing)}"
            )

        if not np.isclose(total_abundance, 1.0, rtol=1e-8):
            util.nebula_exit_with_error(
                f"Elemental abundances must sum to unity (sum = {total_abundance:.8f})"
            )

        # Check for negative abundances
        for elem, val in elemental_abundances.items():
            if val < 0.0:
                util.nebula_exit_with_error(
                    f"Negative abundance detected for element '{elem}'"
                )

        # Initialize CHIANTI object (dummy plasma state)
        chianti_spec = chianti(
            pion_elements=elements,
            temperature=[1.0e7],  # dummy temperature
            ne=[1.0e9],  # dummy density
            verbose=self.verbose
        )
        # Return chianti ion attributes for the species
        self.chianti_species_attributes = chianti_spec.species_attributes_container
        self.elemental_abundances = elemental_abundances

        # do not terminate rather del
        del chianti_spec




    ######################################################################################
    # compute spectrum for a cell or set of cells
    # todo: is this name appropritate for this fucntion
    ######################################################################################
    def compute_emission(self,
                       temperature,
                       density,
                       ne,
                       #ion_fractions,
                       #shell_volume,
                       #dem_indices
                       ):
        """
        # todo:
        compute for a set of cells. the parameters should be arrays of the same length,
         each element corresponds to a cell. the output is also an array of the same length,
          each element is the spectrum of the corresponding cell.

        Parameters:
        - temperature: List or array of temperatures (in K) for the calculation.
        - ne: Electron density (in cm^-3).
        - freefree: Calculate free-free emission if True.
        - freebound: Calculate free-bound emission if True.
        - lines: Calculate line emission if True.
        - twophoton: Calculate two-photon emission if True.
        - multiprocessing: Use multiprocessing for the calculation if True.
        - ncores: Number of processor cores to use in multiprocessing (default is 3).

        Returns:
        # todo: return what?
        """

        # 👉 Numba for inner loops
        # 👉 Multiprocessing for outer cell chunks

        # info: there is no multiprocessing calculation for this method, rather make use of Numba. Will
        #  use multiprocessing to parallelize over cells.

        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature, dtype=np.float64)
        density = np.asarray(density, dtype=np.float64)
        ne = np.asarray(ne, dtype=np.float64)

        # Determine the number of temperature values
        N_temp = len(temperature)

        # Initialize empty arrays for ???.
        bremsstrahlung_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        freebound_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        line_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        twophoton_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)

        for species in self.chianti_species_attributes:

            if species != 'o_8':
                continue


            # find the element the species belong to
            element = self.chianti_species_attributes[species]['Element']
            # position of the corresponding element of the species in silo elements array
            #pos = np.where(self.elements == element)[0][0]
            # charge of the species
            #q = self.species_attributes[species]['Zion']
            # atomic mass of the species
            #Z = self.species_attributes[species]['Z']

            # calculating species density
            # if the species is top ion
            #if q == Z:
            #    top_ion = elemental_abundances[pos]
            #    for i in range(Z):
            #        top_ion -= ion_fractions[pos][i]
            #    species_num_density = density * top_ion / const.mass[element]
            # otherwise
            #else:
            #    species_num_density = density * ion_fractions[pos][q] / const.mass[element]

            # todo: should supply the right ion fraction for the species

            # todo: need to create object for each species and call the corresponding method to calculate
            #  the emission for the species, then sum over all species to get the total emission.
            #  This is because the line emission is calculated for each species separately in CHIANTI,
            #  while the bremsstrahlung and free-bound emissions are calculated for the plasma as a whole.

            species_processes = self.chianti_species_attributes[species]['keys']

            CHIANTI = chianti(
                chianti_ion=species,
                temperature=temperature,
                ne=ne,
                verbose=self.verbose
            )

            # Bremsstrahlung (free-free)
            if self.bremsstrahlung and 'ff' in species_processes:
                bremsstrahlung_emission = CHIANTI.get_bremsstrahlung_emission_rate(
                    wavelength=self.wavelength
                )


            # Line emission
            #if self.lines and 'line' in species_processes:
            #    line_emission = CHIANTI.get_line_emission()

            CHIANTI.terminate()

        '''
        ion_density = [1.0]

        chianti_ion = chianti(chianti_ion='fe_25', temperature=temperature, ne=ne,
                              pion_ion=None, pion_elements=None, verbose=True)

        #bremsstrahlung_emission = chianti_ion.get_bremsstrahlung_emission(self.wavelength, ion_density)

        line_emission = chianti_ion.get_line_emission(self.wavelength)

        #spectrum = {"spectrum": bremsstrahlung_emission, "wavelength": self.wavelength}
        #spectrum = {"spectrum": line_emission, "wavelength": self.wavelength}
        '''


        self.intensity = bremsstrahlung_emission * 1.e+27

