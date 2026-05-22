from NebulaPy.tools import constants as const
#import multiprocessing as mp
#import copy
#from datetime import datetime
#from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np
from .CIE import cieMode
#from ChiantiPy.core import mspectrum
from NebulaPy.tools import util as util
#import ChiantiPy.tools.mputil as mputil
#import ChiantiPy.tools.util as chianti_util
import ChiantiPy.core as ch
import NebulaPy.src.Chianti as nebula_chianti
from tqdm import tqdm
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
            bremsstrahlung=False,
            freebound=False,
            line=False,
            twophoton=False,
            elements=None,
            elemental_abundances=None,
            CIE=False,
            filtername=None,
            filterfactor=None,
            allLines=True,
            user_grid=False,
            grid_size=1000,
            verbose=True
    ):

        # flags and parameters
        self.bremsstrahlung = bremsstrahlung
        self.freebound = freebound
        self.line = line
        self.twophoton = twophoton
        self.filtername = filtername
        self.filterfactor = filterfactor
        self.allLines = allLines
        self.verbose = verbose

        if self.verbose:
            print(" [SPECTRUM MODULE] : Initializing spectrum calculation")

        # Setting up spectrum dictionary to hold results
        self.spectrum_container = {}

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

        if not (bremsstrahlung or freebound or line or twophoton):
            util.nebula_exit_with_error("No emission processes specified")

        # Verbose output
        if self.verbose:
            print(" Active radiative processes :")
            print(f" a) Bremsstrahlung continuum : {'ON' if bremsstrahlung else 'OFF'}")
            print(f" b) Free-bound continuum     : {'ON' if freebound else 'OFF'}")
            print(f" c) Spectral line emission   : {'ON' if line else 'OFF'}")
            print(f" d) Two-photon emission      : {'ON' if twophoton else 'OFF'}")

        # Initialize empty attributes for species and elemental abundances
        self.chianti_species_attributes = {}
        self.build_species_attributes(elements, elemental_abundances)

        # setup wavelength grid #####################################
        self.wavelength_grid = []
        self.setup_wavelength_grid(self.min_wvl, self.max_wvl,
                                   user_grid=user_grid,
                                   grid_size=grid_size)


        self.CIE = CIE
        self.NEQ = not CIE
        if CIE:
            cie = cieMode(verbose=True)
            cie.load_cie()


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
    # SETUP WAVELENGTH GRID
    ######################################################################################
    def setup_wavelength_grid(self, min_wvl, max_wvl, user_grid=False, grid_size=None):

        if self.verbose:
            if not user_grid:
                print(" [WAVELENGTH GRID] : Setting up default CHIANTI grid")
            else:
                print(" [WAVELENGTH GRID] : Setting up uniform grid")

        if min_wvl >= max_wvl:
            util.nebula_exit_with_error(
                " Minimum wavelength must be smaller than maximum wavelength."
            )

        if not self.chianti_species_attributes:
            util.nebula_exit_with_error(
                " Species Attributes Container is not initialized or is empty."
            )

        dummy_temperature = [2.0e6]
        dummy_ne = [1.0e9]
        lines = []

        species_list = list(self.chianti_species_attributes.items())

        for species, attributes in tqdm(
                species_list,
                desc=" Retrieving CHIANTI wavelength grid",
                unit=" species",
                ncols=100,
                leave=True,
                disable=not self.verbose
        ):

            if 'line' not in attributes['keys']:
                continue

            chianti_nebula_object = nebula_chianti.chianti(
                chianti_ion=species,
                temperature=dummy_temperature,
                ne=dummy_ne,
                verbose=False
            )

            try:
                ion_lines = chianti_nebula_object.get_allLineTransitions()['wvl']

                if ion_lines is not None:
                    lines.extend(
                        np.asarray(ion_lines, dtype=np.float64).tolist()
                    )

            finally:
                del chianti_nebula_object

        # ------------------------------------------------------------------
        # CHIANTI line-based wavelength grid
        # ------------------------------------------------------------------
        if not user_grid:

            lines = np.asarray(lines, dtype=np.float64)

            if lines.size == 0:
                util.nebula_exit_with_error(
                    " No CHIANTI lines found for the selected species."
                )

            selected_lines = lines[
                (lines >= min_wvl) &
                (lines <= max_wvl)
                ]

            wavelength_grid = np.unique(
                np.concatenate((
                    selected_lines,
                    np.asarray(
                        [min_wvl, max_wvl],
                        dtype=np.float64
                    )
                ))
            )

        # ------------------------------------------------------------------
        # User-defined uniform wavelength grid
        # ------------------------------------------------------------------
        elif grid_size is not None:

            if grid_size < 2:
                util.nebula_exit_with_error(
                    " grid_size must be at least 2."
                )

            wavelength_grid = np.linspace(
                min_wvl,
                max_wvl,
                int(grid_size),
                dtype=np.float64
            )

        # ------------------------------------------------------------------
        # No valid grid option
        # ------------------------------------------------------------------
        else:
            util.nebula_exit_with_error(
                " Either use_chianti_grid=True or grid_size must be specified."
            )

        wavelength_grid.sort()

        if wavelength_grid.size == 0:
            util.nebula_exit_with_error(
                " No wavelength points found in the requested wavelength range."
            )

        self.wavelength_grid = wavelength_grid
        self.N_wvl = wavelength_grid.size

        # store wavelength grid in the container
        if not hasattr(self, 'spectrum_container'):
            self.spectrum_container = {}

        self.spectrum_container['wavelength_grid'] = wavelength_grid

        if self.verbose:
            print(" Wavelength grid summary")
            print(f" Minimum wavelength   : {self.wavelength_grid[0]:.6f} Å")
            print(f" Maximum wavelength   : {self.wavelength_grid[-1]:.6f} Å")
            print(f" Total spectral points: {self.N_wvl}")


    ######################################################################################
    # Compute DEM
    ######################################################################################
    def compute_DEM_1D(self, temperature, species_density, shell_volume):

        print(" [DEM] : Computing DEM for 1D data, not implemented yet")

    ######################################################################################
    # Compute spectrum for 1D data
    ######################################################################################
    def compute_spectrum_1D(self,
                       temperature,
                       ne,
                       species_density,
                       shell_volume,
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

        if not self.chianti_species_attributes:
            util.nebula_exit_with_error(" Species Attributes Container is not initialized or is empty.")

        # Ensure all arrays have identical shape and size
        required_arrays = {
            "temperature": temperature,
            "ne": ne,
            "species_density": species_density,
            "shell_volume": shell_volume,
        }
        reference_name, reference_array = next(iter(required_arrays.items()))
        reference_shape = np.shape(reference_array)
        reference_size = np.size(reference_array)

        for name, array in required_arrays.items():

            if np.shape(array) != reference_shape:
                util.nebula_exit_with_error(
                    f" Shape mismatch detected: '{name}' has shape {np.shape(array)}, "
                    f"expected {reference_shape}."
                )

            if np.size(array) != reference_size:
                util.nebula_exit_with_error(
                    f" Size mismatch detected: '{name}' has size {np.size(array)}, "
                    f"expected {reference_size}."
                )



        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature, dtype=np.float64)
        ne = np.asarray(ne, dtype=np.float64)

        # Determine the number of temperature values
        N_temp = len(temperature)

        # Initialize empty arrays for ???.
        bremsstrahlung_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        freebound_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        line_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        twophoton_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
        total_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)

        if self.NEQ:
            util.nebula_warning(
                " NEI mode is not implemented yet; using CIE instead."
            )

        elif self.CIE:
            util.nebula_warning(
                " Using collisional ionization equilibrium (CIE)."
            )


        # info: 1) looping over species to calculate the emission rate from each process.
        for species in self.chianti_species_attributes:

            #chianti_species = self.get_chianti_symbol(pion_species, make=False)

            Z = self.chianti_species_attributes[species]['Z']
            ionstage = self.chianti_species_attributes[species]['Ion']
            dielectronic = self.chianti_species_attributes[species]['Dielectronic']



            if species != 'fe_25':
                util.nebula_warning(f"Skipping {species} ...")
                continue
            util.nebula_info(f"Only {species} is calculated in this test...")




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
                    wavelength=self.wavelength_grid
                )

            # bound-free
            if self.freebound and 'fb' in species_processes:
                freebound_emission = CHIANTI.get_freebound_emission_rate(wavelength=self.wavelength_grid)

            # bound-free
            if self.line and 'line' in species_processes:
                # get_line_emission_rate returns emissivity of each lines
                line_emission = CHIANTI.get_line_emission_rate(wavelength=self.wavelength_grid)
                # dividing by ne to obtain the vale same as chianti returns.
                line_emission = line_emission / ne[:, None]

            # two-photon
            if self.twophoton and 'line' in species_processes:
                if (Z - ionstage) in [0, 1] and not dielectronic:
                    twophoton_emission = CHIANTI.get_twophoton_emission_rate(wavelength=self.wavelength_grid)

            CHIANTI.terminate()

            # info: 2) sum all processes
            species_emission = bremsstrahlung_emission + freebound_emission + line_emission + twophoton_emission

        self.spectrum_container['spectrum'] = species_emission
