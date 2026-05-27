#import multiprocessing as mp
#import copy
#from datetime import datetime
#from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np
from .CIE import cieMode
#from ChiantiPy.core import mspectrum
from NebulaPy.src import Utils as utils
#import ChiantiPy.tools.mputil as mputil
#import ChiantiPy.tools.util as chianti_util
import ChiantiPy.core as ch
import NebulaPy.src.Chianti as nebula_chianti
import NebulaPy.src.Constants as const
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
            print(" [ SPECTRUM MODULE ] : Initializing spectrum calculation")

        # Setting up spectrum dictionary to hold results
        self.spectrum_container = {}

        # wavelength and photon energy inputs
        if min_wavelength is not None and max_wavelength is not None:
            self.min_wvl, self.max_wvl = min_wavelength, max_wavelength

        elif min_photon_energy is not None and max_photon_energy is not None:
            self.min_wvl = const.kev2Ang / max_photon_energy
            self.max_wvl = const.kev2Ang / min_photon_energy

        else:
            utils.nebula_exit_with_error(
                "Provide either wavelength range or photon energy range."
            )

        if self.min_wvl <= 0.0 or self.max_wvl <= 0.0:
            utils.nebula_exit_with_error("Wavelength bounds must be positive.")

        if self.min_wvl >= self.max_wvl:
            utils.nebula_exit_with_error(
                "Minimum wavelength must be smaller than maximum wavelength."
            )

        if not (bremsstrahlung or freebound or line or twophoton):
            utils.nebula_exit_with_error("No emission processes specified")

        # Verbose output
        if self.verbose:
            print(" Active radiative processes :")
            print(f" a) Bremsstrahlung continuum : {'ON' if bremsstrahlung else 'OFF'}")
            print(f" b) Free-bound continuum     : {'ON' if freebound else 'OFF'}")
            print(f" c) Spectral line emission   : {'ON' if line else 'OFF'}")
            print(f" d) Two-photon emission      : {'ON' if twophoton else 'OFF'}")

        # Initialize empty attributes for species and elemental abundances
        self.chianti_species_attributes = {}
        self.build_species_attributes(elements)

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
    def build_species_attributes(self, elements):
        """
        Build species attributes.

        Parameters
        ----------
        elements : list
            Elements used in the spectral calculation.
        """
        # Initialize CHIANTI object (dummy plasma state)
        chianti_spec = chianti(
            pion_elements=elements,
            temperature=[1.0e7],  # dummy temperature
            ne=[1.0e9],  # dummy density
            verbose=self.verbose
        )
        # Return chianti ion attributes for the species
        self.chianti_species_attributes = chianti_spec.species_attributes_container
        # do not terminate rather del
        del chianti_spec

    ######################################################################################
    # SETUP WAVELENGTH GRID
    ######################################################################################
    def setup_wavelength_grid(self, min_wvl, max_wvl, user_grid=False, grid_size=None):

        if self.verbose:
            if not user_grid:
                print(" [ WAVELENGTH GRID ] : Setting up default CHIANTI wavelength grid")
            else:
                print(" [ WAVELENGTH GRID ] : Setting up uniform wavelength grid")

        if min_wvl >= max_wvl:
            utils.nebula_exit_with_error(
                " Minimum wavelength must be smaller than maximum wavelength."
            )

        if not self.chianti_species_attributes:
            utils.nebula_exit_with_error(
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
                utils.nebula_exit_with_error(
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
                utils.nebula_exit_with_error(
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
            utils.nebula_exit_with_error(
                " Either use_chianti_grid=True or grid_size must be specified."
            )

        wavelength_grid.sort()

        if wavelength_grid.size == 0:
            utils.nebula_exit_with_error(
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
    # Compute species spectra
    ######################################################################################
    from tqdm import tqdm

    ######################################################################################
    # COMPUTE SPECIES SPECTRA
    ######################################################################################
    from tqdm import tqdm

    ######################################################################################
    # COMPUTE SPECIES SPECTRA
    ######################################################################################
    def compute_species_spectra(self, row_temperature, row_ne, progress_bar=False):

        # Determine the number of temperature values
        N_temp = len(row_temperature)

        species_spectra = {}

        #species_iterator = tqdm(
        #    self.chianti_species_attributes,
        #    desc=" Computing species spectra",
        #    #unit=" species",
        #    ncols=90,
        #    disable=not progress_bar
        #)

        # info: looping over species to calculate the emission rate from each process
        for species in species_iterator:

            Z = self.chianti_species_attributes[species]['Z']
            ionstage = self.chianti_species_attributes[species]['Ion']
            dielectronic = self.chianti_species_attributes[species]['Dielectronic']

            # reset process arrays for each species
            bremsstrahlung_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
            freebound_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
            line_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
            twophoton_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)


            #if species != 'c_5':
            #    util.nebula_warning(f"Skipping {species} ...")
            #    continue

            #util.nebula_info(f"Only {species} is calculated in this test...")


            # show current species name in progress bar
            #species_iterator.set_postfix_str(species)

            species_processes = self.chianti_species_attributes[species]['keys']

            CHIANTI = chianti(
                chianti_ion=species,
                temperature=row_temperature,
                ne=row_ne,
                verbose=False
            )

            # Bremsstrahlung/free-free emission
            if self.bremsstrahlung and 'ff' in species_processes:
                bremsstrahlung_emission = CHIANTI.get_bremsstrahlung_emission_rate(
                    wavelength=self.wavelength_grid
                )

            # Free-bound emission
            if self.freebound and 'fb' in species_processes:
                freebound_emission = CHIANTI.get_freebound_emission_rate(
                    wavelength=self.wavelength_grid
                )

            # Line emission
            if self.line and 'line' in species_processes:
                line_emission = CHIANTI.get_line_emission_rate(
                    wavelength=self.wavelength_grid
                )

                # dividing by ne to obtain the same value returned by CHIANTI
                line_emission = line_emission / row_ne[:, None]

            # Two-photon emission
            if self.twophoton and 'line' in species_processes:
                if (Z - ionstage) in [0, 1] and not dielectronic:
                    twophoton_emission = CHIANTI.get_twophoton_emission_rate(
                        wavelength=self.wavelength_grid
                    )

            CHIANTI.terminate()

            # sum over all processes for this species
            species_emission = (
                    bremsstrahlung_emission
                    + freebound_emission
                    + line_emission
                    + twophoton_emission
            )

            species_spectra[species] = species_emission

        return species_spectra



    def compute_species_spectra_test(self, row_temperature, row_ne, progress_bar=False):

        # Determine the number of temperature values
        N_temp = len(row_temperature)

        species_spectra = {}

        #species_iterator = tqdm(
        #    self.chianti_species_attributes,
        #    desc=" Computing species spectra",
        #    #unit=" species",
        #    ncols=90,
        #    disable=not progress_bar
        #)

        # info: looping over species to calculate the emission rate from each process
        for species in self.chianti_species_attributes.keys():

            Z = self.chianti_species_attributes[species]['Z']
            ionstage = self.chianti_species_attributes[species]['Ion']
            dielectronic = self.chianti_species_attributes[species]['Dielectronic']

            # reset process arrays for each species
            bremsstrahlung_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
            freebound_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
            line_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)
            twophoton_emission = np.zeros((N_temp, self.N_wvl), dtype=np.float64)


            if species != 'fe_25':
                utils.nebula_warning(f"Skipping {species} ...")

                continue

            #util.nebula_info(f"Only {species} is calculated in this test...")


            # show current species name in progress bar
            #species_iterator.set_postfix_str(species)

            species_processes = self.chianti_species_attributes[species]['keys']
            print(species_processes)
            CHIANTI = chianti(
                chianti_ion=species,
                temperature=row_temperature,
                ne=row_ne,
                verbose=True
            )

            # Bremsstrahlung/free-free emission
            if self.bremsstrahlung and 'ff' in species_processes:
                bremsstrahlung_emission = CHIANTI.get_bremsstrahlung_emission_rate(
                    wavelength=self.wavelength_grid
                )

            # Free-bound emission
            if self.freebound and 'fb' in species_processes:
                freebound_emission = CHIANTI.get_freebound_emission_rate(
                    wavelength=self.wavelength_grid
                )

            # Line emission
            if self.line and 'line' in species_processes:
                line_emission = CHIANTI.get_line_emission_rate(
                    wavelength=self.wavelength_grid
                )

                # dividing by ne to obtain the same value returned by CHIANTI
                line_emission = line_emission / row_ne[:, None]

            # Two-photon emission
            if self.twophoton and 'line' in species_processes:
                if (Z - ionstage) in [0, 1] and not dielectronic:
                    twophoton_emission = CHIANTI.get_twophoton_emission_rate(
                        wavelength=self.wavelength_grid
                    )

            CHIANTI.terminate()

            # sum over all processes for this species
            species_emission = (
                    bremsstrahlung_emission
                    + freebound_emission
                    + line_emission
                    + twophoton_emission
            )

            species_spectra[species] = species_emission

        return species_spectra

    ######################################################################################
    # Compute spectrum for 2D data
    ######################################################################################
    def generate_spectrum(self, temperature, ne, species_densities, grid_volume, grid_mask):
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
            utils.nebula_exit_with_error(" Species Attributes Container is not initialized or is empty.")

        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature, dtype=np.float64)
        ne = np.asarray(ne, dtype=np.float64)
        grid_volume = np.asarray(grid_volume, dtype=np.float64)
        grid_mask = np.asarray(grid_mask, dtype=np.float64)

        # Safety check ##################################################
        if (
                temperature.shape != ne.shape or
                temperature.shape != grid_volume.shape or
                temperature.shape != grid_mask.shape
        ):
            utils.nebula_exit_with_error(
                " DEM-2D input arrays have inconsistent shapes."
            )

        for species, species_density in species_densities.items():
            species_density = np.asarray(species_density, dtype=np.float64)
            if temperature.shape != species_density.shape:
                utils.nebula_exit_with_error(
                    f" DEM-2D input arrays have inconsistent shapes for species {species}."
                )

        # which ionisation fraction
        if self.NEQ:
            utils.nebula_warning(
                " NEI mode is not implemented yet; using CIE instead."
            )

        elif self.CIE:
            utils.nebula_warning(
                " Using collisional ionization equilibrium (CIE)."
            )
        N_grid_level = len(temperature)

        # uniform grid
        if N_grid_level == 1:
            utils.nebula_warning("Not implemented")

        # for multiple level grid
        else:
            for level in range(N_grid_level):
                N_row = len(temperature[level])
                for row in range(N_row):
                    row_temperature = temperature[level][row]
                    row_ne = ne[level][row]

                    species_spectra = self.compute_species_spectra_test(row_temperature, row_ne, progress_bar=True)

                    print(species_spectra)
                    exit(1)










