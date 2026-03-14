from NebulaPy.tools import constants as const
import multiprocessing as mp
import copy
from datetime import datetime
from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np
from ChiantiPy.core import mspectrum
from NebulaPy.tools import util as util
import ChiantiPy.tools.mputil as mputil
import ChiantiPy.tools.util as chianti_util

from .ChiantiMultiProc import *

from NebulaPy.src.LineEmission import line_emission


from ChiantiPy.core import mspectrum
import ChiantiPy.tools.filters as chfilters

class spectrum:

    ######################################################################################
    #
    ######################################################################################
    def __init__(
            self,
            min_photon_energy, max_photon_energy, energy_point_count,
            elements,
            bremsstrahlung=False,
            freebound=False,
            lines=False,
            twophoton=False,
            filtername=None,
            filterfactor=None,
            allLines=True,
            verbose=True
            ):

        self.min_energy = min_photon_energy
        self.max_energy = max_photon_energy
        self.N_wvl = energy_point_count
        self.elements = elements
        self.bremsstrahlung = bremsstrahlung
        self.freebound = freebound
        self.lines = lines
        self.twophoton = twophoton
        self.filtername = filtername
        self.filterfactor = filterfactor
        self.allLines = allLines
        self.verbose = verbose

        if self.verbose:
            print(f" ---------------------------")
            print(" initiating X-ray spectrum calculation...")
            print(f" bremsstrahlung emission = {self.bremsstrahlung}")
            print(f" free-bound emission = {self.freebound}")
            print(f" line intensity = {self.lines}")
            print(f" two photon emission = {self.twophoton}")
            if not (self.bremsstrahlung or self.freebound or self.lines):
                util.nebula_exit_with_error(" no emission processes specified")

        self.xray_containter = {
            'min_energy': self.min_energy,
            'max_energy': self.max_energy,
            'energy_unit': 'keV'
        }
        self.setup()

    ######################################################################################
    #
    ######################################################################################
    def setup(self):
        self.min_wvl = const.kev2Ang / self.max_energy
        self.max_wvl = const.kev2Ang / self.min_energy

        self.min_wvl = const.kev2Ang / self.max_energy
        self.max_wvl = const.kev2Ang / self.min_energy
        self.xray_containter['min_wvl'] = self.min_wvl
        self.xray_containter['max_wvl'] = self.max_wvl
        self.xray_containter['wvl_unit'] = 'Angstrom'
        self.wavelength = np.linspace(self.min_wvl, self.max_wvl, self.N_wvl)
        self.xray_containter['wvl_array'] = self.wavelength


        # Initialize the chianti object with the given elements, temperature, and electron density.
        chianti_obj = chianti(
            pion_elements=self.elements,
            temperature=[1.e+7], #dummy temperatute
            ne=[1.e+9], # dummy electron density
            verbose=self.verbose
        )

        # Update the xray_container with the species attributes.
        self.xray_containter.update(chianti_obj.species_attributes_container)
        self.species_attributes = chianti_obj.species_attributes_container



    ######################################################################################
    #
    ######################################################################################
    def xray_spectrum(self,
                       temperature,
                       density,
                       ne,
                       #elemental_abundances,
                       ion_fractions,
                       #shell_volume,
                       #dem_indices
                       ):
        """
        Calculate X-ray intensity for given temperature and electron density (ne).

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
        - Total X-ray intensity from selected emission types as a NumPy array.
        """

        #indices = [i for i, T in enumerate(temperature) if self.Tmin <= T < self.Tmax]
        #temperature = temperature[indices]
        #self.xray_containter['temperature'] = temperature
        #density = density[indices]
        #ne = ne[indices]
        #shell_volume = shell_volume[indices]
        #emission_measure = shell_volume

        #if len(elemental_abundances) != self.elements.size:
        #    util.nebula_exit_with_error('elemental abundance count does not match element count')

        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature)
        N_temp = len(temperature)  # Determine the number of temperature values.

        # Initialize empty arrays for storing X-ray intensity values.
        bremsstrahlung_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        freebound_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        line_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        twophoton_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)

        for species in self.species_attributes:

            # find the element the species belong to
            element = self.species_attributes[species]['Element']
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


        ion_density = [1.0]

        chianti_ion = chianti(chianti_ion='fe_25', temperature=temperature, ne=ne,
                              pion_ion=None, pion_elements=None, verbose=True)

        #bremsstrahlung_emission = chianti_ion.get_bremsstrahlung_emission(self.wavelength, ion_density)

        line_emission = chianti_ion.get_line_emission(self.wavelength)

        #spectrum = {"spectrum": bremsstrahlung_emission, "wavelength": self.wavelength}
        #spectrum = {"spectrum": line_emission, "wavelength": self.wavelength}
        return spectrum
