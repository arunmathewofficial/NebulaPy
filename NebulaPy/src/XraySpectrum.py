from NebulaPy.tools import constants as const
import multiprocessing as mp
import copy
from datetime import datetime
from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np
from ChiantiPy.core import mspectrum

import ChiantiPy.tools.mputil as mputil
import ChiantiPy.tools.util as util

from .ChiantiMultiProc import *

from NebulaPy.src.LineEmission import line_emission


from ChiantiPy.core import mspectrum
import ChiantiPy.tools.filters as chfilters

class xray:

    ######################################################################################
    #
    ######################################################################################
    def __init__(
            self,
            min_photon_energy, max_photon_energy, energy_point_count,
            elements,
            verbose=True
            ):

        self.min_energy = min_photon_energy
        self.max_energy = max_photon_energy
        self.N_wvl = energy_point_count
        self.elements = elements
        self.verbose = verbose
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
        self.xray_containter['min_wvl'] = self.min_wvl
        self.xray_containter['max_wvl'] = self.max_wvl
        self.xray_containter['wvl_unit'] = 'Angstrom'
        self.wavelength = np.linspace(self.min_wvl, self.max_wvl, self.N_wvl)
        self.xray_containter['wvl_array'] = self.wavelength

    ######################################################################################
    #
    ######################################################################################
    def xray_intensity(self, temperature, ne,

                       multiprocessing=False, ncores=None):

        # Initialize the chianti object with the given elements, temperature, and electron density.
        # generate attributes for the specified elements and stored in species attributes container
        chianti_obj = chianti(
            element_list=self.elements,
            temperature=temperature,
            ne=ne,
            verbose=self.verbose
        )

        # Update the xray_container with the species attributes.
        self.xray_containter.update(chianti_obj.species_attributes_container)


        for species in chianti_obj.species_attributes_container:
            print(species)



        '''
        chianti_obj.get_line_spectrum(ion, temperature, ne, wavelength, elemental_abundance=None,
                          ionization_fraction=None, emission_measure=None, filter=None)




        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature)
        N_temp = len(temperature)  # Determine the number of temperature values.

        # wavelength range
        wvl_range = [self.wavelength[0], self.wavelength[-1]]
        '''



