from NebulaPy.tools import constants as const
import multiprocessing as mp

from datetime import datetime
from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np
from ChiantiPy.core import mspectrum

import ChiantiPy.tools.mputil as mputil

from .ChiantiMultiProc import *


from ChiantiPy.core import mspectrum
import ChiantiPy.tools.filters as chfilters

class xray:

    ######################################################################################
    #
    ######################################################################################
    def __init__(
            self,
            min_photon_energy, max_photon_energy,
            elements,
            ncores=3,
            timeout=0.1,
            verbose=True
            ):

        self.min_energy = min_photon_energy
        self.max_energy = max_photon_energy
        self.elements = elements
        self.verbose = verbose
        self.ncores = ncores
        self.timeout = timeout
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
        self.wavelength = np.array([self.min_wvl, self.max_wvl])

    ######################################################################################
    #
    ######################################################################################

    def xray_intensity(self, temperature, ne):
        chianti_obj = chianti(
            element_list=self.elements,
            temperature=temperature,
            ne=ne,
            verbose=self.verbose
        )

        timeout=0.1

        chianti_obj.get_elements_attributes()
        species_attributes = chianti_obj.species_attributes
        self.xray_containter.update(species_attributes)

        temperature = np.array(temperature)  # Convert temperature list to a NumPy array
        self.wavelength = np.array(self.wavelength)  # Convert wavelength list to a NumPy array

        proc = min([self.ncores, mp.cpu_count()])

        freeFree = np.zeros((len(temperature), len(self.wavelength)), np.float64).squeeze()

        print(freeFree)

        ffWorkerQ = mp.Queue()
        ffDoneQ = mp.Queue()


        for species in species_attributes:

            if 'ff' in species_attributes[species]['keys']:
                ffWorkerQ.put((species, temperature, self.wavelength, 1, 1))

        ffWorkerQSize = ffWorkerQ.qsize()

        nCores = mp.cpu_count()
        proc = min(proc, nCores)

        ffProcesses = []
        for i in range(proc):
            p = mp.Process(target=doFfQ, args=(ffWorkerQ, ffDoneQ))
            p.start()
            ffProcesses.append(p)
        #       timeout is not necessary
        for p in ffProcesses:
            if p.is_alive():
                p.join(timeout=timeout)
        #
        for iff in range(ffWorkerQSize):
            thisFreeFree = ffDoneQ.get()
            freeFree += thisFreeFree['intensity']
        for p in ffProcesses:
            if not isinstance(p, str):
                p.terminate()


        return freeFree
