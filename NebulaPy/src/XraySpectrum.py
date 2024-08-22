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
                       freefree=False, freebound=False, lines=False, twophoton=False,
                       multiprocessing=False, ncores=None):
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

        # Initialize the chianti object with the given elements, temperature, and electron density.
        chianti_obj = chianti(
            element_list=self.elements,
            temperature=temperature,
            ne=ne,
            verbose=self.verbose
        )

        # Update the xray_container with the species attributes.
        self.xray_containter.update(chianti_obj.species_attributes_container)
        species_attributes = chianti_obj.species_attributes_container

        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature)
        N_temp = len(temperature)  # Determine the number of temperature values.

        # Wavelength range
        wvl_range = [self.wavelength[0], self.wavelength[-1]]

        # Initialize empty arrays for storing X-ray intensity values.
        freefree_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        freebound_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        line_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        twophoton_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)

        # If multiprocessing is enabled, set up parallel processing.
        if multiprocessing:
            ncores = ncores or 3
            cpu_count = mp.cpu_count()
            proc = min(ncores, cpu_count)
            timeout = 0.1

            # Define worker and done queues for multiprocessing tasks.
            freefree_workerQ, freefree_doneQ = (mp.Queue(), mp.Queue()) if freefree else (None, None)
            freebound_workerQ, freebound_doneQ = (mp.Queue(), mp.Queue()) if freebound else (None, None)
            line_emission_workerQ, line_emission_doneQ = (mp.Queue(), mp.Queue()) if lines else (None, None)

            abundance = None
            em = None
            filter = (chfilters.gaussianR, 1000.)
            allLines = True
            doContinuum = True

            # Populate the worker queues with tasks for the species.
            for species in species_attributes:
                if freefree and 'ff' in species_attributes[species]['keys']:
                    freefree_workerQ.put((species, temperature, self.wavelength, abundance, em))

                if freebound and 'fb' in species_attributes[species]['keys']:
                    freebound_workerQ.put((species, temperature, self.wavelength, abundance, em))

                if lines and 'line' in species_attributes[species]['keys']:
                    line_emission_workerQ.put(
                        (species,
                         temperature,
                         ne,
                         self.wavelength,
                         filter,
                         allLines,
                         abundance,
                         em,
                         doContinuum
                         )
                    )

            # Free-free emission calculation using multiprocessing.
            if freefree:
                freefree_processes = []
                freefree_workerQSize = freefree_workerQ.qsize()
                for i in range(proc):
                    p = mp.Process(target=do_freefree_Q, args=(freefree_workerQ, freefree_doneQ))
                    p.start()
                    freefree_processes.append(p)

                for p in freefree_processes:
                    if p.is_alive():
                        p.join(timeout=timeout)

                for index_freefree in range(freefree_workerQSize):
                    thisFreeFree = freefree_doneQ.get()
                    freefree_spectrum += thisFreeFree['intensity']

                for p in freefree_processes:
                    if p.is_alive():
                        p.terminate()

            # Free-bound emission calculation using multiprocessing.
            if freebound:
                freebound_processes = []
                freebound_workerQSize = freebound_workerQ.qsize()
                for i in range(proc):
                    p = mp.Process(target=mputil.doFbQ, args=(freebound_workerQ, freebound_doneQ))
                    p.start()
                    freebound_processes.append(p)

                for p in freebound_processes:
                    if p.is_alive():
                        p.join(timeout=timeout)

                for index_freebound in range(freebound_workerQSize):
                    thisFreeBound = freebound_doneQ.get()
                    if 'errorMessage' not in thisFreeBound.keys():
                        freebound_spectrum += thisFreeBound['intensity'].squeeze()

                for p in freebound_processes:
                    if p.is_alive():
                        p.terminate()

            # line emission calculation using multiprocessing.
            if lines:
                line_emission_task = []
                line_emission_workerQSize = line_emission_workerQ.qsize()
                for i in range(proc):
                    p = mp.Process(target=do_line_emission_Q, args=(line_emission_workerQ, line_emission_doneQ))
                    p.start()
                    line_emission_task.append(p)

                for task in line_emission_task:
                    if task.is_alive():
                        task.join(timeout=timeout)

                for index_line in range(line_emission_workerQSize):
                    this_line_emission = line_emission_doneQ.get()
                    if 'errorMessage' not in this_line_emission.keys():
                        line_spectrum += this_line_emission['intensity'].squeeze()

                for task in line_emission_task:
                    if task.is_alive():
                        task.terminate()

        # Return the calculated spectra.
        return freefree_spectrum + freebound_spectrum + line_spectrum
        ###############################################################################################
