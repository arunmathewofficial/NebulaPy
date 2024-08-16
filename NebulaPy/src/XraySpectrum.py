from NebulaPy.tools import constants as const
import multiprocessing as mp

from datetime import datetime
from ChiantiPy.base import specTrails
from .Chianti import chianti
import numpy as np

from .ChiantiMultiProc import *

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

    def process_multiprocessing_tasks(self, tasks, result_aggregator, worker_function):
        # Determine the number of cores to use
        ncores = min(self.ncores, mp.cpu_count())

        # Multiprocessing setup
        workerQ = mp.Queue()
        doneQ = mp.Queue()

        # Populate the worker queue
        task_count = len(tasks)
        for task in tasks:
            workerQ.put(task)

        if task_count > 0:
            # Start processes
            processes = [mp.Process(target=worker_function, args=(workerQ, doneQ)) for _ in range(ncores)]
            for p in processes:
                p.start()

            # Wait for processes to complete
            for p in processes:
                p.join(timeout=self.timeout)

            # Aggregate results
            for _ in range(task_count):
                result = doneQ.get()
                result_aggregator(result)

            # No need to explicitly terminate processes as join() waits for them to finish

    def xray_intensity(self, temperature, ne):
        chianti_obj = chianti(
            element_list=self.elements,
            temperature=temperature,
            ne=ne,
            verbose=self.verbose
        )

        chianti_obj.get_elements_attributes()
        species_attributes = chianti_obj.species_attributes
        self.xray_containter.update(species_attributes)

        temperature = np.array(temperature)  # Convert temperature list to a NumPy array
        self.wavelength = np.array(self.wavelength)  # Convert wavelength list to a NumPy array

        print(self.xray_containter)

        freeFree = np.zeros((temperature.size, self.wavelength.size), dtype=np.float64)

        # Define tasks and result aggregation function
        tasks = [(species, temperature, self.wavelength) for species in species_attributes
                 if 'ff' in species_attributes[species]['keys']]

        def aggregate_result(result):
            nonlocal freeFree
            freeFree += result['intensity']  # Assuming 'intensity' is an array matching freeFree's shape

        # Use the common multiprocessing method
        self.process_multiprocessing_tasks(tasks, aggregate_result, doFfQ)

        print(freeFree)



