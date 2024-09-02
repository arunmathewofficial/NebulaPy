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

class xray:

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
            multiprocessing=False,
            ncores=None,
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
        self.multiprocessing = multiprocessing
        self.ncores = ncores
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

        if self.verbose and self.multiprocessing:
            print(f" multiprocessing with {self.ncores} cores")
            self.timeout = 0.1




    ######################################################################################
    #
    ######################################################################################
    def xray_intensity(self,
                       temperature, density, ne, elemental_abundances,
                       ion_fractions, shell_volume, dem_indices
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

        if len(elemental_abundances) != self.elements.size:
            util.nebula_exit_with_error('elemental abundance count does not match element count')

        # Convert the temperature list to a NumPy array for efficient numerical operations.
        temperature = np.array(temperature)
        N_temp = len(temperature)  # Determine the number of temperature values.

        # Initialize empty arrays for storing X-ray intensity values.
        bremsstrahlung_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        freebound_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        line_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)
        twophoton_spectrum = np.zeros((N_temp, self.N_wvl), np.float64)

        # If multiprocessing is enabled, set up parallel processing ##############################
        if self.multiprocessing:
            ncores = self.ncores or 3
            cpu_count = mp.cpu_count()
            proc = min(ncores, cpu_count)

            # Define worker and done queues for multiprocessing tasks.
            bremsstrahlung_workerQ, bremsstrahlung_doneQ = (mp.Queue(), mp.Queue())
            freebound_workerQ, freebound_doneQ = (mp.Queue(), mp.Queue())
            line_emission_workerQ, line_emission_doneQ = (mp.Queue(), mp.Queue())

            # Populate the worker queues with tasks for the species.
            for species in self.species_attributes:

                # find the element the species belong to
                element = self.species_attributes[species]['Element']
                # position of the corresponding element of the species in silo elements array
                pos = np.where(self.elements == element)[0][0]
                # charge of the species
                q = self.species_attributes[species]['Zion']
                # atomic mass of the species
                Z = self.species_attributes[species]['Z']

                # calculating species density
                # if the species is top ion
                if q == Z:
                    top_ion = elemental_abundances[pos]
                    for i in range(Z):
                        top_ion -= ion_fractions[pos][i]
                    species_num_density = density * top_ion / const.mass[element]
                # otherwise
                else:
                    species_num_density = density * ion_fractions[pos][q] / const.mass[element]


                # Processes Work Queue =====================================================
                if self.bremsstrahlung and 'ff' in self.species_attributes[species]['keys']:
                    bremsstrahlung_workerQ.put(
                        (species,
                         temperature,
                         ne,
                         self.wavelength,
                         species_num_density,
                         shell_volume,
                         self.verbose
                         )
                    )

                if self.freebound and 'fb' in self.species_attributes[species]['keys']:
                    freebound_workerQ.put(
                        (species,
                         temperature,
                         ne,
                         self.wavelength,
                         species_num_density,
                         shell_volume,
                         self.verbose
                         )
                    )

                if self.lines and 'line' in self.species_attributes[species]['keys']:
                    line_emission_workerQ.put(
                        (species,
                         temperature,
                         ne,
                         self.wavelength,
                         species_num_density,
                         shell_volume,
                         self.filtername,
                         self.filterfactor,
                         self.allLines
                         )
                    )
                # End of Processes Work Queue =============================================
            #exit(1)
            # Free-free emission calculation using multiprocessing.
            if self.bremsstrahlung:
                bremsstrahlung_processes = []
                bremsstrahlung_workerQSize = bremsstrahlung_workerQ.qsize()
                for i in range(proc):
                    p = mp.Process(target=do_freefree_Q, args=(bremsstrahlung_workerQ, bremsstrahlung_doneQ))
                    p.start()
                    bremsstrahlung_processes.append(p)

                for p in bremsstrahlung_processes:
                    if p.is_alive():
                        p.join(timeout=self.timeout)

                for index_brem in range(bremsstrahlung_workerQSize):
                    thisFreeFree = bremsstrahlung_doneQ.get()
                    bremsstrahlung_spectrum += thisFreeFree['intensity']

                for p in bremsstrahlung_processes:
                    if p.is_alive():
                        p.terminate()

            # Free-bound emission calculation using multiprocessing.
            if self.freebound:
                freebound_processes = []
                freebound_workerQSize = freebound_workerQ.qsize()
                for i in range(proc):
                    p = mp.Process(target=do_freebound_Q, args=(freebound_workerQ, freebound_doneQ))
                    p.start()
                    freebound_processes.append(p)

                for p in freebound_processes:
                    if p.is_alive():
                        p.join(timeout=self.timeout)

                for index_freebound in range(freebound_workerQSize):
                    thisFreeBound = freebound_doneQ.get()
                    freebound_spectrum += thisFreeBound['intensity']

                for p in freebound_processes:
                    if p.is_alive():
                        p.terminate()

            # line emission calculation using multiprocessing #################################
            if self.lines:
                line_emission_task = []
                line_emission_workerQSize = line_emission_workerQ.qsize()
                for i in range(proc):
                    p = mp.Process(target=do_line_emission_Q, args=(line_emission_workerQ, line_emission_doneQ))
                    p.start()
                    line_emission_task.append(p)

                for task in line_emission_task:
                        task.join(timeout=self.timeout)

                for index_line in range(line_emission_workerQSize):
                    this_line_emission = line_emission_doneQ.get()
                    if 'errorMessage' not in this_line_emission.keys():
                        line_spectrum += this_line_emission['intensity'].squeeze()

                for task in line_emission_task:
                    if task.is_alive():
                        task.terminate()
        # end of multiprocessing

        # summing up all processes
        total = bremsstrahlung_spectrum + freebound_spectrum + line_spectrum

        # Initialize an array for the spectrum for each value of Tb
        spectrum = []

        # Sum the values for each set of indices specified in dem_indices
        for indices in dem_indices:
            if len(indices) > 0:  # Ensure indices is not empty
                # sum the total across the specified indices
                bin_spectrum = np.sum(total[indices], axis=0)
                spectrum.append(bin_spectrum)

        return spectrum
        ###############################################################################################
