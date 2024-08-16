from .Chianti import chianti
import NebulaPy.tools.constants as const
import numpy as np
import NebulaPy.tools.util as util

class line_emission():

    ######################################################################################
    # initializing the class emissionline
    ######################################################################################
    def __init__(self, ion, verbose):
        """
        only single ion is considered here
        """
        self.Ion = ion
        self.Verbose = verbose
        self.line_emission_container = {}



    ######################################################################################
    # line luminosity
    ######################################################################################
    def lineluminosity_spherical(self, line, temperature, ne, ns, dV):
        '''

        Parameters
        ----------
        line
        temperature
        ne
        ns
        dV

        Returns
        -------

        '''

        ion = chianti(ion=self.Ion, temperature=temperature, ne=ne, verbose=self.Verbose)
        self.line_emission_container['ion'] = self.Ion
        self.line_emission_container['temperature'] = temperature
        self.line_emission_container['ne'] = ne
        self.line_emission_container['spectroscopicName'] = ion.ChiantiInstant.Spectroscopic

        #print(ion.get_allLines()) # not in use

        # if the line (wavelength) is given in string, get the corresponding
        # float value
        if isinstance(line, str):
            line = const.wvl_dict[line]

        all_emissivity_data = ion.get_emissivity()
        allLines = all_emissivity_data['wvl']
        self.line_emission_container['allLines'] = allLines
        self.line_emission_container['line'] = line

        if self.Verbose:
            print(f' identifying {line} Å from allLines of {ion.ChiantiInstant.Spectroscopic}')
        index = (np.abs(allLines - line)).argmin()
        tolerance = 10 ** -4
        if np.abs(allLines[index] - line) <= tolerance:
            if self.Verbose:
                print(f' line {line} Å found at index {index} in allLines')
        else:
            util.nebula_exit_with_error('line not found in allLines')
        self.line_emission_container['lineIndex'] = index

        if self.Verbose:
            print(f' retrieving cell emissivity values for {ion.ChiantiInstant.Spectroscopic} {line}')

        emissivity = np.asarray(all_emissivity_data['emiss'][index])
        self.line_emission_container['emiss'] = emissivity

        # Calculating line Luminosity
        luminosity = 4.0 * const.pi * np.sum(emissivity * ns * dV)
        self.line_emission_container['luminosity'] = luminosity













