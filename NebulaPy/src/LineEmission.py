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
        self.ion = ion
        self.verbose = verbose
        self.line_emission_container = {}
        self.line_emission_container['ion'] = self.ion




    ######################################################################################
    # line luminosity
    ######################################################################################
    def lineluminosity_spherical(self, lines, temperature, ne, species_density, shell_volume):
        '''

        Parameters
        ----------
        lines
        temperature
        ne
        ns
        dV

        Returns
        -------

        '''

        ion = chianti(pion_ion=self.ion, temperature=temperature, ne=ne, verbose=self.verbose)
        self.line_emission_container['temperature'] = temperature
        self.line_emission_container['ne'] = ne
        self.line_emission_container['spectroscopic'] = ion.chianti_ion.Spectroscopic

        # if the line (wavelength) is given in string, get the corresponding
        # float value
        #if isinstance(line, str):
        #    line = const.wvl_dict[line]

        all_emissivity_data = ion.get_emissivity()
        allLines = all_emissivity_data['wvl']
        self.line_emission_container['allLines'] = allLines
        self.line_emission_container['lines'] = lines

        indices = []
        for line in lines:
            if self.verbose:
                print(f' identifying {line} Å from allLines of {ion.chianti_ion.Spectroscopic}')
            index = (np.abs(allLines - line)).argmin()
            tolerance = 10 ** -4
            if np.abs(allLines[index] - line) <= tolerance:
                if self.verbose:
                    print(f' line {line} Å found at index {index} in allLines')
                indices.append(index)
            else:
                util.nebula_exit_with_error('line not found in allLines')

        self.line_emission_container['line_indices'] = indices

        if self.verbose:
            print(f' retrieving emissivity values for {ion.chianti_ion.Spectroscopic} lines(s): {lines}')

        emissivity = np.asarray(all_emissivity_data['emiss'][indices])
        self.line_emission_container['emiss'] = emissivity

        # Calculating line Luminosity
        if self.verbose:
            print(f' calculating line luminosity for {ion.chianti_ion.Spectroscopic} lines(s): {lines}')
        luminosity = [4.0 * const.pi * np.sum(e * species_density * shell_volume) for e in emissivity]
        self.line_emission_container['luminosity'] = luminosity









