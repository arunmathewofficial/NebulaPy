from .Chianti import chianti
import NebulaPy.tools.constants as const
import numpy as np
import NebulaPy.tools.util as util

class emissionline():

    ######################################################################################
    # get parameter in spherical grid setup
    ######################################################################################
    def __init__(self, ion, verbose):
        """
        only single ion is considered here
        """
        self.Ion = ion
        self.Verbose = verbose


    ######################################################################################
    # line luminosity
    ######################################################################################
    def lineluminosity_spherical(self, line, temperature, ne, ns, dV):
        """
        Calculate the free-free energy loss rate of an ion. The result is returned to the
        `free_free_loss` attribute.
        The free-free radiative loss rate is given by Eq. 5.15a of [107]_. Writing the numerical
        constant in terms of the fine structure constant :math:`\\alpha`,

        .. math::
        \\frac{dW}{dtdV} = \\frac{4\\alpha^3h^2}{3\pi^2m_e}\left(\\frac{2\pi k_B}{3m_e}\\right)^{1/2}Z^2T^{1/2}\\bar{g}_B

        where where :math:`Z` is the nuclear charge, :math:`T` is the electron temperature, and
        math:`\\bar{g}_{B}` is the wavelength-averaged and velocity-averaged Gaunt factor. The
        Gaunt factor is calculated using the methods of [103]_. Note that this expression for the
        loss rate is just the integral over wavelength of Eq. 5.14a of [107]_, the free-free emission, and
        is expressed in units of erg :math:`\mathrm{cm}^3\,\mathrm{s}^{-1}`.

        """

        ion = chianti(ionName=self.Ion, temperature=temperature, ne=ne, verbose=self.Verbose)

        self.LineLuminosity = {'ion':self.Ion, 'temperature':temperature, 'ne':ne}
        self.LineLuminosity['spectroscopicName'] = ion.ChiantiInstant.Spectroscopic

        #print(ion.get_allLines()) # not in use

        # if the line (wavelength) is given in string, get the corresponding
        # float value
        if isinstance(line, str):
            line = const.wvl_dict[line]

        all_emissivity_data = ion.get_emissivity()
        allLines = all_emissivity_data['wvl']
        self.LineLuminosity['allLines'] = allLines
        self.LineLuminosity['line'] = line

        if self.Verbose:
            print(f' identifying {line} Å from allLines of {ion.ChiantiInstant.Spectroscopic}')
        index = (np.abs(allLines - line)).argmin()
        tolerance = 10 ** -4
        if np.abs(allLines[index] - line) <= tolerance:
            if self.Verbose:
                print(f' line {line} Å found at index {index} in allLines')
        else:
            util.nebula_exit_with_error('line not found in allLines')
        self.LineLuminosity['lineIndex'] = index

        if self.Verbose:
            print(f' retrieving cell emissivity values for {ion.ChiantiInstant.Spectroscopic} {line}')

        emissivity = np.asarray(all_emissivity_data['emiss'][index])
        self.LineLuminosity['emiss'] = emissivity
        self.LineLuminosity['ns'] = ns
        self.LineLuminosity['dV'] = dV

        # Calculating line Luminosity
        luminosity = np.sum(emissivity * ns * dV)

        self.LineLuminosity['luminosity'] = luminosity













