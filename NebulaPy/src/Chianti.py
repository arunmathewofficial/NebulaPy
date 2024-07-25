import os
os.environ['XUVTOP']
import ChiantiPy
import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.io as chio
import numpy as np

class chianti:
    """
    The class for calculating emission line spectrum.

    Parameters
    ----------
    Keyword arguments
    -----------------
    Examples
    --------
    >>> temperature = 1.e+9
    >>> ne = 1.e+4
    >>> ion = nebula.chianti('o_4', temperature, ne)
    >>> print(ion.emissivity())
    Notes
    -----
    References
    ----------
    """

    def __init__(self, ionName, temperature, ne, verbose):

        self.IonChiantiName = ionName
        self.Temperature = temperature
        self.electronDensity = ne
        self.Verbose = verbose
        self.setup()


    def setup(self):
        chinati_object = ch.ion(self.IonChiantiName, temperature=self.Temperature,
                                eDensity=self.electronDensity, pDensity='default',
                                radTemperature=None, rStar=None, abundance=None,
                                setup=True, em=None, verbose=self.Verbose)

        self.ChiantiInstant = chinati_object

    def get_allLines(self):
        """
        Retrieve all spectral lines associated with a specified ion
        :return: wave-length array
        """
        if self.Verbose:
            print(' retrieving all spectral lines of ', self.ChiantiInstant.Spectroscopic)
        wvl = np.asarray(self.ChiantiInstant.Wgfa['wvl'], np.float64)
        wvl = np.abs(wvl)
        return wvl

    def get_emissivity(self):
        """
        Retrieve the emissivity values for all spectral lines associated
        with a specified ion.

        :return: Dict Emiss. Emiss has several quantities, namely, ion,
        # wvl(angstrom), emissivity (ergs s^-1 str^-1), pretty1, pretty2.
        """
        if self.Verbose:
            print(' retrieving emissivity values for all spectral lines of', self.ChiantiInstant.Spectroscopic)
        self.ChiantiInstant.emiss(True)
        emissivity = self.ChiantiInstant.Emiss
        return emissivity
