import os
os.environ['XUVTOP']
import ChiantiPy
import ChiantiPy.core as ch
import ChiantiPy.tools.filters as chfilters
import ChiantiPy.tools.io as chio


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

    def __init__(self, ionName, Temp, ne):

        self.IonChiantiName = ionName
        self.Temperature = Temp
        self.electronDensity = ne
        self.setup()


    def setup(self):
        chinati_object = ch.ion(self.IonChiantiName, temperature=self.Temperature,
                                eDensity=self.electronDensity, pDensity='default',
                                radTemperature=None, rStar=None, abundance=None,
                                setup=True, em=None, verbose=0)

        self.ChiantiInstant = chinati_object



    def emissivity(self):
        self.ChiantiInstant.emiss(True)
        emissivity = self.ChiantiInstant.Emiss
        return emissivity
