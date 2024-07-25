from .Chianti import chianti


class emissionline():


    def __init__(self, ion):
        """
        only single ion is considered here
        """
        self.Ion = ion

    #todo: write a fucntion to get all wavelength of a ion

    def Luminosity_1DGrid(self, wavelength, temperature, ne, ns, cellvolume):
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

        ion = chianti(self.Ion, temperature, ne)
        print(ion.emissivity())











