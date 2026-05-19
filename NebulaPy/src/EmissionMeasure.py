import numpy as np
import math
from astropy import units as u
import time


pi = math.pi
m_p = 1.6726219e-24  # g
k = 1.38064852e-16  # erg/K
mu = 0.61  # mH


class emissionMeasure():


    ######################################################################################
    # initializing
    ######################################################################################
    def __init__(self, Tmin, Tmax, Nbins, verbose=False):

        self.verbose = verbose
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.Nbins = Nbins

    ######################################################################################
    # generate differential emission measure indices todo: not verified
    ######################################################################################
    def generate_DEM_indices(self, temperature):

        print(" [EMISSION MEASURE] : Generating DEM indices")

        # Calculate the logarithmic width of each bin
        bin_width = (np.log10(self.Tmax) - np.log10(self.Tmin)) / self.Nbins
        # Half the width of a bin for adjusting bin edges
        half_bin_width = bin_width / 2

        # Generate the logarithmically spaced temperature bin edges.
        # This will create Nbins+1 edges to define the boundaries of Nbins.
        temperature_edges = np.linspace(np.log10(self.Tmin), np.log10(self.Tmax), self.Nbins + 1)

        # Create temperature bins by pairing adjacent edges.
        # Each bin is represented as [bin_min, bin_max].
        temperature_bins = [[temperature_edges[i], temperature_edges[i + 1]] for i in range(self.Nbins)]

        # Calculate the midpoints of each bin for potential further use.
        # These midpoints are the average of the logarithmic bin edges.
        Tb = np.array([(bin[0] + bin[1]) / 2 for bin in temperature_bins])


        DEM_indices = np.zeros_like(temperature, dtype=np.float64)

        # Find corresponding DEM bin for each temperature
        log_temperature = np.log10(temperature)

        # Assign DEM bin indices
        for bin_idx in range(self.Nbins):
            bin_min = temperature_edges[bin_idx]
            bin_max = temperature_edges[bin_idx + 1]

            mask = (log_temperature >= bin_min) & (log_temperature < bin_max)

            DEM_indices[mask] = bin_idx

        self.bin_temperature = Tb
        self.DEM_indices = DEM_indices



    ######################################################################################
    # differential emission measure todo: not verified
    ######################################################################################
    def DEM(self, temperature, ne, species_density, shellvolume):
        """
        Calculate the differential emission measure (DEM) across temperature bins.

        Parameters:
        ----------
        dem_indices : list of numpy.ndarray
            List where each element is an array of indices corresponding to a temperature bin.
        ne : numpy.ndarray
            Array of electron densities corresponding to the temperature values.
        shellvolume : numpy.ndarray
            Array of shell volumes corresponding to the temperature values.

        Returns:
        -------
        DEM : numpy.ndarray
            Array of differential emission measure values for each temperature bin.
        """

        self.generate_DEM_indices(temperature)




    ###############################################################################
    def volume2D(self, xmax, xmin, ngrid):  # Calculates the volume of each cell in the image grid
        xmax = xmax
        xmin = xmin
        ngrid = ngrid

        # Calculate the size of each cell in the x, y, and z-direction:
        delta_z = (xmax[0] - xmin[0]) / ngrid[0]
        delta_R = (xmax[1] - xmin[1]) / ngrid[1]
        # Create a 2D array of zeros with the same dimensions as ngrid:
        v = np.zeros((ngrid[1], ngrid[0]))
        # Loop through each element in the Volume array:
        for ycells in range(ngrid[1]):
            rmin = ycells * delta_R
            rmax = (ycells + 1) * delta_R
            for xcells in range(ngrid[0]):
                v[ycells, xcells] = delta_z * np.pi * (rmax ** 2 - rmin ** 2)
        del delta_z
        del delta_R
        del xmax
        del xmin
        del ngrid
        return v


    ###############################################################################
    def DEM2D(self, density, temperature, ne, mask, ngrid, volume, mesh_edges_min, mesh_edges_max, temp_bin, hw):
        # Function to calculate the differential emission measure of the nebula
        # See Green et al. (2019) - Bubble Nebula - paper for details.
        # temp_bin: array of temperature bins in logspace
        # hw:       half-width of bin

        # todo: this function need to be optimized and rewritten.
        # info: this fucntion is now used inside Xray spectrum.py
        density = density
        temp = temperature
        mask = mask
        ngrid = ngrid
        lim_max = (mesh_edges_max * u.cm).value
        lim_min = (mesh_edges_min * u.cm).value
        vol = volume

        dem_bin_all = np.zeros(len(temp_bin))


        for j in range(len(density)):
            denj = density[j]
            tempj = temp[j]
            maskj = mask[j]
            # info: replaced with my code to calculate the cell volume and ne
            #volj = self.volume2D(lim_max[j], lim_min[j], ngrid)
            volj = vol[j]
            #nei = denj * maskj
            #nei = nei * 1.2 * 0.715 / m_p  # assume 1.2 electrons per H ion
            nei = ne[j]
            # These are the two arrays I need:
            vol_den = volj * nei * nei
            log_masktemp = np.log10(tempj * maskj)
            # If masktemp[i,j,k] is in one bin range w, then add vol_dens[i,j,k] to dem[w]
            for w in range(len(temp_bin)):
              pick = np.zeros_like(log_masktemp)
              pick[(log_masktemp>=(temp_bin[w]-hw))&(log_masktemp<(temp_bin[w]+hw))] = 1
              dem_bin_all[w] += np.sum(vol_den * pick)

        return {'dem_bin': dem_bin_all}

