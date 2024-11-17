

class emission_measure():


    ######################################################################################
    # initializing
    ######################################################################################
    def __init__(self):
        pass

    ######################################################################################
    # differential emission measure
    ######################################################################################
    def DEM(self, dem_indices, ne, shellvolume):
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

        # Calculate ne * ne * dV for all elements
        volume_ne_square = ne * ne * shellvolume

        # Initialize an array for DEM with the same length as dem_indices
        DEM = np.zeros(len(dem_indices))

        # Use list comprehension and NumPy's array indexing to sum the values for each bin
        for i, indices in enumerate(dem_indices):
            if indices.size > 0:
                DEM[i] = np.sum(volume_ne_square[indices])

        # Remove zero values from DEM
        DEM = DEM[DEM > 0]
        return DEM

    ######################################################################################
    # generate differential emission measure indices
    ######################################################################################
    def generate_dem_indices(self, temperature, Tmin, Tmax, Nbins):
        # Calculate the logarithmic width of each bin
        bin_width = (np.log10(Tmax) - np.log10(Tmin)) / Nbins
        # Half the width of a bin for adjusting bin edges
        half_bin_width = bin_width / 2

        # Generate the logarithmically spaced temperature bin edges.
        # This will create Nbins+1 edges to define the boundaries of Nbins.
        temperature_edges = np.linspace(np.log10(Tmin), np.log10(Tmax), Nbins + 1)

        # Create temperature bins by pairing adjacent edges.
        # Each bin is represented as [bin_min, bin_max].
        temperature_bins = [[temperature_edges[i], temperature_edges[i + 1]] for i in range(Nbins)]

        # Calculate the midpoints of each bin for potential further use.
        # These midpoints are the average of the logarithmic bin edges.
        Tb = np.array([(bin[0] + bin[1]) / 2 for bin in temperature_bins])

        # Initialize an empty list to store indices of temperatures within each bin.
        dem_indices = []

        # Loop through each bin to identify temperature values that fall within the bin's range.
        for i in range(Nbins):
            # Find the indices of temperature values that fall within the current bin.
            # The condition checks if the logarithm of the temperature is within the bin range,
            # slightly adjusted by half_bin_width to ensure proper capturing of boundary values.
            indices = np.where((np.log10(temperature) >= temperature_bins[i][0] - half_bin_width) &
                               (np.log10(temperature) < temperature_bins[i][1] + half_bin_width))[0]

            # Append the array of indices to the dem_indices list.
            dem_indices.append(indices)

        # Filter Tb values for which dem_indices[i] is not empty
        Tb = [Tb[i] for i in range(Nbins) if len(dem_indices[i]) > 0]

        # Return the list of indices for further processing or analysis.
        return {'indices': dem_indices, 'Tb': Tb}

