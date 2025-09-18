


import NebulaPy.tools.util as util


class recombination_line():


    def __init__(self, ion, verbose=True):
        """
        only single ion is considered here
        """
        self.ion = ion
        self.verbose = verbose
        self.line_emission_container = {}
        self.line_emission_container['ion'] = self.ion

    ######################################################################################
    # H-alpha recombination line for 1D silo data
    ######################################################################################
    def halpha_reccombination_line_1D(self, temperature, ne, species_density):
        util.nebula_exit_with_error("Not implemented yet!")


    ######################################################################################
    # H-alpha recombination line for 2D silo data
    ######################################################################################
    def halpha_reccombination_line_2D(self, temperature, ne, species_density):
        util.nebula_exit_with_error("Not implemented yet!")