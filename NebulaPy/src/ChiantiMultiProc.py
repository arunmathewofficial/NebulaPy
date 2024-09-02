from .Chianti import chianti


######################################################################################
#
######################################################################################
def do_freefree_Q(inQ, outQ):
    """
    Multiprocessing helper for `ChiantiPy.core.continuum.freeFree`

    Parameters
    -----------
    inQ : `~multiprocessing.Queue`
        Ion free-free emission jobs queued up by multiprocessing module
    outQ : `~multiprocessing.Queue`
        Finished free-free emission jobs
    """
    for inputs in iter(inQ.get, 'STOP'):
        ion = inputs[0]
        temperature = inputs[1]
        ne = inputs[2]
        wavelength = inputs[3]
        ion_density = inputs[4]
        shell_volume = inputs[5]
        verbose = inputs[6]

        chianti_ion = chianti(chianti_ion=ion, temperature=temperature, ne=ne,
                              pion_ion=None, pion_elements=None, verbose=verbose)

        bremsstrahlung_emission = chianti_ion.get_bremsstrahlung_emission(
            wavelength,
            ion_density,
            shell_volume,
        )
        outQ.put({'intensity': bremsstrahlung_emission})
    return

######################################################################################
# do free-bound queue
######################################################################################
def do_freebound_Q(inQ, outQ):
    """
    Multiprocessing helper for `ChiantiPy.core.continuum.freeBound`

    Parameters
    -----------
    inQ : `~multiprocessing.Queue`
        Ion free-bound emission jobs queued up by multiprocessing module
    outQ : `~multiprocessing.Queue`
        Finished free-bound emission jobs
    """
    for inputs in iter(inQ.get, 'STOP'):
        ion = inputs[0]
        temperature = inputs[1]
        ne = inputs[2]
        wavelength = inputs[3]
        ion_density = inputs[4]
        shell_volume = inputs[5]
        verbose = inputs[6]

        chianti_ion = chianti(chianti_ion=ion, temperature=temperature, ne=ne,
                              pion_ion=None, pion_elements=None, verbose=verbose)

        freebound_emission = chianti_ion.get_freebound_emission(
            wavelength,
            ion_density,
            shell_volume,
            verner=True
        )

        outQ.put({'intensity': freebound_emission})
    return
######################################################################################
#
######################################################################################
def do_line_emission_Q(inQ, outQ):
    """
    Processes line emission calculations for a given set of inputs and puts the results into an output queue.

    Args:
        input_arg (queue): A queue containing input data, where each item is a list or tuple with the following:
            - ion (str): The ion species to be analyzed.
            - temperature (float): The temperature at which the ion is found.
            - ne (float): The electron density.
            - wavelength (float or list): The wavelength(s) at which to calculate the line emission.
            - abundance (float): The abundance of the ion species.
            - ion_fraction (float): The fraction of the ion present in the species.
            - em (float): The emission measure.
            - filter_name (str): The name of the filter to apply (optional).
            - filter_factor (float): A factor to apply with the filter (optional).
            - allLines (bool): A flag to include all lines in the spectrum.

        out_arg (queue): A queue for storing the results. Each result is a dictionary with:
            - 'intensity' (array): The calculated line spectrum intensity for the given inputs.

    Workflow:
        - Continuously retrieves input data from the input queue until it encounters a 'STOP' signal.
        - For each set of inputs, initializes a Chianti ion object and calculates the line spectrum using the provided parameters.
        - Places the resulting spectrum intensity into the output queue.
        - Stops processing when the 'STOP' signal is received.

    """

    # Iterate over inputs retrieved from the input_arg queue
    # `input_arg.get` fetches the next set of inputs from the queue
    # The loop will terminate when 'STOP' is retrieved
    for inputs in iter(inQ.get, 'STOP'):
        # Unpack the inputs
        ion = inputs[0]
        temperature = inputs[1]
        ne = inputs[2]
        wavelength = inputs[3]
        ion_density = inputs[4]
        shell_volume = inputs[5]
        filter_name = inputs[6]
        filter_factor = inputs[7]
        allLines = inputs[8]

        # Create a Chianti ion object with the specified ion, temperature, and electron density
        chianti_ion = chianti(chianti_ion=ion, temperature=temperature, ne=ne, verbose=True)


        # Calculate the line spectrum using the Chianti ion object
        line_spectrum = chianti_ion.get_line_spectrum(
            wavelength,
            species_density=ion_density,
            shell_volume=shell_volume,
            filtername=filter_name,
            filterfactor=filter_factor,
            allLines=allLines
        )

        # Put the calculated line spectrum into the output queue
        # The result is stored as a dictionary with the key 'intensity'
        outQ.put({'intensity': line_spectrum})

    # The method returns after processing all inputs and receiving the 'STOP' signal
    return




