"""
Functions needed for standard Python multiprocessing module mspectrum
"""
import numpy as np
import ChiantiPy

from .Chianti import chianti

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
        species = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        ion_abundance = inputs[3]
        em = inputs[4]
        ff = ChiantiPy.core.continuum(species, temperature, abundance=ion_abundance, em=em, verbose=0)
        ff.freeFree(wavelength, includeAbund=False, includeIoneq=False)
        outQ.put(ff.FreeFree)
    return


def do_line_emission_Q(input_arg, out_arg):

    for inputs in iter(input_arg, 'STOP'):
        ion = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        ion_abundance = inputs[3]
        em = inputs[4]
        ion = chianti(ion=ion, temperature=temperature, ne=ne, verbose=self.verbose)
        result = ion.get_line_spectrum()
        out_arg.put(result)
    return



