"""
Functions needed for standard Python multiprocessing module mspectrum
"""
import numpy as np
import ChiantiPy

def doFfQ(inQ, outQ):
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
        #abund = inputs[3]
        #em = inputs[4]
        ff = ChiantiPy.core.continuum(species, temperature, abundance=None, em=None, verbose=0)
        ff.freeFree(wavelength, includeAbund=False, includeIoneq=False)
        outQ.put(ff.FreeFree)
    return


