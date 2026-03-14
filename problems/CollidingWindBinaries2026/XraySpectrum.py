import numpy as np

import NebulaPy.src as nebula
#from NebulaPy.tools import util
#import time
#from pypion.ReadData import ReadData
#import matplotlib.pyplot as plt
#import numpy as np
#import astropy.units as unit
#import os
#from NebulaPy.tools import constants as const
#import pandas as pd

elements = ['Fe']

spectrum = nebula.spectrum(
    min_photon_energy=1.0,  # Minimum photon energy in keV
    max_photon_energy=10.0,  # Maximum photon energy in keV
    energy_point_count=4000,
    elements=elements,
    bremsstrahlung=True,
    freebound=False,
    lines=False,
    twophoton=False,
    verbose=True
)


'''
runtime = 0.0
# Loop over each time instant in the batched silo files
for step, silo_instant in enumerate(batched_silos):
    silo_instant_start_time = time.time()  # Record the start time

    # Read the data from the current silo file
    dataio = ReadData(silo_instant)
    basic = dataio.get_1Darray('Density')  # Retrieve basic simulation data, such as density
    sim_time = (basic['sim_time'] * unit.s).to(unit.kyr)  # Convert simulation time to kiloyears
    dataio.close()  # Close the data file

    # Print the current time instant
    print(f" ---------------------------")
    print(f" step: {step} | simulation time: {sim_time:.6e}")

    temperature = pion.get_parameter('Temperature', silo_instant)
    density = pion.get_parameter('Density', silo_instant)
    ne = pion.get_ne(silo_instant)
    elemental_mass_fraction = pion.get_elemental_mass_frac(silo_instant)
    ion_fractions = pion.get_tracer_values(silo_instant)

    Tmin = 1e+5
    Tmax = 1e+9

    ##
    #hack
    relevant_grid_points = np.where((temperature >= Tmin) & (temperature <= Tmax))[0]
    temperature = temperature[relevant_grid_points]
    density = density[relevant_grid_points]
    ne = ne[relevant_grid_points]
    shell_volume = shell_volume_original[relevant_grid_points]
    filtered_elem_mass_frac = []
    for e in range(len(elemental_mass_fraction)):
        # Ensure the slice operation is valid for each elemental_mass_fraction[e]
        filtered_elem_mass_frac.append(elemental_mass_fraction[e][relevant_grid_points])

    filtered_ion_frac = []
    for e in range(len(ion_fractions)):
        element_wise = []
        for i in range(len(ion_fractions[e])):
            element_wise.append(ion_fractions[e][i][relevant_grid_points])
        filtered_ion_frac.append(element_wise)
    ##

    dem = pion.generate_dem_indices(temperature=temperature, Tmin=Tmin, Tmax=Tmax, Nbins=100)
    dem_indices = dem['indices']

    spectrum = xray_emission.xray_intensity(
        temperature=temperature,
        density=density, ne=ne,
        elemental_abundances=filtered_elem_mass_frac,
        ion_fractions=filtered_ion_frac,
        shell_volume=shell_volume,
        dem_indices=dem_indices
    )


    silo_instant_finish_time = time.time()  # Record the finish time
    # Calculate the time spent on the current step
    dt = silo_instant_finish_time - silo_instant_start_time
    # Update the runtime with the time spent on the current step
    runtime += dt
    print(f" runtime: {runtime:.4e} s | dt: {dt:.4e} s")

    # calculate differential emission measure
    DEM = pion.DEM(
        dem_indices=dem_indices,
        ne=ne,
        shellvolume=shell_volume,
    )

    # making DEM plot
    dem_filename = filebase + f"_t{int(sim_time.value)}_dem.png"
    dem_file = os.path.join(output_path, dem_filename)
    plt.plot(dem['Tb'], np.log10(DEM), linestyle='-', color='b', label=f'time = {sim_time.value:.4e} kyr')
    plt.xlabel(r'log(T$_b$) K', fontsize=14)
    plt.ylabel(r'log(DEM) cm$^{-3}$', fontsize=14)
    plt.legend(fontsize=14, frameon=False)
    plt.savefig(dem_file)  # Save as a PNG file
    plt.close()  # Close the plot to free memory
'''


'''
temperature = [3.e+7]
density = [1.e+9]
ne = [100]
ion_fraction = [1.0]




spectrum = spectrum.xray_spectrum(temperature=temperature, density=density, ne=ne,
                                        # #elemental_abundances=filtered_elem_mass_frac,
                                        ion_fractions=ion_fraction,
                                        #shell_volume=shell_volume,
                                        #dem_indices=dem_indices
                                        )


wavelength = spectrum["wavelength"]
emission = spectrum["spectrum"]

print(emission)

import matplotlib.pyplot as plt
# Plot
plt.figure(figsize=(8, 5))

plt.plot(wavelength, emission, linewidth=2)

plt.xlabel("Wavelength (Å)")
plt.ylabel("Bremsstrahlung Emission (erg cm$^{-3}$ s$^{-1}$ Å$^{-1}$)")
plt.title("Free-Free (Bremsstrahlung) Spectrum")

plt.yscale("log")   # spectra usually span many orders
plt.grid(True)

plt.tight_layout()
plt.show()

'''


import ChiantiPy.core as ch
import ChiantiPy.tools.data as chdata



fe25 = ch.ion('fe_25', temperature=[2.e+6], eDensity=1.e+9, em=1.e+27)




wvl = 200 + 0.125 * np.arange(801)


fe25.populate()
fe25.intensity()
fe25.spectrum(wavelength=wvl)

plotimport ChiantiPy.tools.data as chdata


fe25.intensityPlot(wvlRange=[210, 220.0])









