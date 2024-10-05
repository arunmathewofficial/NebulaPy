# NebulaPy - Version 0.0.1

 [![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
 [![PyPI version](https://badge.fury.io/py/package-name.svg)](https://pypi.org/project/package-name/) 

## Overview
NebulaPy is a Python library that generates synthetic emission maps from PION simulation data, which utilize an advanced multi-ion non-equilibrium solver for ionized astrophysical plasmas with arbitrary elemental abundances. NebulaPy also supports energy binning of stellar atmosphere models, including Atlas, Potsdam, and CMF, across a wide range of metal abundances for radiative sources used in PION simulations.

 ## Key Features
 - **Feature 1**: 
- **Feature 2**: [Brief description of feature 2 and its advantage]
- **Feature 3**: [Any performance or scalability advantage]
- **Feature 4**: [Highlight if it supports multiprocessing or is highly customizable] 

## Installation You can install the package using
 `pip`: ```bash pip install package-name


# Wolf-Rayet model grids

| Metallicity | Composition | log L | vfinal | Dmax | XH   | XHe  | XC     | XN     | XO    | XNe    | XFe    |
|-------------|-------------|-------|--------|------|------|------|--------|--------|-------|--------|--------|
| MW          | WNE         | 5.3   | 1600   | 4    | None | None | 0.98   | 1.0E-4 | 0.015 | None   | 1.4E-3 |
| MW          | WNL-H20     | 5.3   | 1000   | 4    | 0.2  | 0.78 | 1.0E-4 | 0.015  | None  | None   | 1.4E-3 |
| MW          | WNL-H50     | 5.3   | 1000   | 4    | 0.5  | 0.48 | 1.0E-4 | 0.015  | None  | None   | 1.4E-3 |
| MW          | WC          | 5.3   | 2000   | 10   | None | 0.55 | 0.4    | None   | 0.05  | None   | 1.6E-3 |
| LMC         | WNE         | 5.3   | 1600   | 10   | None | 0.995| 7.0E-5 | 4.0E-3 | None  | None   | 7.0E-4 |
| LMC         | WNL-H20     | 5.3   | 1000   | 10   | 0.2  | 0.795| 7.0E-5 | 4.0E-3 | None  | None   | 7.0E-4 |
| LMC         | WNL-H40     | 5.3   | 1000   | 10   | 0.4  | 0.595| 7.0E-5 | 4.0E-3 | None  | None   | 7.0E-4 |
| LMC         | WC          | 5.3   | 2000   | 10   | None | 0.55 | 0.4    | None   | 0.05  | 1.0E-3 | 7.0E-4 |
| SMC         | WNE         | 5.3   | 1600‡  | 4    | None | 0.998| 2.5E-5 | 1.5E-3 | None  | None   | 3.0E-4 |
| SMC         | WNL-H20     | 5.3   | 1600‡  | 4    | 0.2  | 0.798| 2.5E-5 | 1.5E-3 | None  | None   | 3.0E-4 |
| SMC         | WNL-H40     | 5.3   | 1600‡  | 4    | 0.4  | 0.598| 2.5E-5 | 1.5E-3 | None  | None   | 3.0E-4 |
| SMC         | WNL-H60     | 5.3   | 1600‡  | 4    | 0.6  | 0.398| 2.5E-5 | 1.5E-3 | None  | None   | 3.0E-4 |
| SMC         | WC          | 5.3   | 2000   | 10   | None | 0.547| 0.4    | None   | 0.05  | 2.4E-3 | 3.0E-4 |
| Z0.07       | WNE         | 5.3   | 1600   | 10   | None | 0.999| 1.0E-5 | 6.1E-4 | 1.0E-5| None   | 9.2E-5 |
| Z0.07       | WNL-H20     | 5.3   | 1600   | 10   | 0.2  | 0.799| 1.0E-5 | 6.1E-4 | 1.0E-5| None   | 9.2E-5 |
| Z0.07       | WNL-H40     | 5.3   | 1600   | 10   | 0.4  | 0.599| 1.0E-5 | 6.1E-4 | 1.0E-5| None   | 9.2E-5 |
| Z0.07       | WC          | 5.3   | 2000   | 10   | None | 0.549| 0.4    | None   | 0.05  | 8.3E-4 | 9.2E-5 |



---

# ATLAS9 Model Parameters

The **ATLAS9** model grid provides stellar atmosphere models over a comprehensive range of metallicity, surface gravity, and effective temperature. This allows for the modeling of various stellar types, from cool dwarf stars to hot, massive stars, and everything in between. 

### Overview
The ATLAS9 model includes grids for the following parameters:

- **Metallicity [M/H]**: Metallicity reflects the abundance of elements heavier than hydrogen and helium relative to the solar composition. The grid spans a range of metallicity values, from highly metal-poor stars to stars with higher-than-solar metal content.
- **Gravity log(g)**: Surface gravity, a measure of the gravitational acceleration at the surface of a star, is represented in the grid with a wide range of values suitable for different types of stars, from giants (low gravity) to dwarfs (high gravity).
- **Effective Temperature (Teff)**: The grid covers a wide range of effective temperatures, from relatively cool stars to very hot stars, with uneven spacing to capture relevant changes in stellar atmospheres across the temperature spectrum.

### Model Parameters

| **Parameter**           | **Values**                                                                                                                                  |
|-------------------------|--------------------------------------------------------------------------------------------------------------------------------------------|
| **Metallicity [M/H]**    | -0.5, -1.0, -1.5, -2.0, -2.5, -0.0, 0.0, +0.0, +0.2, 0.2, +0.5, 0.5                                                                       |
| **Gravity log(g)**       | 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, +0.0, +0.5, +1.0, +1.5, +2.0, +2.5, +3.0, +3.5, +4.0, +4.5, +5.0                  |
| **Effective Temperature**| 3500 K to 50000 K (with an uneven grid)                                                                                                    |

### Detailed Descriptions

- **Metallicity [M/H]**: 
    - This parameter defines the abundance of metals (all elements heavier than hydrogen and helium) relative to solar metallicity.
    - Available values: **-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, +0.2, +0.5**
    
- **Gravity log(g)**: 
    - The logarithm of surface gravity (g) measured in cm/s².
    - Available values: **0.0 to 5.0** in steps of **0.5**, where lower values represent giant stars and higher values correspond to dwarf stars.
    
- **Effective Temperature (Teff)**: 
    - The surface temperature of a star, with the grid covering temperatures from **3500 K** to **50000 K**. The grid is uneven, with more models at certain temperature ranges.

### Usage
The ATLAS9 model grid provides a flexible framework for studying stellar atmospheres across a wide range of conditions, making it suitable for different astrophysical research and applications, including the study of stellar evolution, stellar populations, and synthetic photometry.

