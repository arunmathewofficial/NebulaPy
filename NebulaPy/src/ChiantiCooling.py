from .Chianti import chianti

"""
===================================
PION-Chianti - Cooling Rate Tables
===================================
Python code to calculate cooling table for Mpv10 using ChiantiPy.
CHIANTI provides a database of atomic data that can be used to interpret
the emission of spectral lines and continua emitted from
high-temperature, optically-thin astrophysical sources.
see: https://chiantipy.readthedocs.io/en/latest/api/ChiantiPy.core.html?highlight=FreeFree#id11
"""
__author__ = "Arun Mathew"

"""
Setting up the environment variable in `.bashrc` file.
------------------------------------------------------
mathew@gamma2021:~$ touch file.bashrc
mathew@gamma2021:~$ chmod +x file.bashrc
mathew@gamma2021:~$ XUVTOP="/home/mathew/Chianti_database"
mathew@gamma2021:~$ export XUVTOP
mathew@gamma2021:~$ bash
mathew@gamma2021:~$ echo $XUVTOP
/home/mathew/Chianti_database/
"""

# ChiantiPy packages
import ChiantiPy.core as ch
import ChiantiPy.tools as util
import ChiantiPy.tools.util as util
import ChiantiPy.tools.data as chdata
from ChiantiPy.base import specTrails

# convert mpv10 package
from chianti2pion import *

# Other packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
import sys
import argparse
import os

print('*************** PION-CHIANTI: Cooling Rate Tables *****************')

# Edit Here
ion_list = ['fe_15']

#############################################################
# Making Output directory
parser = argparse.ArgumentParser(
    description='Process some values.',
    usage='script.py <output-path> ')

parser.add_argument('out_path', help='give output path')
args = parser.parse_args()

out_path = args.out_path
if not out_path.endswith('/'):
    out_path += "/"
# make the directory
outputDir = out_path + 'pion-chianti/'
try:
    # Create target Directory
    os.mkdir(outputDir)
    print("Output directory:", outputDir, "created.")
# If the directory already exist, the pass.
except FileExistsError:
    print("Output directory", outputDir, "already exists.")
    pass


# End of make directory ########################################


# todo_keys ####################################################
def todo_keys(ion_name, temperature):
    '''
     Function to obtain the todo keys for the ion.
    if neutral, no free-free contribution, the ket exclude 'ff'
    if the ion contribute free-bound emission, then the
    key include 'fb'.
    if the ion contributes line emissions, key include
    keywords 'line'

    :param ion_name:
    :param temperature:
    :return: keys
    '''

    ion_info = util.convertName(ion_name)
    ionstage = ion_info['Ion']
    charge = 1 - ionstage

    ion_list = []
    ion_list.append(ion_name)

    AbundanceName = 'unity'
    abundAll = chdata.Abundance[AbundanceName]['abundance']

    species = specTrails()  # species of the object of the class specTrails
    species.AbundAll = abundAll
    species.Temperature = temperature

    species.ionGate(ionList=ion_list, minAbund=None, doLines=True,
                    doContinuum=True, doWvlTest=0, doIoneqTest=0, verbose=False)

    keys = species.Todo[ion_list[0]]
    if charge == 0: keys = '_line_'
    return keys


# End of todo_keys ####################################################


# Cooling Table Function ###############################################
def cooling_table(ion_name, ion_info, keys, temperature, eDensity):
    '''
    Calculate the radiative loss rate as a function of temperature and density.
    '''

    # object for the class continuum
    if 'ff' in keys or 'fb' in keys:
        cont = ch.continuum(ion_name, temperature)

    if 'line' in keys:
        thisIon = ch.ion(ion_name, temperature, eDensity, abundance='unity')

    # set rate tables to zero
    total_rate = np.zeros_like(temperature)
    freefree_rate = np.zeros_like(temperature)
    freebound_rate = np.zeros_like(temperature)
    boundbound_rate = np.zeros_like(temperature)
    twophoton_rate = np.zeros_like(temperature)

    # get relevant ion information
    Z = ion_info['Z']
    ionstage = ion_info['Ion']
    dielectronic = ion_info['Dielectronic']
    experimental_name = ion_info['experimental']

    # Free-free loss rate
    # -------------------
    '''
    In this section we calculate the free-free energy loss 
    rate of an ion using the free-free radiative loss rate
    equation given by Eq. 5.15a of Ref[Rybicki and Lightman,
    1979,  Radiative Processes in Astrophysics]

    One may skip the neutral ions since the expression of 
    the rate goes like Z^2. 

    The calculated rates are returned to the FreeFreeLoss
    attributes with title temperature, rate, gf, prefactor
    '''
    if 'ff' in keys:
        cont.freeFreeLoss(includeAbund=False, includeIoneq=False)
        freefree_rate = cont.FreeFreeLoss['rate']
        total_rate += freefree_rate
    # End of Free-free loss rate

    # Free-bound loss rate
    # --------------------
    '''
    Calculate free-bound loss rate using Eq.1a of 106 is 
    integrated over wavelength to get the free-bound loss rate,
    where the free-bound Gaunt factor is given by Eq. 15 of 
    106 and is the numerical constant C_ff from Eq. 4 of 108.

    The rate is returned to FreeBoundLoss with attribute 'rate'.
    '''
    if 'fb' in keys:
        cont.freeBoundLoss(includeAbund=False, includeIoneq=False)
        freebound_rate = cont.FreeBoundLoss['rate']
        total_rate += freebound_rate
    # End of free-bound loss rate

    # line emission rate
    # ------------------
    '''
    In the class ion, we reset the ionization equilibrium array as
    a function of temperature as all ones.
    added a new line at 2863
    '''
    if 'line' in keys:
        # thisIon = ch.ion(ion_name, temperature, eDensity, abundance='unity')
        thisIon.intensity(allLines=0)
        thisIon.boundBoundLoss()
        boundbound_rate = thisIon.BoundBoundLoss['rate']
        total_rate += boundbound_rate

        if (Z - ionstage) in [0, 1] and not dielectronic:
            thisIon.twoPhotonLoss()
            twophoton_rate = thisIon.TwoPhotonLoss['rate']
            total_rate += twophoton_rate
    # End of line emission rate

    radiative_loss = {'temperature': temperature, 'density': eDensity,
                      'ff-rate': freefree_rate, 'fb-rate': freebound_rate,
                      'bb-rate': boundbound_rate,
                      'twophoton-rate': twophoton_rate,
                      'total-rate': total_rate,
                      'experimental_name': experimental_name,
                      'Z': Z, 'Ion_stage': ionstage,
                      'Dielectronic': dielectronic
                      }
    del freefree_rate, freebound_rate, \
        boundbound_rate, twophoton_rate, total_rate
    return radiative_loss


# end of cooling table function ###################################


# plot cooling table ##############################################
def plot_cooling_table(radiative_loss, keys, density_string, out_path):
    '''
    plot cooling table for a specific ion with particular
    electron number density

    :param radiative_loss:
    :param keys:
    :param density_string:
    :return:
    '''

    temperature = radiative_loss['temperature']
    ff_rate = radiative_loss['ff-rate']
    fb_rate = radiative_loss['fb-rate']
    bb_rate = radiative_loss['bb-rate']
    twophoton_rate = radiative_loss['twophoton-rate']
    total_rate = radiative_loss['total-rate']
    experimental_name = radiative_loss['experimental_name']
    Z = radiative_loss['Z']
    ionisation_stage = radiative_loss['Ion_stage']
    dielectronic = radiative_loss['Dielectronic']

    # print(f'Plotting {experimental_name} cooling rate for log(ne) '
    #      f'cm^-3 = {density_string}')
    warnings.filterwarnings("ignore")
    plt.figure()
    plt.title('Species: ' + exp_name + ',  log(n_e): ' + density_string)
    plt.xlabel('log(T) K')
    plt.ylabel('log($\Lambda$) erg cm^3 s^-1')
    plt.xlim([1, 9])
    plt.ylim([-30, -15])
    if 'ff' in keys:
        plt.plot(np.log10(temperature), np.log10(ff_rate), label='free-free')
    if 'fb' in keys:
        plt.plot(np.log10(temperature), np.log10(fb_rate), label='free-bound')
    if 'line' in keys:
        plt.plot(np.log10(temperature), np.log10(bb_rate), label='bound-bound')
        if (Z - ionisation_stage) in [0, 1] and not dielectronic:
            plt.plot(np.log10(temperature), np.log10(twophoton_rate), label='two-photon')

    plt.plot(np.log10(temperature), np.log10(total_rate), label='total')
    plt.legend()
    plt.savefig(out_path + '/' + experimental_name + '_' + density_string + '.png')
    plt.close()


# end of plot cooling table ##############################################


def print_to_file(data, output_path, filename):
    print(f'Writing to file {filename}.txt ...')
    table = pd.DataFrame(data)
    # table = table.to_string(index=False)
    np.savetxt(output_path + filename + '.txt', table.values, fmt='%f')
    return table.values


# Main ############################################################

exp_names = []
for ion in ion_list:
    exp_names.append(util.convertName(ion)['experimental'])

if len(ion_list) == 1:
    exp_name = util.convertName(ion_list[0])['experimental']
    print('Species: ', exp_names[0])
else:
    print('Species List: ', end=' ')
    print(exp_names)

# Making temperature and electron number density array
# identical to the one in Mellema cooling table
temperature_array = []
for i in range(81):
    log_temp = 1.0 + i * 0.1
    temperature_array.append(pow(10, log_temp))

ne = []
for i in range(13):
    log_ne = i * 0.5
    ne.append(pow(10, log_ne))

for ion_name in ion_list:
    # getting ion info
    ion_info = util.convertName(ion_name)
    exp_name = ion_info['experimental']

    # creating list for this ion
    data = {}
    # print(f'Writing temperature to {exp_name} table')
    data.update({'Temperature': np.log10(temperature_array)})

    # Obtain the todo_keys
    ion_keys = todo_keys(ion_name, temperature_array)

    # loop over electron number density
    for i in range(len(ne)):
        print_statement = f"Computing {exp_name} with log(ne) {np.log10(ne[i])}"
        print(print_statement, end="\r")
        radiative_loss = cooling_table(ion_name, ion_info, ion_keys, temperature_array, ne[i])
        ne_string = str(np.log10(ne[i]))
        # plot_cooling_table(radiative_loss, ion_keys, ne_string, outputDir)

        # making data to print to file
        total_rate = radiative_loss['total-rate']
        filename = radiative_loss['experimental_name']
        total_rate[total_rate < pow(10, -50)] = pow(10, -50)
        data.update({ne_string: np.log10(total_rate)})

    sys.stdout.write('\x1b[2K')
    sys.stdout.write(f"Computing {exp_name}: Done\n")
    table = print_to_file(data, outputDir, filename)
    mpv10_format(table, outputDir, filename)






