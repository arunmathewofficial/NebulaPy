import os
import glob
import re
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import NebulaPy.tools.constants as const

nebula_database = '/Users/tony/Desktop/NebulaPy/Nebula_Database'


class sed:

    def __init__(self, energy_bins, plot=None, pion=None, verbose=False):
        self.EnergyBins = energy_bins
        self.Plot = plot
        self.Pion = pion
        self.Verbose = verbose
        self.container = {'energy_bins': self.EnergyBins,
                          'plot': self.Plot, 'pion': self.Pion}
        self.AtlasDatabase = nebula_database + '/SED/Atlas/'
        self.PoWRDatabase = nebula_database + '/SED/PoWR/'
        self.setup_lambda_bin()

    ##############################################################################
    def setup_lambda_bin(self):
        '''
        making wavelenght bins from the given energy bins
        :return:
        '''
        # making wave-lenght bins from the given energy bins
        given_lambdabins = []
        # Loop through the rows (outer loop)
        for Ebin in self.EnergyBins:
            # Loop through the columns (inner loop)
            lambdaBin = []
            for energy in Ebin:
                Lambda = const.ev2Ang / energy  # Waveleghts are given in Angstrom
                lambdaBin.append(Lambda)
            lambdaBin.reverse()
            given_lambdabins.append(lambdaBin)
        given_lambdabins.reverse()
        self.given_lambdabins = given_lambdabins
        self.container['lambda_bins'] = given_lambdabins
        self.container['wvl_unit'] = 'Angstrom'

    #########################################################################################
    # progress bar
    #########################################################################################
    def progress_bar(self, iteration, Nmodels, prefix='', fill='█'):
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total_files - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            fill       - Optional  : bar fill character (Str)
        """
        length = 20  # length of progress bar
        percent = ("{0:." + str(1) + "f}").format(100 * (iteration / float(Nmodels)))
        filled_length = int(length * iteration // Nmodels)
        bar = fill * filled_length + '-' * (length - filled_length)

        sys.stdout.write(f'\r {prefix}: |{bar}| {percent}% complete')
        sys.stdout.flush()
        # Print New Line on Complete
        if iteration == Nmodels:
            sys.stdout.write('\r' + ' ' * (len(f'{prefix}: |{bar}| {percent}% complete')) + '\r')
            sys.stdout.flush()
            sys.stdout.write(f'\r {prefix}: completed\n')



    ##############################################################################
    def bundle_up_atlas_models(self, metallicity):
        '''
        TODO: no completed
        1.get the model and return it
        2.if the parameter are wrong return some advice

        :param metallicity:
        :param gravity:
        :return:
        '''

        sign = None
        if metallicity < 0:
            sign = 'm'
        if metallicity >= 0:
            sign = 'p'

        metallicity_str = str(metallicity).replace('.', '')

        # Construct model directory name
        model_dirname = 'ck' + sign + metallicity_str

        model_dir = os.path.join(self.AtlasDatabase, model_dirname)

        # Define the fits file pattern.
        model_file_pattern = 'ck' + sign + metallicity_str + '_*.fits'

        # Collecting all fits files with similar names.
        model_set = glob.glob(os.path.join(model_dir, model_file_pattern))

        self.model_name = 'ck' + sign + metallicity_str

        return model_set







    ######################################################################################
    def sort_modelset_by_temperature(self, model_set):
        '''
        This will extract temperature values from the all the file name
         in the model set
        :param model_set:
        :return:
        '''
        Teff = []
        TeffFileSet = []
        for file in model_set:
            FileMatch = re.search(self.model_name + r'_(\d+)\.fits', file)
            temperature = FileMatch.group(1)
            if FileMatch:
                Teff.append(float(temperature))
                TeffFileSet.append((float(temperature), file))
        # Sorting.
        TeffFileSet.sort(key=lambda x: x[0])
        # Sorted finalized file set.
        model_set = [file for _, file in TeffFileSet]
        del TeffFileSet
        self.Nmodels = len(model_set)
        self.Teff = sorted(Teff)
        self.container['Teff'] = self.Teff
        return model_set

    ##############################################################################
    def CastelliKuruczAtlas(self, metallicity, gravity):
        '''

        :param metallicity:
        :param gravity:
        :return:
        '''
        self.container['model'] = 'Atlas'
        self.Model = 'Atlas'
        self.container['metallicity'] = metallicity
        self.metallicity = metallicity
        self.container['gravity'] = gravity
        self.gravity = gravity

        # model set corresponding to a specific metallicity
        model_set = self.bundle_up_atlas_models(metallicity)
        # model set sorting according to accending temperature of SED
        model_set = self.sort_modelset_by_temperature(model_set)


        # make column name from the input gravity parameter
        # TODO: do more, cross reference with existing info
        column_name = 'g' + str(gravity).replace('.', '')

        self.sed_set_name = 'atlas_' + self.model_name + column_name

        prefix_comment = self.sed_set_name + ' binning'
        model_flux_set = []
        model_lambda_set = []
        binned_flux_set = []
        total_flux_set = []


        for model_index, model in enumerate(model_set, start=1):

            self.progress_bar(model_index, self.Nmodels,
                              prefix=prefix_comment)

            # Fits File Reading #####################################################
            # Atlas model in PyMicroPION database are in fits file format.
            # The datas are obtained form: https://www.stsci.edu/hst/instrumentation/
            # reference-data-for-calibration-and-tools/astronomical-catalogs/castelli
            # -and-kurucz-atlas
            # Open the FITS file to read wavelength and flam
            with fits.open(model) as hdul:
                # The open function returns an object called an HDUList which is a
                # list-like collection of HDU objects. An HDU (Header Data Unit) is
                # the highest level component of the FITS file structure, consisting
                # of a header and (typically) a data array or table.

                # Files opened for reading can be verified and fixed with method
                # HDUList.verify. This method is typically invoked after opening the
                # file but before accessing any headers or data:
                hdul.verify('fix')

                # hdul[0] is the primary HDU, hdul[1] is the first extension HDU
                # The header in hdul[0] of the Atlas fits file contains comprehensive
                # data details, while the actual data is stored in hdul[1].

                # Retrieve the wavelength and FLAM from the designated FITS file.
                model_lambda = hdul[1].data['Wavelength']
                model_flux = hdul[1].data[column_name]

                # Binning the wavelength and flux of the original model.
                # Binning lambda values and flam values according to the Lambda
                # Bins (obtained from the input energy bin).
                binned_lambda = []
                binned_flux = []
                for bin in self.given_lambdabins:
                    sub_binned_lambda = []
                    sub_binned_flux = []
                    for lam, flux in zip(model_lambda, model_flux):
                        if bin[0] <= lam <= bin[1]:
                            sub_binned_lambda.append(lam)
                            sub_binned_flux.append(flux)
                    binned_lambda.append(sub_binned_lambda)
                    binned_flux.append(sub_binned_flux)
                    del sub_binned_lambda
                    del sub_binned_flux

                # calculating the normalization factor, perform integration across
                # the entire wavelength domain to obtain the total flux.
                total_flux = np.trapz(np.asarray(model_flux), np.asarray(model_lambda))
                # Append the total Flux into TotalFlux_BundledGrids
                total_flux_set.append(total_flux)

                model_lambda_set.append(model_lambda)
                model_flux_set.append(model_flux)

                # calculating flux in each binned flux by integrating within the bin
                # interval
                flux_bin = []
                for i in range(len(binned_lambda)):
                    flux_bin.append(np.trapz(np.asarray(binned_flux[i]),
                                             np.asarray(binned_lambda[i])))

                # Reverse the order of the flux bins since we are interested in
                # obtaining flux in energy bins.
                flux_bin.reverse()
                # Normalizing flux within each bin.
                if total_flux == 0.0:
                    norm_flux_bin = np.zeros_like(flux_bin)
                else:
                    norm_flux_bin = flux_bin / total_flux
                # Removing binned flux array
                del flux_bin

                # Appending the result fractional flux for each sed effective temperature
                # to final binned flux set
                binned_flux_set.append(norm_flux_bin)
                # BinnedFracFluxSet.append([f"{x:.4e}" for x in NormBinFlux])
            # End of Fits File Reading ##################################################

        self.binned_flux_set = binned_flux_set
        self.container['binned_flux'] = binned_flux_set
        self.container['total_flux'] = total_flux_set
        self.model_lambda_set = model_lambda_set
        self.model_flux_set = model_flux_set
        # End of binning model set ######################################################


        #if self.Plot is not None:
        #    self.plotter(self.Plot, binned_flux_set)

        # Gathering SED model info
        model_info = f"{'atlas_' + self.sed_set_name + '.info'} = " \
                     f"\"#COMMENT: SED -' + {self.Model} 'Model Atmospheres\\n\""\
                     f"\"#COMMENT: SED parameters:\\n\""\
                     f"\"#COMMENT: LOG Z = {self.metallicity}\\n\""\
                     f"\"#COMMENT: LOG G = {self.gravity}\\n\";"
        if self.Pion is not None:
            self.pion_format(self.Pion, binned_flux_set, model_info)


    ######################################################################################
    # plotter
    ######################################################################################
    def plotter(self, plot_path, binned_flux_set):
        '''

        :param PlotDir:
        :return: None
        '''

        energy_bins = self.EnergyBins
        min_plot_energy = 5.0
        max_plot_energy = 80.0

        plot_dir = os.path.join(plot_path, self.sed_set_name)
        os.makedirs(plot_dir, exist_ok=True)
        if not plot_dir.endswith('/'):
            plot_dir += '/'

        # Calculate bin center value and bin width for each energy bin
        bin_widths = [bin_max - bin_min for bin_min, bin_max in energy_bins]
        bin_centers = [(bin_min + bin_max) / 2 for bin_min, bin_max in energy_bins]

        prefix_comment = self.sed_set_name + ' plotting'

        for model_index, Teff in enumerate(self.Teff, start=0):

            self.progress_bar(model_index, self.Nmodels, prefix=prefix_comment)

            fig, axs = plt.subplots(2, 1, figsize=(12, 6))
            # converting model lambda to electron volt unit
            model_energy = [const.ev2Ang / Lambda for Lambda in self.model_lambda_set[model_index]]

            # SubPlot 1: Plot the original model flux data
            axs[0].plot(model_energy, self.model_flux_set[model_index], label="Original Flux",
                        color='black', linestyle='-', linewidth=2)
            axs[0].set_xlabel("Energy, eV")
            axs[0].set_ylabel(r'$\rm log \ F_{\lambda} \  (ergs \, cm^{-2} s^{-1} \AA^{-1})$')
            axs[0].set_yscale('log')
            if self.Model == 'Atlas':
                axs[0].set_title(f'Model: {self.Model } {self.model_name},  T_eff: {Teff} K,'
                                 f'  Log Z: {self.metallicity}, Log g: {self.gravity}', fontsize=12)
            if self.Model == 'PoWR':
                axs[0].set_title(f'Model: {self.model_name},  T_eff: {Teff} K,  Log Z: {self.metallicity}, '
                                 f' Log g: {self.gravity}', fontsize=12)

            axs[0].tick_params(axis="both", direction="inout", which="both",
                               bottom=True, top=True, left=True, right=True, length=3)
            axs[0].legend(loc='upper right')
            axs[0].set_xlim(min_plot_energy, max_plot_energy)
            axs[0].grid(True, linestyle='--', alpha=0.5)

            # SubPlot 2: Plot (bar plot) the binned data calculated by PyMicroPion
            axs[1].bar(bin_centers, binned_flux_set[model_index], width=bin_widths,
                       align='center', color='orange',
                       alpha=0.5, label="NebulaPy")
            axs[1].set_xlabel("Energy, eV")
            axs[1].set_ylabel("log Fractional Binned Flux")
            axs[1].set_yscale('log')
            axs[1].tick_params(axis="both", direction="inout", which="both", bottom=True, top=True,
                               left=True, right=True, length=3)
            axs[1].legend(loc='upper right')
            axs[1].set_xlim(min_plot_energy, max_plot_energy)
            axs[1].grid(True, linestyle='--', alpha=0.5)

            plt.tight_layout()
            image_file = plot_dir + self.sed_set_name + '_' + str(Teff) + ".png"
            plt.savefig(image_file, dpi=100)
            plt.close(fig)

    ######################################################################################
    # pion format
    ######################################################################################

    def pion_format(self, pion_format_path, binned_flux_set, model_info):

        if not pion_format_path.endswith('/'):
            pion_format_path += '/'
        pion_format_filename = self.sed_set_name + '.txt'
        pion_format_file = os.path.join(pion_format_path, pion_format_filename)
        print(pion_format_file)

        '''
        # **************************************************************************
        # Converting Input Energy Bins PION cpp 2D-Vector Format
        EnergyBins = self.SEDParameter['EnergyBins']
        #EnergyBins = EnergyBins.values.tolist()
        EnergyBins_PIONFormatName = 'energy_bins'

        EnergyBins_num_rows = len(EnergyBins)
        EnergyBins_num_cols = len(EnergyBins[0]) if EnergyBins_num_rows > 0 else 0

        # Convert the Python 2D list to a C++ 2D array string representation
        EnergyBins_PIONFormat = f"{EnergyBins_PIONFormatName}= {{\n"
        for row in EnergyBins:
            EnergyBins_PIONFormat += "    {"
            EnergyBins_PIONFormat += ', '.join(map(lambda x: f"{x:.6e}", row))
            EnergyBins_PIONFormat += "},\n"
        EnergyBins_PIONFormat += "};"
        '''

        with open(pion_format_file, 'w') as outfile:
            outfile.write(model_info)

        '''
        # Converting Binned SED Data Frame into PION cpp 2D-Vector Format
        SED = self.DataSet_DF.drop(['MODEL', 'TOTFLUX'], axis=1)
        SED = SED.values.tolist()
        SED_PIONFormatName = 'atlas_' + self.GridName + '.data'

        SED_num_rows = len(SED)
        SED_num_cols = len(SED[0]) if SED_num_rows > 0 else 0

        SED_PIONFormat = f"{SED_PIONFormatName} = {{\n"
        for row in SED:
            SED_PIONFormat += "    {"
            SED_PIONFormat += ', '.join(map(lambda x: f"{x:.6e}", row))
            SED_PIONFormat += "},\n"
        SED_PIONFormat += "};"
        # **************************************************************************

        # PIONFormat text file content
        PIONFormat = ModelInfoText + '\n\n' + EnergyBins_PIONFormat + '\n\n' + SED_PIONFormat

        #f'\n\n# PION-Formating ***\n' \
        #f'# PION-Format File: {PIONFormatFile}\n' \
        #f'# File Content: Model Configurations, Energy Bins, Binned SED\n' \
        #f'# Energy Bins - 2D Array with {EnergyBins_num_rows} bins (in unit of eV)\n' \
        #f'# Binned Spectral Energy Distribution - 2D Array with {SED_num_rows} Atlas Models\n' \
        #f'# The initial element in each model represent its effective temperature in K.\n' \
        #f'# The subsequent elements represent fractional dimensionless flux in each bin.\n'

        PIONFormatFile = PlotFormatPathDir + 'atlas_' + self.GridName + '_pionformat.txt'
        logger.info(f"Writing PION-formatted PoWR Data into: " + ' ' + PIONFormatFile)



        with open(PIONFormatFile, 'w') as outfile:
            #outfile.write()
            outfile.write(PIONFormat)
        '''
    ################################################################################





