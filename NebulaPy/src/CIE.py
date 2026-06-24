import numpy as np
from pathlib import Path
import os
from NebulaPy.src import Utils as utils
from tqdm import tqdm

class cieMode:

    def __init__(self, verbose=False):
        self.verbose = verbose

        # get database
        database = os.environ.get("NEBULAPYDB")

        # Check if the database exists, exit if missing
        if database is None:
            utils.nebula_exit_with_error(
                "required database dir missing, install database to proceed"
            )

        self.cie_file = os.path.join(database, "IonBalance", "CIE.txt")

    ######################################################################################
    # Load CIE ion fraction table from the database directory
    ######################################################################################
    def loadCIEFile(self):

        cie_file = self.cie_file

        if not os.path.exists(cie_file):
            utils.nebula_exit_with_error(f"CIE file not found: {cie_file}")

        header = None
        data = []

        if self.verbose:
            print(" [ CIE ]: Loading collision ionization equilibrium table")

        with open(cie_file, "r") as f:

            lines = f.readlines()

            iterator = tqdm(
                lines,
                desc=" Importing CIE grid",
                unit=" lines",
                ncols=90,
                disable=not self.verbose
            )

            for line in iterator:

                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                parts = line.split()

                if header is None:
                    header = parts
                    continue

                data.append([float(x) for x in parts])

        if header is None:
            utils.nebula_exit_with_error("CIE file header not found.")

        if not data:
            utils.nebula_exit_with_error("CIE file contains no data.")

        self.data = np.array(data, dtype=np.float64)

        # Keep the first column as log_T and convert ion columns to PION symbols
        self.col_index = {
            ("log_T" if i == 0 else utils.getPionSymbol(name)): i
            for i, name in enumerate(header)
        }

        self.AllSpecies = np.array(
            [utils.getPionSymbol(name) for name in header[1:]],
            dtype=str
        )

        if self.verbose:
            print(" CIE table summary")
            print(f" Temperature grid points : {self.data.shape[0]}")
            print(f" Ion fraction columns    : {self.data.shape[1] - 1}")

    ######################################################################################
    # Interpolate ion fraction for a given ion and temperature(s)
    ######################################################################################
    def getCIEFraction(self, ion, Temperature):
        """
        Interpolate ion fraction for scalar or array Temperature (Kelvin).
        """
        if ion not in self.col_index:
            raise KeyError(f"Ion '{ion}' not found.")

        Temperature = np.asarray(Temperature, dtype=float)
        scalar = Temperature.ndim == 0
        Temperature = np.atleast_1d(Temperature)

        logT = np.log10(Temperature)

        Tgrid = self.data[:, self.col_index["log_T"]]
        fgrid = self.data[:, self.col_index[ion]]

        logT = np.clip(logT, Tgrid[0], Tgrid[-1])
        frac = np.interp(logT, Tgrid, fgrid)

        return frac[0] if scalar else frac

    ######################################################################################
    # Interpolate ion fraction for a given ion and temperature(s)
    ######################################################################################
    def buildCIENumberDensities(self, ElementMassFraction, Temperature):
        """
        Build a dictionary of ion number densities for all ions in the CIE table.
        ElementMassFraction: dict of element mass fractions (e.g. {'H': 0.7, 'He': 0.28, ...})
        Temperature: scalar or array of temperatures (Kelvin)
        Returns: dict of {ion: number density array}
        """

        pass