import numpy as np
from pathlib import Path
import os
from NebulaPy.tools import util as util

class database:
    import os

    def __init__(self, verbose=False):
        self.verbose = verbose

        # get database
        database = os.environ.get("NEBULAPYDB")

        # Check if the database exists, exit if missing
        if database is None:
            util.nebula_exit_with_error(
                "required database dir missing, install database to proceed"
            )

        self.cie_file = os.path.join(database, "IonBalance", "CIE.txt")



    ######################################################################################
    # Load CIE ion fraction table from the database directory
    ######################################################################################
    def load_cie(self):

        cie_file = self.cie_file

        if not os.path.exists(cie_file):
            util.nebula_exit_with_error(f"CIE file not found: {cie_file}")

        header = None
        data = []

        if self.verbose:
            print("[CIEIonFractions] Reading CIE ion fraction table")
            print(f"  File: {cie_file}")

        with open(cie_file, "r") as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith("#"):
                    continue

                parts = line.split()

                if header is None:
                    header = parts
                    continue

                data.append([float(x) for x in parts])

        if header is None:
            util.nebula_exit_with_error("CIE file header not found.")

        if not data:
            util.nebula_exit_with_error("CIE file contains no data.")

        self.data = np.array(data, dtype=float)
        self.col_index = {name: i for i, name in enumerate(header)}

        if self.verbose:
            print(f"  Loaded grid points : {self.data.shape[0]}")
            print(f"  Loaded ion columns : {self.data.shape[1] - 1}")
            print("  Status            : ready")

    ######################################################################################
    # Interpolate ion fraction for a given ion and temperature(s)
    ######################################################################################
    def get_ion_fraction(self, ion, Temperature):
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

        if self.verbose:
            print(f"[CIEIonFractions] Interpolating")
            print(f"  Ion          : {ion}")
            print(f"  N(Temperature): {Temperature.size}")

        logT = np.clip(logT, Tgrid[0], Tgrid[-1])
        frac = np.interp(logT, Tgrid, fgrid)

        return frac[0] if scalar else frac
