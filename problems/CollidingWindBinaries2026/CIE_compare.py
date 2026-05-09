#import numpy as np
import NebulaPy.src as nebula
from NebulaPy.tools import util
import time
#from pypion.ReadData import ReadData
import matplotlib.pyplot as plt
import numpy as np


# info: code to test emission measure calculation.
import ChiantiPy.core as ch


# constants
cm2au = 6.68459e-14  # cm to au conversion factor

# Macbook
OutputDir = '/Users/tony/Desktop/CWBs-NEMOv1/Post-Processing/XraySpectrum'  # Output image directory


# fixing the following issues
# info: CIE calculation match for Nebulapy.
database = nebula.database(verbose=True)
database.load_cie()
Temperature = np.logspace(np.log10(1.e6), np.log10(6.e6), 100)
plt.figure()
plt.title("Fe Ionisation Fractions (CIE)")
# Loop over all Fe ions (Fe I to Fe XXVI → 1 to 26)
stages=[13, 14, 15]
for i in range(stages[0], stages[-1] + 1):
    chianti_ion = f'fe_{i}'

    CIE_data = database.get_cie_fraction(chianti_ion, Temperature)

    plt.plot(np.log10(Temperature), CIE_data,
             linewidth=2,
             label=chianti_ion.upper())
plt.ylim(1.e-2, 0.4)  # Adjust y-axis limits for better visibility
plt.xlabel("log10(Temperature [K])")
plt.ylabel("Ion Fraction")
plt.legend(ncol=2, fontsize=8)  # better layout for many ions
plt.yscale('log')
outfile = OutputDir + "/cie_fe_all.png"
plt.savefig(outfile, dpi=300)
plt.close()


fe = ch.ioneq(26)
fe.load()
plt.figure()
fe.plot(stages=stages, tRange=[1.e+6, 6.e+6], yr=[1.e-2, 0.4])
plt.tight_layout()
plt.show()